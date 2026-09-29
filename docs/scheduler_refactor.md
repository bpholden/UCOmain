# Scheduler refactor: turning UCOScheduler into a class

Status: steps 1-7 done and committed on `scheduler_object`. Step 8 (bug fixes)
is partly done: four bugs are fixed and five remain (see "Outstanding bugs").
This document describes the branch as built. Where the build departed from the
original plan, the reason is given.

## Goal

Before this refactor, `Main/UCOScheduler.py` was a 1147-line module of free
functions with a module-level global. Every call to `get_next()` threaded
`star_table`, `rank_table`, `hour_table`, `observed`, `focval`, `owner`,
`outdir` and `start_time` through by hand.

After it:

* `UCOScheduler` is a class. `Main.sin` builds one instance, `Observe` holds it
  and calls `self.scheduler.get_next(...)`.
* Target-table state and the files those tables live in are in one object,
  `UCOTargetTables`, separate from the scheduler's own run-state.
* The pure helpers (scriptobs line formatting, visibility and condition cuts)
  are in their own modules so they can be tested without a scheduler.
* Steps 1-7 were verified to pick the same targets and write the same files as
  the code they replaced (see "How this was verified"). The step 8 bug fixes
  change target selection on purpose.

## Decisions

| Question | Decision |
| --- | --- |
| Class relationship | Composition: `UCOScheduler` holds a `UCOTargetTables` (the old `UCOTargets`, renamed) |
| Tables class name | `UCOTargetTables`, keeping the `UCO` prefix the rest of the package uses (rather than `TargetTables`) |
| Entry point name | `get_next()` (snake_case, matches the repo) |
| File layout | Split into four modules plus a test module |
| Pure selection helpers | `compute_priorities`, `need_cal_star`, `make_result`, `last_attempted` stay module-level functions in `UCOScheduler.py`, not methods; they use no scheduler state and are easier to test that way |
| `start_time` | Lives on the scheduler (`scheduler.start_time`); `Observe` reads and writes it there |
| `do_templates` / `do_too` defaults | `False`, as the old `get_next` had; both can also be overridden per call |
| Failed-object tracking | The list stays on the scheduler (`last_objs_attempted`). Whether it is used is a constructor option, `track_failures`, default `False`. `Main.sin` passes `track_failures=False`, matching `main`, which had turned tracking off |
| Template budget | Lives on the scheduler: `tot_temps` (constructor, default `None` = no limit) and the count `n_temps`. `get_next` counts each template it returns and stops offering templates once `n_temps >= tot_temps`. `Main.sin` passes `do_templates=True, tot_temps=4`, the old `Observe` defaults. A template counts when `get_next` returns it, even if `Observe` then fails to write it to scriptobs (previously such a template did not count) |
| `utils/` | Port the live scripts; leave the already-dead ones alone |

## Why composition and not one merged class

`UCOTargets` was not owned by the scheduler. `Main.sin` built it and handed the
*same* object to two places (before the refactor):

```
Main.sin:488   uco_targets = UCOTargets.UCOTargets(opt)
Main.sin:489   getUCOTargets.getUCOTargets(uco_targets, ...)   # background thread
Main.sin:490   observe = Observe.Observe(apf, tel, opt, uco_targets, ...)
```

`getUCOTargets` is a `threading.Thread` that waits on sun elevation and then calls
`make_hour_table()` / `make_star_table()`, i.e. it downloads the Google sheets and
writes them to disk on its own schedule. Merging the tables into `UCOScheduler`
would mean that download thread mutates the scheduler object directly, and it
would drag the selection algorithm into `getUCOTargets`'s import graph. Keeping
them as two objects preserves the separation: one object owns the tables and
the files, the other owns target selection.

## Lifetime: the scheduler is built in Main.sin, not in Observe

`Main.sin` **recreates** `Observe` when its thread dies (`Main.sin:533`):

```
if observe.is_alive() is False:
    observe = Observe.Observe(apf, tel, opt, scheduler, task=parent)
    observe.start()
```

`last_objs_attempted`, the list of objects that just failed to observe, used to
survive that restart because it was a module global. If the scheduler were
constructed inside `Observe.__init__`, that list would be silently wiped every
time the thread restarted, and the scheduler would immediately re-select a
target it had already failed on. So `Main.sin` constructs the scheduler
alongside the tables and passes it to `Observe` (`Main.sin:488-501`):

```python
uco_targets = UCOTargetTables.UCOTargetTables(opt)
_ = getUCOTargets.getUCOTargets(uco_targets, task=parent, wait_time=target_time)

start_time = None
if opt.start:
    try:
        start_time = float(opt.start)
    except ValueError as e:
        apflog("ValueError: %s" % (e), echo=True, level='error')
scheduler = UCOScheduler.UCOScheduler(uco_targets, owner=opt.owner if opt.owner else 'public',
                                      start_time=start_time, track_failures=False,
                                      do_templates=True, tot_temps=4)
observe = Observe.Observe(apf, tel, opt, scheduler, task=parent)
```

`Observe` reaches the tables through `self.scheduler.targets`, so it takes one
argument instead of two.

**Behavior change:** because `start_time` and the template budget now live on
the scheduler, they also survive an `Observe` restart. Previously
`Observe.__init__` re-read `start_time` from `--start` and reset the template
count to zero on every restart. This was a deliberate choice; the alternative
(keep `start_time` on `Observe` and pass it per call) was considered and
rejected.

## File layout

| File | Contents |
| --- | --- |
| `Main/SchedulerConsts.py` (existing) | gained `ACQUIRE`, `BLANK`, `FIRST`, `LAST`, `BUFFERSEC`, `BUFFER` |
| `Main/UCOTargetTables.py` (renamed from `UCOTargets.py`, 314 lines) | `class UCOTargetTables`: rank/hour/star tables and all disk bookkeeping |
| `Main/ScriptobsLine.py` (new, 216 lines) | pure string generation for scriptobs lines |
| `Main/Observability.py` (new, 256 lines) | pure array/astro filters, no scheduler state |
| `Main/UCOScheduler.py` (rewritten, 565 lines) | the pure selection helpers and `class UCOScheduler` |
| `Main/test_UCOScheduler.py` (new, 154 lines) | the `test_*` functions formerly at the bottom of `UCOScheduler.py` |

`UCOTargets.py` was renamed with `git mv`, so its history carries over.

## Where every symbol went

### `Main/ScriptobsLine.py`

| From `UCOScheduler.py` | Note |
| --- | --- |
| `make_scriptobs_line` | unchanged; `utils/make_scriptobsline.py` and `gen_template_entry.sin` call it |
| `num_template_exp` | unchanged |
| `config_defaults` | unchanged |
| `make_obs_block` | **dead**: its only caller is commented out in `make_result` (`UCOScheduler.py:184`). Moved; see outstanding bug 2, "Disabled obsblock path" |

### `Main/Observability.py`

All pure functions of `(star_table, moon, apf_obs, dt, ...)`, no instance state:

`compute_datetime`, `tot_exp_times`, `time_check`, `condition_cuts`,
`behind_moon`, `template_conditions`, `find_closest`, `find_Bstars`,
`enough_time_templates`.

`compute_preferred_el` had no live callers and was deleted (commit `f4ae698`)
rather than moved.

These are the functions worth having unit tests for, which is the main reason
they were pulled out.

### `Main/UCOTargetTables.py`

Everything that reads or writes the bookkeeping files:

| Method | From |
| --- | --- |
| `make_rank_table` | `UCOTargets` (unchanged) |
| `make_hour_constraints` | `UCOTargets` (unchanged) |
| `make_hour_table` | `UCOTargets` (unchanged) |
| `make_star_table` | `UCOTargets` (unchanged) |
| `append_too_column` | `UCOTargets` (unchanged) |
| `copy_backup` | `UCOTargets` (unchanged) |
| `update_hour_table(observed, dt, outfn='hour_table', outdir=None)` | **moved** from the free function `UCOScheduler.update_hour_table`. It writes `hour_table` to disk, so it is bookkeeping. Now a no-op when `hour_table` is `None` (that check used to be in `get_next`) |
| `update_from_observed(ptime, outfn=None, toofn='too.dat')` | **new**: wraps `ParseUCOSched.update_local_starlist`, sets `self.star_table`, returns the `ObservedLog`. `outfn` defaults to `self.star_table_name`. If the star table file is missing, it rebuilds it with `make_star_table` and re-applies the observations (the observations only reach the table through the file). If the rebuild fails, it keeps the previous in-memory table |
| `gen_stars()` | **new**: wraps `ParseUCOSched.gen_stars` and keeps the result on `self.stars`. It regenerates on every call; there is no caching logic |
| `check_files(outfn=None)` | **moved** from `Observe.check_files`: restores `googledex.dat` from its `.1` backup if it is missing |

Files this object owns, unchanged in name and format:

* `rank_table` (ascii)
* `hour_table` (ascii)
* `googledex.dat` + `googledex.dat.1` backup (ecsv)
* `too.dat` (ecsv)
* `observed_targets` (read via `ObservedLog`; written by `Observe`)
* the time-left CSV named by `opt.time_left` (read only)

Nothing about these files changed. That is the constraint that makes this
refactor safe to deploy: the on-disk contract is identical, so a night can be
resumed from files written by the old code.

### `Main/UCOScheduler.py`

Module-level functions, kept as pure helpers:

* `need_cal_star(star_table, observed, priorities)`
* `compute_priorities(star_table, cur_dt, observed=None, hour_table=None, rank_table=None, do_templates=False)`
* `make_result(stars, star_table, totexptimes, final_priorities, dt, idx, focval=0, bstar=False, mode='')`
* `last_attempted()`: reads the last scriptobs result from `ktl`

The class:

```python
class UCOScheduler(object):

    def __init__(self, targets, owner='public', outdir=None,
                 do_templates=False, do_too=False, start_time=None,
                 outfn='googledex.dat', toofn='too.dat', track_failures=False,
                 tot_temps=None):
        self.targets   = targets          # UCOTargetTables
        self.owner     = owner
        self.outdir    = outdir or os.getcwd()
        self.outfn     = outfn
        self.toofn     = toofn
        self.do_templates = do_templates
        self.do_too       = do_too
        self.start_time   = start_time
        self.track_failures = track_failures
        self.tot_temps = tot_temps

        # run-state, was a module global
        self.last_objs_attempted = []
        # templates returned so far, counted against tot_temps
        self.n_temps = 0

        # per-call scratch, kept for logging and for Observe to inspect
        self.observed  = None             # ObservedLog from the last refresh
        self.apf_obs   = None
        self.moon      = None
        self.stars     = None
        self.result    = None
        self.template_conditions_met = False

    # --- public -------------------------------------------------------
    def get_next(self, ctime, seeing, slowdown, bstar=False, focval=0,
                 do_templates=None, do_too=None)
    def zero_last_objs_attempted(self)
    def record_last_attempt(self)         # calls last_attempted(), appends failures

    # --- pipeline stages, called in order by get_next ------------------
    def _previous_obs_time(self, dt)      # apfguide midptfin, falls back to dt
    def _refresh_tables(self, ptime)      # update_from_observed, hour table, star table, gen_stars
    def _sky_state(self, dt)              # apf_obs + moon
    def _available(self, dt, seeing, slowdown, bstar, bstars, totexptimes,
                   do_too, do_templates)  # -> (available, cur_elevations, scaled_elevations) or None
    def _select(self, available, final_priorities, bstar,
                cur_elevations, scaled_elevations)   # -> idx or None
    def _add_template(self, res, idx, dt, bstars)
```

`get_next` is about 80 lines, down from 240. It runs the stages above, plus
`compute_priorities` and `make_result`, in order. The early-exit checks and the
template decision stay inline in `get_next` rather than being split further.
Every stage keeps its original `apflog` calls, in the original order, so the
night log reads the same.

Constructor arguments (`owner`, `outdir`, `outfn`, `toofn`, `start_time`,
`do_templates`, `do_too`) used to be per-call `get_next` keyword arguments.
They moved to `__init__` because `Observe` passed the same value on every call.
`do_templates` and `do_too` are also accepted per call as overrides, because
`Observe` changes `self.do_temp` and `self.do_too` during a night. `outdir` is
stored but, as before the refactor, not used by `get_next`.

**Failed-object tracking.** With `track_failures=True`, `get_next` calls
`record_last_attempt()`, which asks `last_attempted()` for the last object
attempted and, if it failed, appends it to `last_objs_attempted`; `_available`
then excludes every object in that list until `zero_last_objs_attempted()` is
called (`Observe` calls it after power-cycling the telescope). With
`track_failures=False`, `last_attempted()` is never called and the list stays
empty.

**Template budget.** At the start of each call `get_next` turns templates off
if `tot_temps` is set and `n_temps >= tot_temps`, before the "Will attempt
templates" log line, so the log reflects the cap. When it returns a template
(`isTemp`), it adds one to `n_temps`. This replaces the counting `Observe` and
`sim_night.py` each did themselves.

Removed from the module: the `last_objs_attempted` global, the module-level
`get_next` and `zero_last_objs_attempted` (a temporary shim during step 6, see
"Steps as done"), the free function `update_hour_table` (moved to
`UCOTargetTables`), and the now-unused `ParseUCOSched` and `UCOTargets`
imports.

### `Main/test_UCOScheduler.py`

`test_basic_ops`, `test_failure`, `test_templates` and `test_main` moved here
unchanged apart from:

* `test_main` reads `time_left.csv` from the current directory instead of
  `/home/holden/time_left.csv`, so it uses the fixture copy.
* `test_basic_ops` and `test_failure` take a `UCOScheduler`, which `test_main`
  builds once. The "nonsensical start time" case sets `scheduler.start_time = 1`
  instead of passing `start_time=1`.

They stay runnable as a script; converting them to pytest is a separate job.
Run them from a *copy* of the fixtures, because `get_next` rewrites
`googledex.dat` and creates `observed_targets` and `hour_table` in the current
directory:

```
cp /Users/holden/src/test_main/{googledex.dat,rank_table,time_left.csv} <scratch>/ && cd <scratch>
PYTHONPATH=/Users/holden/src/UCOmain/Main python3 /Users/holden/src/UCOmain/Main/test_UCOScheduler.py
```

The tests take their start time from the clock (`time.time()`,
`datetime.now()`), so their output depends on when they run. Two runs can only
be compared if they run at the same time. The tests build the scheduler with
the default `track_failures=False`.

## Observe.py changes

| Line (now) | Before | After |
| --- | --- | --- |
| `Observe.py:21-22` | `import UCOScheduler as ds` / `import UCOTargets` | `import UCOScheduler` / `import UCOTargetTables` (both used only by the `__main__` block) |
| `Observe.py:33` | `def __init__(self, apf, tel, opt, uco_targets, tot_temps=4, task='master')` | `def __init__(self, apf, tel, opt, scheduler, task='master')` |
| `Observe.py:496` | `self.check_files()` | `self.scheduler.targets.check_files()` |
| `Observe.py:498` | `ds.get_next(time.time(), seeing, slowdown, self.uco_targets, bstar=..., do_too=..., owner=..., do_templates=..., focval=..., start_time=...)` | `self.scheduler.get_next(time.time(), seeing, slowdown, bstar=self.obs_B_star, focval=self.focval, do_too=self.do_too)` |
| `Observe.py:624` | `ds.zero_last_objs_attempted()` | `self.scheduler.zero_last_objs_attempted()` |

Also:

* `Observe.check_files()` and the `shutil` import it needed were removed; the
  method now lives on `UCOTargetTables`.
* `self.owner` was removed. Passing it to `get_next` was its only use; the owner
  is now set on the scheduler in `Main.sin`.
* The `start_time` initialisation from `opt.start` moved to `Main.sin`, and all
  14 uses of `self.start_time` became `self.scheduler.start_time`, including
  `should_start_list()` (`Observe.py:284`), which clears it an hour after the
  start time, the `MASTER_WHENSTARTLIST` check and the fixed-list branches.
* The template budget moved to the scheduler: `self.do_temp`, `self.n_temps`,
  `self.tot_temps`, the `tot_temps` constructor argument and the counting after
  each template were removed.
* The `__main__` test block builds a `UCOScheduler` (with `do_templates=True,
  tot_temps=4`) and passes it in.

`Observe.py` also has unrelated changes from `main` (the power-cycle limit and
`do_not_open`), brought in by the merge `afd3357`. They are not part of this
refactor.

## Main.sin changes

* `import UCOScheduler as ds` became `import UCOScheduler`; `import UCOTargets`
  became `import UCOTargetTables`.
* The scheduler is built alongside the tables, with `owner`, the parsed
  `start_time`, `track_failures=False`, `do_templates=True` and `tot_temps=4`, and passed to `Observe` both at
  startup and on thread restart (see "Lifetime" above).

## getUCOTargets.py changes

Import and type name only: `UCOTargets.UCOTargets` became
`UCOTargetTables.UCOTargetTables`. The thread does not touch the scheduler.

## utils/ changes

Live scripts, ported:

* `utils/sim_night.py`, `utils/sim_nights.py`: build one `UCOScheduler` before
  the loop (in `sim_nights`, before the loop over nights, so the failed-object
  list carries across nights as the old global did) and call
  `scheduler.get_next(...)` inside it. `outfn`, `outdir` and `start_time` go to
  the constructor, along with `do_templates=True`. `sim_night.py` passes
  `tot_temps=2` (its old "two per night" rule) and no longer counts templates
  itself; `sim_nights.py` has no cap, as before. Neither passes
  `track_failures`, so tracking is off in the sims; pass `track_failures=True`
  to simulate with it on. They still call
  `ParseUCOSched.gen_stars` directly for their own `stars` list.
* `utils/make_scriptobsline.py`: `ds.make_scriptobs_line` became
  `ScriptobsLine.make_scriptobs_line`.
* `utils/gen_template_entry.sin`: `ds.find_Bstars` became
  `Observability.find_Bstars`, `ds.make_scriptobs_line` became
  `ScriptobsLine.make_scriptobs_line`, `UCOTargets` became `UCOTargetTables`.
* `utils/download_googledex.py`: `UCOTargets` became `UCOTargetTables`. It
  still has an unused `import UCOScheduler as ds`, which is harmless.

Left alone, already broken against the deployed code:

* `utils/calc_cadence.py` calls `ds.parseGoogledex` and `ds.DS_APFPRI`
* `utils/calc_etime_precision.py` calls `ds.get_speadsheet`, `ds.getI`, `ds.DS_BV`
* `utils/calc_precision_mag.py` calls `ds.getI`, `ds.getEXPMeter`
* `utils/gen_template_entry.py` calls `ds.makeResult`, `ds.makeScriptobsLine`
  (the `.sin` is the live version)

None of those names existed in `UCOScheduler.py` before this refactor, so these
scripts were already dead. (`gen_template_entry.py` also calls
`ds.find_Bstars`, which did exist until step 3 moved it to `Observability`; the
script was already broken by the other two calls.)

## Outstanding bugs

These are pre-existing. They are fixed one commit each, outside the move-only
steps, because each can change which target gets picked. Line numbers are
current as of `c308960`. Bugs already fixed are listed under "Steps as done".

1. **`star_table['too'] is False`** at `SunPos.py:80` (`sun_el_check`). `is` on
   a numpy array is always `False`, so `faint &= False` makes `faint`
   all-`False` and faint stars are never rejected when the sun is above `-18`.
   Intended: `star_table['too'] == False`, as already done for the same line in
   `Observability.time_check` (`594aa50`).

2. **Disabled obsblock path.** `make_obs_block` (`ScriptobsLine.py:153`) has no
   live callers; the only call is commented out in `make_result` at
   `UCOScheduler.py:184`, part of the commented-out `obsblock` lines at
   `:162-163` and `:182-184`. Proposal: keep `make_obs_block`, since
   re-enabling obsblocks is a plausible future want, and add a comment saying
   the calling path is disabled.

3. **`config_defaults` result is almost entirely unused**: `get_next` builds
   `config` at `UCOScheduler.py:288` and reads only `config['mode']`, which is
   `''`, at `:349`. `config_defaults` stays for
   `ParseUCOSched.parse_UCOSched`'s `config=` parameter, but the `get_next`
   call can drop it.

4. **`sim_nights.py` skips the first and last night of the range.**
   `gen_datelist` (`utils/sim_nights.py:40`) advances `cur` by a day *before*
   appending it, and loops `while cur < end`, so both the start and end dates
   are dropped. `sim_nights.py 2026/09/28 2026/10/01` simulates only 09/29 and
   09/30 (two "sun rose" lines, not four). Intended: include both endpoints.
   This changes the baseline's night count, so fix it before recording the
   `sim_nights.py` baseline that the other fixes are diffed against.

5. **`compute_datetime` silently uses the current time for an `int`.**
   `Observability.compute_datetime` (`Observability.py:22`) accepts a `float`,
   `datetime` or `ephem.Date`, and falls through to `utcnow()` for anything
   else, including an `int` Unix timestamp such as `calendar.timegm(...)`
   returns. No error or log line, just the wrong time. Production is not
   affected (`Observe` passes `time.time()`, a `float`), but it is an easy trap
   in the sim scripts and tests. Intended: accept any real number
   (`isinstance(ctime, (int, float))`, excluding `bool`).

**Order.** Bug 4 first, then record the `sim_nights.py` baseline. Then bug 1,
which changes target selection, diffed against that baseline. Bugs 2, 3 and 5
last; they should not change target selection.

## Steps as done

Each step left the tree runnable, so a bad step can be bisected.

| Step | What | Commit(s) |
| --- | --- | --- |
| (pre) | Deleted `compute_preferred_el`; fixed cal-star boosting by passing `observed=` to `compute_priorities` | `f4ae698` |
| 1 | Constants `ACQUIRE`/`BLANK`/`FIRST`/`LAST`/`BUFFERSEC`/`BUFFER` moved to `SchedulerConsts.py` | `d092459` |
| 2 | `ScriptobsLine.py`; updated `utils/make_scriptobsline.py` and `utils/gen_template_entry.sin` | `d092459` |
| 3 | `Observability.py`; pure code motion | `9fe7518` |
| 4 | `test_UCOScheduler.py`; tests moved out, `time_left` path made local | `ec5fd9a` |
| 5 | `UCOTargets` renamed to `UCOTargetTables` (`git mv`), then `update_hour_table`, `update_from_observed`, `gen_stars`, `check_files` added. All importers updated, including `Observe.py`, the sim scripts, `gen_template_entry.sin` and the tests, which the original plan did not list | `c1805f5`, `3a7cf23` |
| 6 | `UCOScheduler` class; `get_next` split into stages | `7d5af65` |
| 7 | `Observe.py`, `Main.sin`, `sim_night.py`, `sim_nights.py`, `test_UCOScheduler.py` moved to the object API; temporary shim deleted | `8f0f51f` |
| merge | Parallel line (see below) merged in | `8ac7dc7` |
| 8 | Fixes lost by the merge reapplied to the class code; `track_failures` option added | `c308960` |
| 8 | Template budget moved into the scheduler (`tot_temps`, `n_temps`) | uncommitted |
| 8 | Remaining bug fixes | outstanding |

**Step 6 shim.** While `Observe` and the sim scripts still called the module
functions, step 6 kept a module-level `get_next(ctime, seeing, slowdown,
targets, **kw)` and `zero_last_objs_attempted()`. The original plan had the
shim construct a new scheduler per call. That would have emptied
`last_objs_attempted` on every call in production, so the shim instead kept
one shared module-level scheduler, which behaved exactly like the old global.
Step 7 deleted it.

**Parallel line and merge.** A second line of work branched from step 3
(`9fe7518`). It did its own version of steps 4-5 (`75581d0`, with the plan
update `af2ae93`), merged `main` (`afd3357`), and fixed bugs against the old
module-level `get_next` and the old `update_from_observed` (`1507451`,
`994c692`, `ae57130`). The merge `8ac7dc7` kept the class code from the steps
4-7 line, which silently dropped `1507451` and `994c692`, and `main`'s
`3010802` (which had turned failed-object tracking off). `ae57130` survived
because it touched `Observability.py`, which both lines shared. `c308960`
reapplied the two lost fixes to the class code; tracking became the
`track_failures` option instead.

**Bugs fixed:**

| Bug | Fix | Commit(s) |
| --- | --- | --- |
| Cal-star boosting never ran: `get_next` did not pass `observed=` to `compute_priorities` | pass it | `f4ae698` |
| `star_table['too'] is False` in `time_check` made the faint-star exposure limit never apply | `== False` (`Observability.py:108`); the same bug in `SunPos.py` is still open (outstanding bug 1) | `594aa50` |
| `np.any(...) is False` made the "no B stars" and "not enough time" returns unreachable; the "no B stars" `apflog` call passed `label=`, which would have raised `TypeError` | `not np.any(...)`, `level='error'` | `1507451`, reapplied in `c308960` |
| `vstack(too_table, star_table)` raised whenever `too.dat` existed | `vstack([too_table, star_table])` | `bc1262d`, `3c2f7bd` |
| A missing `googledex.dat` replaced the in-memory star table with `None` and lost the night's observations for that call | rebuild the file and re-apply the observations; keep the previous table if the rebuild fails | `994c692`, reapplied in `c308960` |
| `dt.strftime('%s')` in `time_check` was non-portable and an hour off across a DST change | `calendar.timegm(dt.utctimetuple())` | `ae57130` |
| Failed-object tracking | made optional: `track_failures`, off in `Main.sin` | `3010802` (turned off), `c308960` (option) |

## How this was verified

There is no test suite, so verification is differential: for each step, the
previous commit and the new tree were run side by side, each on its own fresh
copy of the fixtures, and their printed output and written files compared.

**Fixtures.** `/Users/holden/src/test_main` holds `googledex.dat`,
`rank_table` and `time_left.csv`. It is outside the repo. Every run copies
these first, because the scheduler rewrites `googledex.dat` and creates
`observed_targets` and `hour_table` in the current directory.

**What was run:**

* **Steps 4-7:** `test_UCOScheduler.py`, old and new started at the same moment
  (the tests read the clock). Printed output and the written `googledex.dat`,
  `observed_targets`, `hour_table` and `rank_table` were identical at every
  step.
* **Step 6, a wider sweep:** 176 `get_next` calls, hourly through two nights
  (2026-10-10, near new moon, and 2026-09-27), with every combination of
  `bstar`, `do_templates` and `do_too`. It faked a failed observation every few
  calls and called `zero_last_objs_attempted` partway through. Output
  (4,855 lines), written files and the final `last_objs_attempted` list were
  identical. A second pass with every target set to `Template=N` ran the
  template branch 8 times; also identical.
* **Step 7, the sim scripts:** run through a small wrapper that seeds numpy's
  random generator identically, since `NightSim` draws seeing, clouds and
  count rates from `np.random`. `sim_night.py -d 2026-09-28`: `.simout`
  (57 observations), `observed_targets`, `googledex.dat`, `hour_table` and
  735 lines of output identical. `sim_nights.py 2026/09/28 2026/10/01`:
  `.simout` (159 lines), every written file and 2,554 lines of output
  identical.
* **`c308960`:** `test_UCOScheduler.py` output identical to the merge. Each
  reapplied path was run on its own: no B stars with `bstar=True` returns
  `None` without a `TypeError`; `time_check` rejecting everything returns
  `None` at the "not enough time" check; a missing `googledex.dat` is rebuilt
  and the night's observation applied (`nobs` 25 to 26), and a failed rebuild
  keeps the previous table. With `last_attempted()` faked to report the first
  pick as failed, `track_failures=True` skips it on the next call and
  `track_failures=False` never calls `last_attempted()`.

* **Template budget:** seeded `sim_night.py`, `sim_nights.py` and
  `test_UCOScheduler.py` identical to `c308960`. With every target set to
  `Template=N`, seeded `sim_night.py` on 2026-10-10 and 2026-10-12 returned
  exactly 2 templates each, old and new, with identical output and `.simout`.
  Over 30 calls on one night, `tot_temps=4` returned 4 templates and
  `tot_temps=None` returned 18.

**Not executed:** `Observe.py`, `Main.sin` and `getUCOTargets.py` need `ktl`,
so they were compiled or parsed, and grepped for leftover references, but never
run. The first night on the telescope is their real test.

**Not in the repo:** the sweep script and the numpy-seeding wrapper were ad hoc
and were not committed. No saved baseline exists yet. Because the tests read
the clock, a repeatable baseline needs either a fixed start time in
`test_main` or the seeded `sim_nights.py` run, recorded after outstanding
bug 4 is fixed.

## Open questions

1. Should `get_next` keep returning a plain `dict`, or become a small
   `Target` result class? The dict is consumed in `Observe.py` by string key
   (`self.target['NAME']`, `self.target["SCRIPTOBS"]`) and in the sim scripts.
   A result class is nicer but widens the diff. **As built:** still a dict.
   The plan for replacing it is in `docs/scheduler_target_class.md`.

Settled: the tables class name (`UCOTargetTables`), failed-object tracking
(the `track_failures` option) and the template budget (on the scheduler).
