# Scheduler refactor: turning UCOScheduler into a class

Status: steps 1-5 and the bug fixes are done; steps 6-7 (the class itself and
porting the call sites) are next.
Branch: `scheduler_object`, with `main` merged in at `afd3357`.

## Progress

| Step | State | Commit |
| --- | --- | --- |
| 1. Constants to `SchedulerConsts.py` | done | `d092459` |
| 2. `ScriptobsLine.py` | done | `d092459` |
| 3. `Observability.py` | done | `9fe7518` |
| 4. `test_UCOScheduler.py` | done | `75581d0` |
| 5. `UCOTargetTables.py` | done | `75581d0` |
| 6. `UCOScheduler` class | not started | |
| 7. Port the call sites | not started (sim scripts renamed only) | |
| 8. Bug fixes | done except bug 8, which is optional | see "Bugs found while reading" |

## Goal

`Main/UCOScheduler.py` was a 1147-line module of free functions with a
module-level global. Steps 1-5 took it down to 478 lines, but it is still free
functions. Every call to `get_next()` threads `ucotargets`, `focval`, `owner`,
`outdir` and `start_time` through by hand.

After this refactor:

* `UCOScheduler` is a class. `Observe` holds an instance and calls
  `self.scheduler.get_next(...)`.
* Target-table state and the files those tables live in stay in one object,
  separate from the scheduler's own run-state.
* The pure helpers (scriptobs line formatting, visibility and condition cuts)
  move to their own modules so they can be tested without a scheduler. (Done.)

## Decisions taken

| Question | Decision |
| --- | --- |
| Class relationship | Composition: `UCOScheduler` holds a `UCOTargetTables` (the old `UCOTargets`, renamed) |
| Entry point name | `get_next()` (snake_case, matches the repo) |
| Tables class name | `UCOTargetTables` (keeps the package's `UCO` prefix) |
| File layout | Split into four modules plus a test module |
| `utils/` | Port the live scripts; leave the already-dead ones alone |

## Why composition and not one merged class

The tables object is not owned by the scheduler. `Main.sin` builds it and
hands the *same* object to two places:

```
Main.sin:488   uco_targets = UCOTargetTables.UCOTargetTables(opt)
Main.sin:489   getUCOTargets.getUCOTargets(uco_targets, ...)   # background thread
Main.sin:490   observe = Observe.Observe(apf, tel, opt, uco_targets, ...)
```

`getUCOTargets` is a `threading.Thread` that waits on sun elevation and then calls
`make_hour_table()` / `make_star_table()`, i.e. it downloads the Google sheets and
writes them to disk on its own schedule. Merging the tables into `UCOScheduler`
would mean that download thread mutates the scheduler object directly, and it
would drag the selection algorithm into `getUCOTargets`'s import graph. Keeping
them as two objects preserves the current separation: one object owns the
tables and the files, the other owns target selection.

## Lifetime: build the scheduler in Main.sin, not in Observe

This matters and is easy to get wrong. `Main.sin:515-522` **recreates**
`Observe` when its thread dies:

```
if observe.is_alive() is False:
    observe = Observe.Observe(apf, tel, opt, uco_targets, task=parent)
    observe.start()
```

`uco_targets` survives that restart, and so does `last_objs_attempted`, because
it is a module global. If the scheduler were constructed inside
`Observe.__init__`, the list of objects that just failed to observe would be
wiped every time the thread restarted, and the scheduler would immediately
re-select a target it had already failed on.

`main` has currently turned the failure tracking off: the `last_attempted()`
call in `get_next` (`UCOScheduler.py:297-299`) is commented out, so nothing is
appended to `last_objs_attempted` at the moment. The lifetime argument still
holds for when it is turned back on.

So: `Main.sin` constructs the scheduler alongside the tables and passes it to
`Observe`, exactly as it does with `uco_targets` now.

```python
# Main.sin
targets   = UCOTargetTables.UCOTargetTables(opt)
_         = getUCOTargets.getUCOTargets(targets, task=parent, wait_time=target_time)
scheduler = UCOScheduler.UCOScheduler(targets, opt)
observe   = Observe.Observe(apf, tel, opt, scheduler, task=parent)
```

`Observe` reaches the tables through `self.scheduler.targets` on the rare
occasions it needs them, so it takes one argument instead of two.

## File layout

| File | Contents | State |
| --- | --- | --- |
| `Main/SchedulerConsts.py` | gained `ACQUIRE`, `BLANK`, `FIRST`, `LAST`, `BUFFERSEC`, `BUFFER` | done |
| `Main/UCOTargetTables.py` (from `UCOTargets.py`) | `class UCOTargetTables`: rank/hour/star tables and all disk bookkeeping | done |
| `Main/ScriptobsLine.py` | pure string generation for scriptobs lines | done |
| `Main/Observability.py` | pure array/astro filters, no scheduler state | done |
| `Main/UCOScheduler.py` | `class UCOScheduler` only | still free functions |
| `Main/test_UCOScheduler.py` | the `test_*` functions from the bottom of `UCOScheduler.py` | done |

`UCOTargets.py` is gone (`git mv` to `UCOTargetTables.py`, so history follows).

## Where every symbol went

### `Main/ScriptobsLine.py` (done)

| From `UCOScheduler.py` | Note |
| --- | --- |
| `make_scriptobs_line` | unchanged; `utils/make_scriptobsline.py` and `gen_template_entry.sin` call it |
| `num_template_exp` | unchanged |
| `config_defaults` | unchanged |
| `make_obs_block` | **dead**, kept with a comment saying so. Its only caller is the commented-out obsblock path in `make_result` (`UCOScheduler.py:197`) |

### `Main/Observability.py` (done)

All pure functions of `(star_table, moon, apf_obs, dt, ...)`, no instance state:

`compute_datetime`, `tot_exp_times`, `time_check`, `condition_cuts`,
`behind_moon`, `template_conditions`, `find_closest`, `find_Bstars`,
`enough_time_templates`.

(`compute_preferred_el` was on this list; it had no callers and was deleted in
`f4ae698`.)

These are the functions worth having unit tests for, which is the main reason to
pull them out.

### `Main/UCOTargetTables.py` (done)

Everything that reads or writes the bookkeeping files:

| Method | From |
| --- | --- |
| `make_rank_table` | `UCOTargets` (unchanged) |
| `make_hour_constraints` | `UCOTargets` (unchanged) |
| `make_hour_table` | `UCOTargets` (unchanged) |
| `make_star_table` | `UCOTargets` (unchanged) |
| `append_too_column` | `UCOTargets` (unchanged) |
| `copy_backup` | `UCOTargets` (unchanged) |
| `update_hour_table(observed, dt)` | **moved** from `UCOScheduler.update_hour_table`. It writes `hour_table` to disk, so it is bookkeeping. Does nothing if there is no hour table |
| `update_from_observed(ptime)` | **new**: wraps `ParseUCOSched.update_local_starlist`, sets `self.star_table`, returns the `ObservedLog`. If `googledex.dat` is missing it rebuilds it and then applies the observations (bug 5) |
| `gen_stars()` | **new**: wraps `ParseUCOSched.gen_stars` and keeps the result as `self.stars` |
| `check_files()` | **moved** from `Observe.check_files`. Restores `googledex.dat` from its `.1` backup if it is missing |

`get_next` already calls `update_from_observed`, `update_hour_table` and
`gen_stars` on the tables object.

Files this object owns, unchanged in name and format:

* `rank_table` (ascii)
* `hour_table` (ascii)
* `googledex.dat` + `googledex.dat.1` backup (ecsv)
* `too.dat` (ecsv)
* `observed_targets` (read via `ObservedLog`; written by `Observe`)
* the time-left CSV named by `opt.time_left` (read only)

Nothing about these files changes. That is the constraint that makes this
refactor safe to deploy: the on-disk contract is identical, so a night can be
resumed from files written by the old code.

### `Main/UCOScheduler.py` (step 6, not started)

What is left in the module today: the `last_objs_attempted` global,
`zero_last_objs_attempted`, `need_cal_star`, `compute_priorities`,
`make_result`, `last_attempted` and `get_next` (`UCOScheduler.py:242`, about
235 lines). The target shape:

```python
class UCOScheduler(object):

    def __init__(self, targets, opt=None, owner='public', outdir=None,
                 do_templates=True, do_too=True, start_time=None,
                 outfn='googledex.dat', toofn='too.dat'):
        self.targets   = targets          # UCOTargetTables
        self.owner     = owner
        self.outdir    = outdir or os.getcwd()
        self.outfn     = outfn
        self.toofn     = toofn
        self.do_templates = do_templates
        self.do_too       = do_too
        self.start_time   = start_time

        # run-state, was a module global
        self.last_objs_attempted = []

        # per-call scratch, kept for logging and for Observe to inspect
        self.observed  = None             # ObservedLog from the last refresh
        self.apf_obs   = None
        self.moon      = None
        self.result    = None
        self.template_conditions_met = False

    # --- public -------------------------------------------------------
    def get_next(self, ctime, seeing, slowdown, bstar=False, focval=0,
                 do_templates=None, do_too=None):
    def zero_last_objs_attempted(self):
    def record_last_attempt(self):        # was module-level last_attempted(); call is disabled on main

    # --- pipeline stages, one per current block of get_next ------------
    def _previous_obs_time(self, dt):     # apfguide midptfin, falls back to dt
    def _refresh_tables(self, ptime):     # update_from_observed + hour table
    def _sky_state(self, dt):             # apf_obs + moon
    def _available(self, dt, seeing, slowdown, bstar)   -> bool array
    def compute_priorities(self, dt)      # uses self.observed for need_cal_star
    def _need_cal_star(self, priorities)
    def _select(self, available, priorities, bstar)     -> idx
    def _make_result(self, idx, totexptimes, priorities, dt, focval, bstar)
    def _add_template(self, res, idx, dt, bstars)
```

The ephem star list lives on `self.targets.stars` (set by `gen_stars()`), so the
scheduler does not keep its own copy.

`get_next` becomes roughly forty lines calling the stages in order. Each stage
keeps its existing `apflog` calls so the night log reads the same.

Constructor arguments that are currently per-call `get_next` keyword arguments
(`owner`, `outdir`, `outfn`, `toofn`, `start_time`, `do_templates`, `do_too`)
move to `__init__`, since `Observe` passes the same value on every call.
`do_templates` and `do_too` are also accepted per call as overrides, because
`Observe` mutates `self.do_temp` and `self.do_too` during a night.

### `Main/test_UCOScheduler.py` (done)

`test_basic_ops`, `test_failure`, `test_templates` and `test_main` moved here
verbatim, except that they now call `UCOScheduler.get_next(...)`. They stay
runnable as a script (`python test_UCOScheduler.py`); converting them to pytest
is a separate job. Step 7 changes them to the object API.

## Observe.py changes

Done so far: the `UCOTargetTables` rename, and `check_files()` moved to
`UCOTargetTables.check_files()`. `Observe` now calls
`self.uco_targets.check_files()` (`Observe.py:511`). That becomes
`self.scheduler.targets.check_files()` in step 7.

Remaining for step 7, all mechanical:

| Line | Now | After |
| --- | --- | --- |
| `Observe.py:21-22` | `import UCOScheduler as ds` / `import UCOTargetTables` | `import UCOScheduler` |
| `Observe.py:33` | `def __init__(self, apf, tel, opt, uco_targets, ...)` | `def __init__(self, apf, tel, opt, scheduler, ...)` |
| `Observe.py:511` | `self.uco_targets.check_files()` | `self.scheduler.targets.check_files()` |
| `Observe.py:513` | `ds.get_next(time.time(), seeing, slowdown, self.uco_targets, bstar=..., do_too=..., owner=..., do_templates=..., focval=..., start_time=...)` | `self.scheduler.get_next(time.time(), seeing, slowdown, bstar=self.obs_B_star, focval=self.focval, do_templates=self.do_temp, do_too=self.do_too)` |
| `Observe.py:646` | `ds.zero_last_objs_attempted()` | `self.scheduler.zero_last_objs_attempted()` |
| `Observe.py:1132-1134` (test `main`) | builds `UCOTargetTables` and passes it | builds the tables and a scheduler, passes the scheduler |

`Observe.start_time` is read by `get_next` today via the `start_time=` keyword.
It moves to the scheduler constructor. `Observe` clears `self.start_time` in
several places: `__init__` (`:89`, `:91`), `should_start_list()` (`:299`,
`:311`, after an hour), and `:716`, `:755`, `:767`, `:877`, `:883`. Every one
of those must clear `self.scheduler.start_time` instead, or the scheduler keeps
limiting exposure times to a start window that has passed. Easiest is to make
`start_time` the scheduler's attribute only, and have `Observe` read and write
`self.scheduler.start_time`.

## getUCOTargets.py changes (done)

Import and type name only: `UCOTargets.UCOTargets` became
`UCOTargetTables.UCOTargetTables`. The thread does not touch the scheduler.

## utils/ changes

Done:

* `utils/make_scriptobsline.py` calls `ScriptobsLine.make_scriptobs_line` and
  no longer imports `UCOScheduler`.
* `utils/gen_template_entry.sin` calls `Observability.find_Bstars` and
  `ScriptobsLine.make_scriptobs_line`, and no longer imports `UCOScheduler`.
* `utils/download_googledex.py` uses `UCOTargetTables`.
* `utils/sim_night.py` and `utils/sim_nights.py` use `UCOTargetTables` (rename
  only).

Remaining for step 7:

* `utils/sim_night.py:144`, `utils/sim_nights.py:230`: build a `UCOScheduler`
  once before the loop and call `scheduler.get_next(...)` inside it.

Left alone, already broken against the deployed code:

* `utils/calc_cadence.py` calls `ds.parseGoogledex` and `ds.DS_APFPRI`
* `utils/calc_etime_precision.py` calls `ds.get_speadsheet`, `ds.getI`, `ds.DS_BV`
* `utils/calc_precision_mag.py` calls `ds.getI`, `ds.getEXPMeter`
* `utils/gen_template_entry.py` calls `ds.makeResult`, `ds.makeScriptobsLine`
  and `ds.find_Bstars` (the `.sin` is the live version)

None of those names exist in `UCOScheduler.py`. They are dead, and this
refactor does not make them any deader.

## Bugs found while reading

These were pre-existing. Each fix went in as its own commit, separate from the
code moves.

1. **Cal-star boosting never ran.** *Fixed in `f4ae698`.* `get_next` called
   `compute_priorities(...)` without `observed=`, so `need_cal_star` returned
   the priorities unchanged. `get_next` now passes `observed=observed`
   (`UCOScheduler.py:406-410`). Checked by hand: with an observed program
   marked `need_cal`, its two cal stars are raised to max priority + 1, and are
   no longer raised once one of them has been observed. The 2025A fixtures have
   no `need_cal` or `cal_star` rows, so `test_main` does not exercise this.

2. **`star_table['too'] is False`** in `time_check`. *Fixed on `main` in
   `594aa50`, applied to `Observability.py:108` in the merge `afd3357`.* Now
   `star_table['too'] == False`. Before, the faint-star limit (`maxfaintexptime`,
   the -18° sunrise) never applied. On the fixtures it now covers 519 faint
   targets; at 13:30 UT none of them pass while 290 bright ones still do.

3. **`np.any(...) is False`**. *Fixed in `1507451`.* Both early returns in
   `get_next` (`UCOScheduler.py:321`, `:372`) now use `not np.any(...)`. The
   no-B-stars branch also called `apflog(..., label='Error')`, which neither
   `apflog` nor `fake_apflog` accepts, so the first time the branch was reached
   it would have raised `TypeError`. That is now `level='error'`.

4. **`vstack` called with two positional arguments** in
   `ParseUCOSched.update_local_starlist`. *Fixed on `main` in `594aa50`,
   merged in `afd3357`.* Now `vstack([too_table, star_table])`.

5. **Missing `googledex.dat` lost the night's observations.** *Fixed in
   `994c692`.* `update_local_starlist` returns `None` for the star table when
   the file is missing. `get_next` then rebuilt the table, but observations are
   only applied through the file, so the rebuilt table had none of them. If the
   rebuild failed too, `get_next` raised `AttributeError`.
   `UCOTargetTables.update_from_observed` now rebuilds the file and then
   applies the observations. If the rebuild fails, it keeps the previous table.

   The fix first proposed here, keeping the in-memory table whenever the file is
   missing, was not used. Nothing would recreate the file, so the table would
   stop picking up observations for the rest of the night.

6. **`dt.strftime('%s')`** in `time_check`. *Fixed in `ae57130`.* Now
   `calendar.timegm(dt.utctimetuple())`. The old arithmetic was an hour off
   whenever `dt` and the current time were on opposite sides of a DST change.
   That happens on the night the clocks change, and whenever a sim replays a
   night from the other half of the year. For a January night run from
   California with a start list 30 minutes away, the old code ignored the
   window completely.

7. **Dead code.** *Done.* `compute_preferred_el` was deleted in `f4ae698`.
   `make_obs_block` is kept in `ScriptobsLine.py` with a comment saying its
   calling path is disabled, since re-enabling obsblocks is a plausible future
   want.

8. **`config_defaults` result is almost entirely unused.** *Not changed.*
   `get_next` builds `config` (`UCOScheduler.py:260`) and reads only
   `config['mode']` (`:455`), which is `''`. `config_defaults` stays for
   `ParseUCOSched.parse_UCOSched`'s `config=` parameter. The `get_next` call can
   go when step 6 rewrites that code.

Found during the merge, not a bug fix: `main` commented out the
`last_attempted()` call (`3010802`, "turn off last_obj_attempted"), so the
scheduler no longer avoids a target that just failed. The merge keeps that.

## Implementation order

Each step leaves the tree runnable, so a bad step can be bisected.

1. ~~**Constants.**~~ Done.
2. ~~**`ScriptobsLine.py`.**~~ Done.
3. ~~**`Observability.py`.**~~ Done.
4. ~~**`test_UCOScheduler.py`.**~~ Done.
5. ~~**`UCOTargetTables.py`.**~~ Done.
6. **`UCOScheduler` class.** Write the class and split `get_next` into the
   stages above. Keep a module-level `get_next(ctime, seeing, slowdown,
   targets, **kw)` that builds a scheduler and delegates to it, *temporarily*,
   so `sim_night.py` still runs unchanged and can be diffed against the
   baseline. Drop the unused `config` from `get_next` (bug 8) while there.
7. **Port the call sites.** `Observe.py`, `Main.sin`, `sim_night.py`,
   `sim_nights.py` and `test_UCOScheduler.py` to the object API, including the
   `start_time` handling above. Delete the temporary module-level shim.

## How this gets verified

There is no test suite, so verification is differential: the old and new code
run on the same inputs, and the outputs are compared.

What was done for steps 1-5 and the bug fixes:

* **Inputs**: the 2025A tables in `~/Dropbox/src/test_UCOmain` (`rank_table`,
  `googledex.dat` and their `.1` backups, `hour_table`, `time_left.csv`).
  Each run copies them into a fresh scratch directory.
* **Run**: `test_UCOScheduler.test_main()` with `time.time()` pinned to
  2025-03-11 06:00 UT. This uses the backup-file path; there is no Google
  certificate in this checkout.
* **Refactor steps**: `observed_targets`, every table file and the log were
  byte-identical to the previous code. The one exception is the uth/utm on the
  template lines, which `test_templates` takes from the wall clock.
* **Bug fixes**: `test_main` was unchanged by every fix, because none of them
  change the fixtures' normal path. Each fix was checked separately on the
  case it affects, as described under each bug above.

Still to do:

* **`sim_nights.py` baseline**: run it over a fixed set of nights on the
  current branch head and save the full sequence of selected scriptobs lines.
  Steps 6-7 must reproduce it exactly. It has not been run yet.
* **Fixtures in the repo**: copy the pinned inputs and the time-pinning runner
  into a fixtures directory, so the comparison doesn't depend on a Dropbox path.
* **Google-sheet path**: `test_UCOScheduler.py` has only run on the backup path.
  Running it where the certificate is available (e.g. the deployed machine)
  would cover the download path.
* **`Observe` and `Main.sin`** need `ktl`, so nothing above exercises them.
  Step 7's changes there have to be checked on the telescope machine or in
  `--test` mode.

## Open questions

1. Should `get_next` keep returning a plain `dict`, or become a small
   `Target` result class? The dict is consumed in `Observe.py` by string key
   (`self.target['NAME']`, `self.target["SCRIPTOBS"]`) and in the sim scripts.
   A result class is nicer but widens the diff; the plan above keeps the dict.
2. `Observe` mutates `self.do_temp` and `self.n_temps` during the night to cap
   template observations at `tot_temps`. That template budget arguably belongs
   to the scheduler. The plan leaves it in `Observe` for now and passes
   `do_templates` per call.
3. Should failure tracking (`last_attempted()`) come back on? If so, the
   class's `record_last_attempt` should be called from `get_next` again.
