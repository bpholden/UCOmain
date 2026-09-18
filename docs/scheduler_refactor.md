# Scheduler refactor: turning UCOScheduler into a class

Status: draft plan, not yet implemented.
Branch: clean branch off the deployed version.

## Goal

`Main/UCOScheduler.py` is a 1147-line module of free functions with a module-level
global. Every call to `get_next()` threads `star_table`, `rank_table`, `hour_table`,
`observed`, `focval`, `owner`, `outdir` and `start_time` through by hand.

After this refactor:

* `UCOScheduler` is a class. `Observe` holds an instance and calls
  `self.scheduler.get_next(...)`.
* Target-table state and the files those tables live in stay in one object,
  separate from the scheduler's own run-state.
* The pure helpers (scriptobs line formatting, visibility and condition cuts)
  move to their own modules so they can be tested without a scheduler.

## Decisions taken

| Question | Decision |
| --- | --- |
| Class relationship | Composition: `UCOScheduler` holds a `TargetTables` (today's `UCOTargets`, renamed) |
| Entry point name | `get_next()` (snake_case, matches the repo) |
| File layout | Split into four modules plus a test module |
| `utils/` | Port the live scripts; leave the already-dead ones alone |

## Why composition and not one merged class

`UCOTargets` is not owned by the scheduler today. `Main.sin` builds it and hands
the *same* object to two places:

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
them as two objects preserves the current separation: one object owns the
tables and the files, the other owns target selection.

## Lifetime: build the scheduler in Main.sin, not in Observe

This matters and is easy to get wrong. `Main.sin:522` **recreates** `Observe` when
its thread dies:

```
if observe.is_alive() is False:
    observe = Observe.Observe(apf, tel, opt, uco_targets, task=parent)
    observe.start()
```

`uco_targets` survives that restart; `last_objs_attempted` also survives today
because it is a module global. If the scheduler were constructed inside
`Observe.__init__`, the list of objects that just failed to observe would be
silently wiped every time the thread restarted, and the scheduler would
immediately re-select a target it had already failed on.

So: `Main.sin` constructs the scheduler alongside the tables and passes it to
`Observe`, exactly as it does with `uco_targets` now.

```python
# Main.sin
targets   = TargetTables.TargetTables(opt)
_         = getUCOTargets.getUCOTargets(targets, task=parent, wait_time=target_time)
scheduler = UCOScheduler.UCOScheduler(targets, opt)
observe   = Observe.Observe(apf, tel, opt, scheduler, task=parent)
```

`Observe` reaches the tables through `self.scheduler.targets` on the rare
occasions it needs them, so it takes one argument instead of two.

## New file layout

| File | Contents |
| --- | --- |
| `Main/SchedulerConsts.py` (existing) | gains `ACQUIRE`, `BLANK`, `FIRST`, `LAST`, `BUFFERSEC`, `BUFFER` |
| `Main/TargetTables.py` (new, from `UCOTargets.py`) | `class TargetTables` — rank/hour/star tables and all disk bookkeeping |
| `Main/ScriptobsLine.py` (new) | pure string generation for scriptobs lines |
| `Main/Observability.py` (new) | pure array/astro filters — no scheduler state |
| `Main/UCOScheduler.py` (rewritten) | `class UCOScheduler` only |
| `Main/test_UCOScheduler.py` (new) | the `test_*` functions currently at the bottom of `UCOScheduler.py` |

`UCOTargets.py` is deleted; `TargetTables.py` replaces it. (Naming is the one
detail still open — `TargetTables` reads well but drops the `UCO` prefix the rest
of the package uses. `UCOTargetTables` is the alternative.)

## Where every current symbol goes

### `Main/ScriptobsLine.py`

| From `UCOScheduler.py` | Note |
| --- | --- |
| `make_scriptobs_line` | unchanged; `utils/` and `gen_template_entry.sin` call it |
| `num_template_exp` | unchanged |
| `config_defaults` | unchanged |
| `make_obs_block` | **dead** — its only caller is a commented-out block at `UCOScheduler.py:696`. Moved but marked, or dropped; see "Dead code" below |

### `Main/Observability.py`

All pure functions of `(star_table, moon, apf_obs, dt, ...)`, no instance state:

`compute_datetime`, `tot_exp_times`, `time_check`, `condition_cuts`,
`behind_moon`, `template_conditions`, `find_closest`, `find_Bstars`,
`enough_time_templates`, `compute_preferred_el` (**dead**, see below).

These are the functions worth having unit tests for, which is the main reason to
pull them out.

### `Main/TargetTables.py`

Everything that reads or writes the bookkeeping files:

| Method | From |
| --- | --- |
| `make_rank_table` | `UCOTargets` (unchanged) |
| `make_hour_constraints` | `UCOTargets` (unchanged) |
| `make_hour_table` | `UCOTargets` (unchanged) |
| `make_star_table` | `UCOTargets` (unchanged) |
| `append_too_column` | `UCOTargets` (unchanged) |
| `copy_backup` | `UCOTargets` (unchanged) |
| `update_hour_table(observed, dt)` | **moved** from `UCOScheduler.update_hour_table` — it writes `hour_table` to disk, so it is bookkeeping |
| `update_from_observed(ptime)` | **new** — wraps `ParseUCOSched.update_local_starlist`, sets `self.star_table`, returns the `ObservedLog` |
| `gen_stars()` | **new** — thin wrapper on `ParseUCOSched.gen_stars`, caches the ephem objects alongside the table that produced them |

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

### `Main/UCOScheduler.py`

```python
class UCOScheduler(object):

    def __init__(self, targets, opt=None, owner='public', outdir=None,
                 do_templates=True, do_too=True, start_time=None,
                 outfn='googledex.dat', toofn='too.dat'):
        self.targets   = targets          # TargetTables
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
        self.stars     = None
        self.result    = None
        self.template_conditions_met = False

    # --- public -------------------------------------------------------
    def get_next(self, ctime, seeing, slowdown, bstar=False, focval=0,
                 do_templates=None, do_too=None):
    def zero_last_objs_attempted(self):
    def record_last_attempt(self):        # was module-level last_attempted()

    # --- pipeline stages, one per current block of get_next ------------
    def _previous_obs_time(self, dt):     # apfguide midptfin, falls back to dt
    def _refresh_tables(self, ptime):     # update_from_observed + hour table
    def _sky_state(self, dt):             # apf_obs + moon
    def _available(self, dt, seeing, slowdown, bstar)   -> bool array
    def compute_priorities(self, dt)      # uses self.observed, so need_cal_star works
    def _need_cal_star(self, priorities)
    def _select(self, available, priorities, bstar)     -> idx
    def _make_result(self, idx, totexptimes, priorities, dt, focval, bstar)
    def _add_template(self, res, idx, dt, bstars)
```

`get_next` becomes roughly forty lines calling the stages in order, instead of
the current 240-line function. Each stage keeps its existing `apflog` calls so
the night log reads the same.

Constructor arguments that are currently per-call `get_next` keyword arguments
(`owner`, `outdir`, `outfn`, `toofn`, `start_time`, `do_templates`, `do_too`)
move to `__init__`, since `Observe` passes the same value on every call.
`do_templates` and `do_too` are also accepted per call as overrides, because
`Observe` mutates `self.do_temp` and `self.do_too` during a night.

### `Main/test_UCOScheduler.py`

`test_basic_ops`, `test_failure`, `test_templates`, `test_main` move here
verbatim apart from the call-site changes. They stay runnable as a script
(`python test_UCOScheduler.py`); converting them to pytest is a separate job.

## Observe.py changes

Three call sites, all mechanical:

| Line | Now | After |
| --- | --- | --- |
| `Observe.py:22-23` | `import UCOScheduler as ds` / `import UCOTargets` | `import UCOScheduler` |
| `Observe.py:33` | `def __init__(self, apf, tel, opt, uco_targets, ...)` | `def __init__(self, apf, tel, opt, scheduler, ...)` |
| `Observe.py:526` | `ds.get_next(time.time(), seeing, slowdown, self.uco_targets, bstar=..., do_too=..., owner=..., do_templates=..., focval=..., start_time=...)` | `self.scheduler.get_next(time.time(), seeing, slowdown, bstar=self.obs_B_star, focval=self.focval, do_templates=self.do_temp, do_too=self.do_too)` |
| `Observe.py:659` | `ds.zero_last_objs_attempted()` | `self.scheduler.zero_last_objs_attempted()` |

`Observe.check_files()` (`Observe.py:295`) restores `googledex.dat` from its
`.1` backup and duplicates `TargetTables.copy_backup`. It moves to
`TargetTables.check_files()` and `Observe` calls
`self.scheduler.targets.check_files()`.

`Observe.start_time` is read by `get_next` today via the `start_time=` keyword.
It moves to the scheduler constructor, and `Observe.should_start_list()`
(`Observe.py:312`) — which clears `self.start_time` after an hour — must clear
`self.scheduler.start_time` instead, or the scheduler keeps constraining
exposure times to a start window that has passed. Easiest is to leave
`start_time` as the scheduler's attribute and have `should_start_list` read and
write `self.scheduler.start_time`.

## getUCOTargets.py changes

Import and type name only: `UCOTargets.UCOTargets` becomes
`TargetTables.TargetTables`. The thread does not touch the scheduler.

## utils/ changes

Live scripts, ported:

* `utils/sim_night.py:144`, `utils/sim_nights.py:230` — build a `UCOScheduler`
  once before the loop, call `scheduler.get_next(...)` inside it. These are the
  main regression harness for this refactor, so they get ported first.
* `utils/make_scriptobsline.py:44` — `ds.make_scriptobs_line` becomes
  `ScriptobsLine.make_scriptobs_line`.
* `utils/gen_template_entry.sin:77-83` — `ds.find_Bstars` becomes
  `Observability.find_Bstars`, `ds.make_scriptobs_line` becomes
  `ScriptobsLine.make_scriptobs_line`.
* `utils/download_googledex.py` — `UCOTargets` to `TargetTables`.

Left alone, already broken against the deployed code:

* `utils/calc_cadence.py` calls `ds.parseGoogledex` and `ds.DS_APFPRI`
* `utils/calc_etime_precision.py` calls `ds.get_speadsheet`, `ds.getI`, `ds.DS_BV`
* `utils/calc_precision_mag.py` calls `ds.getI`, `ds.getEXPMeter`
* `utils/gen_template_entry.py` calls `ds.makeResult`, `ds.makeScriptobsLine`
  (the `.sin` is the live version)

None of those names exist in `UCOScheduler.py` today. They are dead and this
refactor does not make them any deader.

## Bugs found while reading

These are pre-existing. I would **not** fix them silently inside a move-only
refactor — each one changes which target gets picked. Listing them so you can
decide which to take, and in which commit.

1. **Cal-star boosting never runs.** `get_next` calls `compute_priorities(...)`
   without `observed=`, so `need_cal_star` (`UCOScheduler.py:45`) hits its
   `if observed is None: return priorities` guard and returns unchanged
   priorities. The whole `need_cal` / `cal_star` mechanism is dead in
   production. On the class this is fixed by construction — `compute_priorities`
   reads `self.observed`, which `_refresh_tables` has already set. **This is a
   real behavior change** and should be its own commit, verified with
   `sim_nights.py`.

2. **`star_table['too'] is False`** at `UCOScheduler.py:272`. `is` on a numpy
   array is always `False`, so `faint &= False` makes `faint` all-`False` and
   the faint-star exposure-time branch (`maxfaintexptime`, the `-18` horizon)
   never applies. Intended: `~star_table['too']`.

3. **`np.any(...) is False`** at `UCOScheduler.py:842` and `:901`. `np.any`
   returns `np.bool_`, never the `False` singleton, so both early returns are
   unreachable — "no B stars listed" and "not enough time left to observe any
   targets" never fire. Intended: `not np.any(...)`.

4. **`vstack` called with two positional arguments** at
   `ParseUCOSched.py:944`: `astropy.table.vstack(too_table, star_table)`.
   `vstack`'s second positional is `join_type`, so this raises whenever a
   `too.dat` exists. Intended: `vstack([too_table, star_table])`.

5. **`update_local_starlist` returns `None` for the star table** when
   `googledex.dat` is missing, and `get_next` assigns that straight onto
   `ucotargets.star_table` (`UCOScheduler.py:827`), discarding a perfectly good
   in-memory table before rebuilding it. In `TargetTables.update_from_observed`
   the assignment becomes conditional.

6. **`dt.strftime('%s')`** at `UCOScheduler.py:238` is a glibc extension, not
   portable, and the surrounding UTC-offset arithmetic is fragile (its own
   comment says so). `calendar.timegm(dt.utctimetuple())` is the correct
   spelling and needs no offset correction.

7. **Dead code**: `compute_preferred_el` (`:383`) and `make_obs_block` (`:569`)
   have no live callers. The obsblock path in `make_result` is commented out at
   `:674-696`. Proposal: delete `compute_preferred_el`, and keep
   `make_obs_block` in `ScriptobsLine.py` with a comment saying the calling path
   is disabled, since re-enabling obsblocks is a plausible future want.

8. **`config_defaults` result is almost entirely unused** — `get_next` builds
   `config` and reads only `config['mode']`, which is `''`. It stays for
   `ParseUCOSched.parse_UCOSched`'s `config=` parameter, but the `get_next` call
   can drop it.

## Implementation order

Each step leaves the tree runnable, so a bad step can be bisected.

1. **Constants.** Move `ACQUIRE`/`BLANK`/`FIRST`/`LAST`/`BUFFERSEC`/`BUFFER` to
   `SchedulerConsts.py`. No behavior change.
2. **`ScriptobsLine.py`.** Move the four formatting functions out. Update
   `utils/make_scriptobsline.py` and `utils/gen_template_entry.sin`.
3. **`Observability.py`.** Move the pure filters out. Pure code motion.
4. **`test_UCOScheduler.py`.** Move the `test_*` functions out, still passing
   the module-level `get_next`. Run them — this is the baseline.
5. **`TargetTables.py`.** Rename `UCOTargets`, absorb `update_hour_table`, add
   `update_from_observed`, `gen_stars`, `check_files`. Update
   `getUCOTargets.py`, `Main.sin`, `utils/download_googledex.py`.
6. **`UCOScheduler` class.** Write the class, decompose `get_next` into the
   stages above. Keep a module-level `get_next(ctime, seeing, slowdown, targets,
   **kw)` that constructs a scheduler and delegates, *temporarily*, so
   `sim_night.py` still runs unchanged and can be diffed against step 4.
7. **Port the call sites.** `Observe.py`, `Main.sin`, `sim_night.py`,
   `sim_nights.py` to the object API. Delete the temporary module-level shim.
8. **Bug fixes**, one commit each, each one re-run through `sim_nights.py` and
   diffed against the step-4 baseline. Items 2, 3, 4, 5 first (clear fixes),
   then item 1 (cal stars) last, since it is the one that visibly changes target
   selection.

## How this gets verified

There is no test suite, so verification is differential.

* **Baseline**: run `utils/sim_nights.py` on the current deployed code over a
  fixed set of nights with a pinned `googledex.dat`, `rank_table`,
  `hour_table` and `observed_targets`, and save the full sequence of selected
  scriptobs lines.
* **After steps 1-7** (pure refactor), the same run must produce a
  byte-identical sequence. Any difference is a refactor bug, not an
  improvement.
* **After step 8**, differences are expected. Each one is attributed to a
  specific fix before the commit lands.
* `test_UCOScheduler.py` runs against a real rank table in `--test` mode as a
  smoke test that the Google-sheet path and the backup path both still work.

The pinned inputs live in a fixtures directory so the comparison is repeatable;
they are not generated fresh per run.

## Open questions

1. `TargetTables` vs `UCOTargetTables` for the class and module name.
2. Should `get_next` keep returning a plain `dict`, or become a small
   `Target` result class? The dict is consumed in `Observe.py` by string key
   (`self.target['NAME']`, `self.target["SCRIPTOBS"]`) and in the sim scripts.
   A result class is nicer but widens the diff; the plan above keeps the dict.
3. `Observe` mutates `self.do_temp` and `self.n_temps` during the night to cap
   template observations at `tot_temps`. That template budget arguably belongs
   to the scheduler. The plan leaves it in `Observe` for now and passes
   `do_templates` per call.
