# Replacing the scheduler's result dict with a `Target` class

Status: all four steps done. Steps 1 and 2 are committed (`057172f`,
`314cfc1`); steps 3 and 4 are uncommitted. The three decisions are
made (see "Decisions").
Branch: `scheduler_object`, after the `UCOScheduler` class refactor
(see `docs/scheduler_refactor.md`).

## Goal

`UCOScheduler.get_next()` returns a plain `dict` with 21 string keys
(`'NAME'`, `'SCRIPTOBS'`, `'isTemp'`, ...). Callers read and change it by key,
so a misspelled key fails only at run time, and the dict has no single
definition anywhere. The goal is to return a small `Target` class instead,
without changing which targets are picked or what is written to disk.

## Who builds and uses the dict today

| Where | Role | Keys / operations |
| --- | --- | --- |
| `make_result` (`UCOScheduler.py:139-186`) | **builds it** | all 21 keys; values are numpy scalars (`np.float64`, `np.str_`) |
| `UCOScheduler._add_template` (`:567-579`) | **changes it** | reads `owner`; replaces `SCRIPTOBS`; sets `isTemp`, `DECKER` |
| `UCOScheduler.get_next` (`:360-371`) | **changes it** | reads `isTemp` (template budget) and `isTOO` (ToO tracking); sets `template_conditions_met`; stores it on `self.result` |
| `Observe.get_target` (`Observe.py:497-554`) | reads it | `NAME`, `VMAG`, `BV`, `DECKER`, `mode`, `PRI`, `COUNTS`, `EXP_TIME`, `NEXP`, `isTemp`, `isTOO` (to write `MASTER_OBSTOO`), `template_conditions_met`; **pops** `SCRIPTOBS` |
| `Observe` fixed-list path (`:672-686`, `:855-856`) | **builds its own dicts** | `self.fixed_target = {'SCRIPTOBS': [...]}`, then `self.target = {'SCRIPTOBS': ...}`, holding no target data at all |
| `pop_next`, `empty_queue` (`:388-413`) | treat both as a queue | `'SCRIPTOBS' in self.target.keys()`, then pop |
| `NightSim.compute_simulation` (`utils/NightSim.py:103-133`) | reads it | `VMAG`, `BV`, `DECKER`, `COUNTS`, `EXP_TIME`, `NAME` |
| `utils/sim_night.py`, `utils/sim_nights.py` | read it | `isBstar`, `NAME`, `NEXP`, `owner`; pop `SCRIPTOBS` |
| `Main/test_UCOScheduler.py` | reads it | pops `SCRIPTOBS`, reads `TOTEXP_TIME`, and **prints the whole dict** at lines 81 and 86 |

Six keys are never read by any caller: `RA`, `DEC`, `PM_RA`, `PM_DEC`, `I2`,
`BINNING`.

## The design problem: `SCRIPTOBS` does two jobs

`SCRIPTOBS` is both:

* **A result:** the scriptobs lines the scheduler chose for one target
  (one line normally, two for a B star, three for a template sequence).
* **Observe's work queue:** `Observe` pops lines off it as it sends them to
  scriptobs, and the fixed-list path puts a whole starlist into the same slot
  with no target attached.

A `Target` class that also has to act as "a list of lines with no target"
would be awkward. Decision 1 below splits the two.

## Proposed class

A plain class in a new `Main/Target.py`, matching the repo, which does not use
`dataclasses` anywhere. Attributes are snake_case versions of today's keys:

```python
class Target(object):
    def __init__(self, name, owner, ra, dec, pm_ra, pm_dec, vmag, bv, counts,
                 exp_time, nexp, totexp_time, pri, decker, i2, binning,
                 is_bstar=False, is_too=False, mode=''):
        ...
        self.is_temp = False
        self.scriptobs = []                  # lines, last one sent first (pop order)
        self.template_conditions_met = False

    @classmethod
    def from_star_table(cls, star_table, idx, star, totexptime, pri, bstar)

    def set_template_lines(self, lines, decker)   # the _add_template bookkeeping
    def to_dict(self)                             # the old dict, same keys and values
    def __repr__(self)                            # repr(self.to_dict())
```

| Old key | Attribute |
| --- | --- |
| `NAME`, `owner` | `name`, `owner` |
| `RA`, `DEC`, `PM_RA`, `PM_DEC` | `ra`, `dec`, `pm_ra`, `pm_dec` |
| `VMAG`, `BV` | `vmag`, `bv` |
| `COUNTS`, `EXP_TIME`, `NEXP`, `TOTEXP_TIME` | `counts`, `exp_time`, `nexp`, `totexp_time` |
| `PRI` | `pri` |
| `DECKER`, `I2`, `BINNING`, `mode` | `decker`, `i2`, `binning`, `mode` |
| `isTemp`, `isBstar`, `isTOO` | `is_temp`, `is_bstar`, `is_too` |
| `SCRIPTOBS` | `scriptobs` |
| `template_conditions_met` | `template_conditions_met` |

* `make_result` stays a module-level function in `UCOScheduler.py` and returns
  a `Target`. It builds the object with `Target.from_star_table` and then adds
  the scriptobs lines, which depend on `ScriptobsLine` and on `bstar`/`focval`.
* `__repr__` returns `repr(self.to_dict())`, so log lines and the two
  whole-result prints in `test_UCOScheduler.py` stay byte-identical. That keeps
  the differential checks used throughout the refactor working.
* Values keep their numpy types (decision 3).

## Impact, per file

| File | Change |
| --- | --- |
| `Main/Target.py` | new |
| `Main/UCOScheduler.py` | `make_result` returns a `Target`; `_add_template` calls `set_template_lines`; `get_next` uses `is_temp` and `is_too` and sets `template_conditions_met` as attributes |
| `Main/Observe.py` | `get_target` uses attributes. Per decision 1: a new `self.pending_lines` list; `pop_next`, `empty_queue` and the two fixed-list sites use it, and the `{'SCRIPTOBS': ...}` dicts go away. This is the only part that is not mechanical |
| `utils/NightSim.py` | six key reads become attributes |
| `utils/sim_night.py`, `utils/sim_nights.py` | about four reads each |
| `Main/test_UCOScheduler.py` | a few reads; `print(result)` output is unchanged thanks to `__repr__` |
| `Main/ScriptobsLine.py`, `Main/Observability.py`, `Main/UCOTargetTables.py`, `Main/Main.sin` | no change |

## Steps

Each step should leave the seeded `sim_night.py` / `sim_nights.py` runs and
`test_UCOScheduler.py` identical to the step before, checked the same way as
the class refactor (old and new side by side, on fresh copies of the
`test_main` fixtures).

1. **Add `Target` with a temporary dict interface.** `__getitem__`,
   `__setitem__`, `keys()` and `__contains__` map the old keys onto the
   attributes. `make_result` returns a `Target`. No caller changes.
   **Done:** `Main/Target.py`; `make_result` builds the `Target` with
   `Target.from_star_table` and still fills the scriptobs lines through the
   old keys. Verified against `14a3510`: `test_UCOScheduler.py` (which prints
   whole results), seeded `sim_night.py` on 2026-09-28 (B star and a ToO) and
   on 2026-10-10 with every target `Template=N` (two templates), and seeded
   `sim_nights.py` 2026/09/28-10/01 all identical. Old and new `make_result`
   on the same inputs (three rows, with and without `bstar`, one row a ToO)
   gave the same keys in the same order, with identical values and types.
2. **Move the callers to attributes,** one commit each: the scheduler
   internals, then `NightSim` and the sim scripts, then the tests.
   **Done**, in three file groups that can be committed separately:
   (a) `UCOScheduler.py` (`make_result`, `get_next`, `_add_template`) and
   `Target.set_template_lines`; (b) `utils/NightSim.py`, `utils/sim_night.py`,
   `utils/sim_nights.py`; (c) `Main/test_UCOScheduler.py`. Verified against
   `14a3510` with the same four runs as step 1, all identical, using a copy of
   the tree in which the temporary dict interface raises an error: none of
   these files touches it any more. `Observe.py` is now its only user.
3. **Split Observe's queue** (decision 1). This is the
   riskiest step and, like the rest of `Observe`, cannot be run without `ktl`;
   review it by reading, then watch the first night's log.
   **Done:** `self.fixed_target` became `self.fixed_lines` (a list or
   `None`), and the old `self.target['SCRIPTOBS']` queue became
   `self.pending_lines`. `self.target` is now only ever a `Target` or `None`.
   After `get_next`, `pending_lines` takes a *copy* of `target.scriptobs`, so
   the `Target` keeps a full record of what was chosen. Every place that used
   to empty the queue by setting `self.target = None` (an empty starlist, a
   finished fixed list, `get_next` returning `None`) now also empties
   `pending_lines`, so the same lines are dropped at the same moments. Running
   a fixed list now sets `self.target = None` and moves the list into
   `pending_lines`, instead of building a dict with only `SCRIPTOBS`.
   Checked by loading `Observe.py` with stand-ins for `ktl`, `APFTask` and the
   hardware modules (it imports and the class builds), by confirming it uses
   no undefined names, and that every `self.target.<attr>` it reads exists on
   `Target`. The logic itself was checked by reading only.
4. **Remove the temporary dict interface.** Any caller still using a key then
   fails loudly instead of silently working.
   **Done:** `__getitem__`, `__setitem__`, `__contains__`, `keys()` and the
   key-to-attribute map they used were removed; `to_dict()` and `__repr__`
   stay. A search of `Main/` and `utils/` finds no remaining old-key access on
   a result. `test_UCOScheduler.py`, seeded `sim_night.py` (2026-09-28, and
   2026-10-10 with every target `Template=N`) and seeded `sim_nights.py`
   remain identical to `14a3510`, the code from before `Target`.

## Decisions

### 1. Split Observe's queue from the target?

**Decided: split.**

* **Split (recommended).** `Target.scriptobs` holds the lines chosen for that
  target. `Observe` keeps its own `self.pending_lines`, filled from
  `target.scriptobs` or from the fixed list; `pop_next` and `empty_queue` work
  on that list. `Target` then always means a real target, and the fixed-list
  code stops pretending to be one.
* **Keep them together.** `Observe` keeps popping `target.scriptobs`, and the
  fixed-list path builds a `Target` with only `scriptobs` set. Smaller change
  to `Observe`, but `Target` needs every field to be optional, and the
  awkwardness stays.

### 2. Drop the six fields nobody reads?

**Decided: keep them.**

`RA`, `DEC`, `PM_RA`, `PM_DEC`, `I2` and `BINNING` are set by `make_result` but
never read by a caller.

* **Keep them (recommended).** They cost nothing, `ra`/`dec` are useful in
  logs, and `to_dict()` then reproduces today's dict exactly, which keeps the
  whole-result prints identical.
* **Drop them.** A smaller class, but the printed results change, so the
  tests' output can no longer be compared byte-for-byte across that commit.

### 3. Convert numpy scalars to plain Python types?

**Decided: not in this change.**

Values are numpy scalars today (`np.float64(8.43)`, `np.str_('HD66485')`).

* **Not in this change (recommended).** Converting changes how results print,
  which would break the byte-identical comparison in the same commits that
  restructure the code. If wanted, do it as its own commit after step 4, where
  the only expected difference is formatting.
* **Convert in `from_star_table`.** Cleaner values from the start, at the cost
  of losing the identical-output check for step 1.
