# Replacing the scheduler's result dict with a `Target` class

Status: plan, not started. Three decisions below are still open.
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
| `UCOScheduler._add_template` (`:562-574`) | **changes it** | reads `owner`; replaces `SCRIPTOBS`; sets `isTemp`, `DECKER` |
| `UCOScheduler.get_next` (`:357-366`) | **changes it** | reads `isTemp` (template budget); sets `template_conditions_met`; stores it on `self.result` |
| `Observe.get_target` (`Observe.py:498-557`) | reads it | `NAME`, `VMAG`, `BV`, `DECKER`, `mode`, `PRI`, `COUNTS`, `EXP_TIME`, `NEXP`, `isTemp`, `isTOO`, `template_conditions_met`; **pops** `SCRIPTOBS` |
| `Observe` fixed-list path (`:675-689`, `:856-859`) | **builds its own dicts** | `self.fixed_target = {'SCRIPTOBS': [...]}`, then `self.target = {'SCRIPTOBS': ...}`, holding no target data at all |
| `pop_next`, `empty_queue` (`:389-414`) | treat both as a queue | `'SCRIPTOBS' in self.target.keys()`, then pop |
| `NightSim.compute_simulation` (`utils/NightSim.py:103-133`) | reads it | `VMAG`, `BV`, `DECKER`, `COUNTS`, `EXP_TIME`, `NAME` |
| `utils/sim_night.py`, `utils/sim_nights.py` | read it | `isBstar`, `isTOO`, `NAME`, `NEXP`, `owner`; pop `SCRIPTOBS` |
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
would be awkward. This is decision 1 below.

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
* Values keep their numpy types in this change (decision 3).

## Impact, per file

| File | Change |
| --- | --- |
| `Main/Target.py` | new |
| `Main/UCOScheduler.py` | `make_result` returns a `Target`; `_add_template` calls `set_template_lines`; `get_next` uses `is_temp` and sets `template_conditions_met` as attributes |
| `Main/Observe.py` | `get_target` uses attributes. With decision 1 as recommended: a new `self.pending_lines` list; `pop_next`, `empty_queue` and the two fixed-list sites use it, and the `{'SCRIPTOBS': ...}` dicts go away. This is the only part that is not mechanical |
| `utils/NightSim.py` | six key reads become attributes |
| `utils/sim_night.py`, `utils/sim_nights.py` | about five reads each |
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
2. **Move the callers to attributes,** one commit each: the scheduler
   internals, then `NightSim` and the sim scripts, then the tests.
3. **Split Observe's queue** (if decision 1 goes that way). This is the
   riskiest step and, like the rest of `Observe`, cannot be run without `ktl`;
   review it by reading, then watch the first night's log.
4. **Remove the temporary dict interface.** Any caller still using a key then
   fails loudly instead of silently working.

## Decisions to make

### 1. Split Observe's queue from the target?

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

`RA`, `DEC`, `PM_RA`, `PM_DEC`, `I2` and `BINNING` are set by `make_result` but
never read by a caller.

* **Keep them (recommended).** They cost nothing, `ra`/`dec` are useful in
  logs, and `to_dict()` then reproduces today's dict exactly, which keeps the
  whole-result prints identical.
* **Drop them.** A smaller class, but the printed results change, so the
  tests' output can no longer be compared byte-for-byte across that commit.

### 3. Convert numpy scalars to plain Python types?

Values are numpy scalars today (`np.float64(8.43)`, `np.str_('HD66485')`).

* **Not in this change (recommended).** Converting changes how results print,
  which would break the byte-identical comparison in the same commits that
  restructure the code. If wanted, do it as its own commit after step 4, where
  the only expected difference is formatting.
* **Convert in `from_star_table`.** Cleaner values from the start, at the cost
  of losing the identical-output check for step 1.
