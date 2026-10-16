# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Static, pre-execution validation of a task pipeline.

Each task declares a :class:`TaskContract` -- what pipeline "currencies" it
**requires** on entry and what it **provides** on exit -- via
:meth:`BaseTask.pipeline_contract`.  :func:`validate_pipeline` walks a task list
once, before any work is done, and reports malformed hand-offs (a task consuming
a currency no earlier task produced, an origin task that would silently discard an
already-built system, or anything scheduled after the terminal task) with an
explanation the user can act on.

The currencies are the things that actually flow between tasks:

- ``SOURCE`` (``base_coordinates``): raw fetched coordinates -- produced by ``fetch``,
  consumed by ``psfgen``.
- ``STATE`` (``state``): the built PSF/PDB[/XSC] fileset -- produced by ``psfgen``,
  ``continuation``, ``merge``, and every transform; consumed by transforms and
  ``terminate``.
- ``MOLECULE`` (``base_molecule``): the in-memory :class:`Molecule` (sequence,
  biological assembly, chains, loops) -- produced only by the ``psfgen`` family;
  consumed by ``cleave``/``ligate``/``pdb2pqr``.  ``continuation``/``merge`` give
  files but *not* this object.
- ``MD_OUTPUT`` (``md_output``): molecular-dynamics timeseries -- produced by ``md``,
  consumed by ``mdplot``.
- ``SOLVATED`` (``solvated``): the system carries bulk solvent -- a soft flag used
  only to warn about re-solvation.
"""
import logging

from ..core.errors import PestiferBuildError

logger = logging.getLogger(__name__)

# pipeline currencies
SOURCE = 'base_coordinates'
STATE = 'state'
MOLECULE = 'base_molecule'
MD_OUTPUT = 'md_output'
SOLVATED = 'solvated'


class TaskContract:
    """A task's pipeline input/output contract.

    Parameters
    ----------
    requires : iterable of str
        Currencies that must already be available when the task runs (hard).
    provides : iterable of str
        Currencies available after the task runs.
    discards_state : bool
        True if the task builds a fresh system from scratch; running it when a
        ``STATE`` already exists silently discards that prior system.
    terminal : bool
        True if the task closes the build (``terminate``); nothing may follow it.
    warn_if_present : iterable of str
        Currencies whose prior presence is suspicious (a warning, not an error).
    standalone : bool
        True if the task is a self-contained utility that operates on explicit inputs and
        does not participate in the build data-flow (e.g. ``desolvate``); it is an error to
        place it inside a build pipeline.
    """

    def __init__(self, requires=(STATE,), provides=(STATE,), discards_state=False,
                 terminal=False, warn_if_present=(), standalone=False):
        self.requires = frozenset(requires)
        self.provides = frozenset(provides)
        self.discards_state = discards_state
        self.terminal = terminal
        self.warn_if_present = frozenset(warn_if_present)
        self.standalone = standalone


# per-currency explanation used when a required currency is missing
_HINT = {
    STATE: ("requires an existing molecular system, but no earlier task establishes one -- "
            "begin the pipeline with `fetch` + `psfgen`, `continuation`, or `merge`"),
    SOURCE: ("builds a PSF from fetched coordinates (`psfgen`), but no `fetch` task precedes it -- "
             "add a `fetch`, or, if the system is already built (e.g. by `continuation`), drop the `psfgen` "
             "(it would rebuild the system from scratch)"),
    MOLECULE: ("needs the in-memory molecule that `psfgen` builds (chains, sequence, loops); "
               "`continuation`/`merge` provide files but not that object, so this task must follow a "
               "`psfgen`-family task, not a `continuation`"),
    MD_OUTPUT: "needs molecular-dynamics output, but no `md` task precedes it",
}


def validate_pipeline(tasks):
    """Statically check a task list for malformed hand-offs before execution.

    Raises :class:`~pestifer.core.errors.PestiferBuildError` listing every hard
    problem found; emits warnings for softer, occasionally-legitimate issues.
    ``tasks`` is any sequence of objects exposing ``index``, ``taskname``,
    ``specs``, and ``pipeline_contract(specs)``.
    """
    import inspect
    available = set()
    errors, warnings = [], []
    terminal_task = None
    for task in tasks:
        # Most contracts depend only on the task's own specs, but a pipeline-aware contract may also
        # need to know what currencies precede it (e.g. psfgen preserves an incoming STATE vs. builds
        # from a fetched SOURCE).  Pass `available` only to contracts that declare the parameter, so
        # the common `pipeline_contract(cls, specs)` signature is untouched.
        try:
            _accepts_available = 'available' in inspect.signature(task.pipeline_contract).parameters
        except (TypeError, ValueError):
            _accepts_available = False
        if _accepts_available:
            contract = task.pipeline_contract(task.specs, available=frozenset(available))
        else:
            contract = task.pipeline_contract(task.specs)
        label = f"task {task.index:02d} '{task.taskname}'"
        if contract.standalone:
            errors.append(f"{label} is a standalone utility that operates on explicit file inputs "
                          f"(run it on its own, e.g. `pestifer {task.taskname} ...`); it does not "
                          f"participate in the build data-flow and cannot be placed in a `tasks:` list")
            continue  # contributes nothing to the pipeline
        if terminal_task is not None:
            errors.append(f"{label} is scheduled after the terminal task '{terminal_task}', which "
                          f"packages and closes the build -- nothing may run after it")
        missing = contract.requires - available
        for currency in sorted(missing):
            errors.append(f"{label} {_HINT.get(currency, f'requires {currency}, which no earlier task provides')}")
        # only meaningful if the task could otherwise run (its inputs are satisfied)
        if not missing and contract.discards_state and STATE in available:
            errors.append(f"{label} builds a new molecular system from scratch and would silently discard "
                          f"the system an earlier task already established -- remove it, or run it in a "
                          f"separate pipeline if a rebuild is intended")
        for currency in sorted(contract.warn_if_present & available):
            if currency == SOLVATED:
                warnings.append(f"{label} acts on an already-solvated system; its solvation may duplicate "
                                f"or conflict with the existing bulk water")
        available |= contract.provides
        if contract.terminal:
            terminal_task = task.taskname
    errors.extend(_check_unminimized_dynamics(tasks))

    for w in warnings:
        logger.warning(f'pipeline check: {w}')
    if errors:
        bullets = '\n  - '.join(errors)
        raise PestiferBuildError(f'Malformed task pipeline ({len(errors)} problem(s) found before execution):'
                                 f'\n  - {bullets}')


#: Tasks that MOVE ATOMS without integrating equations of motion, so whatever they displace is
#: left unrelaxed.
#:
#: Only ``ring_check`` has an observed failure behind it: it rotates a pendant branch or an
#: aromatic side chain to unthread a piercing, and example 17 ran 410 K dynamics directly on its
#: output for several releases.  ``cleave``/``ligate`` are included because the same reasoning
#: applies -- they cut and reconnect a chain -- and because the rule is free there: every shipped
#: use of both already minimizes immediately afterwards (example 15), so nothing is being
#: retro-fitted to a working pipeline on a hunch.
#:
#: The list is deliberately short and explicit.  A task added later is simply not covered, which
#: is a visible gap rather than a silent misclassification; `solvate`, `psfgen` and
#: `make_membrane_system` also hand over unminimized coordinates, but every shipped pipeline
#: minimizes after them and no failure has been traced to them, so they are left out rather than
#: guessed in.
_MOVES_ATOMS = frozenset({'ring_check', 'cleave', 'ligate'})


def _is_dynamics(task):
    """True for a task that integrates equations of motion (not a minimization)."""
    name = getattr(task, 'taskname', '')
    if name == 'membrane_equilibrate':
        return True
    if name != 'md':
        return False
    return str((task.specs or {}).get('ensemble', '') or '').casefold() != 'minimize'


def _is_minimize(task):
    return (getattr(task, 'taskname', '') == 'md'
            and str((task.specs or {}).get('ensemble', '') or '').casefold() == 'minimize')


def _check_unminimized_dynamics(tasks):
    """Refuse dynamics on coordinates an earlier task moved and nothing minimized.

    This is the static half of a rule ``make_membrane_system._guard_piercings`` already
    enforces dynamically for the membrane relaxation protocol it generates, and whose docstring
    states the consequence: a resolution "frees the other lipid into a strained pose", and a
    threading "survives further minimization (it is topological, not just strained) and then
    RATTLE-fails the first dynamics step".  That guard never saw the user's own top-level task
    list, so example 17 shipped the same sequence with the minimize on the WRONG SIDE of its
    ``ring_check`` -- minimize, ring_check, then 410 K dynamics -- for several releases.

    A RATTLE failure costs the whole build and is reported by NAMD as an atom serial, hours in.
    This costs a second, before anything runs, and names the fix.
    """
    problems = []
    moved_by = None
    for task in tasks:
        name = getattr(task, 'taskname', '')
        if name in _MOVES_ATOMS:
            moved_by = task
            continue
        if _is_minimize(task):
            moved_by = None
            continue
        if _is_dynamics(task) and moved_by is not None:
            problems.append(
                f"task {task.index:02d} '{name}' integrates dynamics on coordinates that task "
                f"{moved_by.index:02d} '{moved_by.taskname}' moved, with no minimization in "
                f"between -- `{moved_by.taskname}` displaces atoms to resolve a problem and "
                f"leaves them unrelaxed, which RATTLE-fails on the first dynamics step. "
                f"Insert `- md: {{ensemble: minimize, minimize: 1000}}` between them")
            moved_by = None     # report each moving task once, not once per following stage
    return problems
