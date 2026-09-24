# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""The step ceilings of the self-terminating equilibrations, against what real sweeps measured.

A ceiling is a judgment call, not a constant: it binds only on a run that has *not* converged, so
setting it too low does not save time on a healthy build -- it converts a slow replica into a
system reported as possibly unsettled.  The 3.22.1 example sweep (93 builds, 111 adaptive stages)
measured how far real runs go, and these tests pin the ceilings to stay clear of it.  Raise them
again rather than relaxing a tolerance if a later sweep finds another.
"""
import os
import unittest

import yaml

import pestifer

SCHEMA = os.path.join(os.path.dirname(pestifer.__file__), 'schema', 'base.yaml')
EXAMPLES = os.path.join(os.path.dirname(pestifer.__file__), 'resources', 'examples')

#: Longest converging density_equilibrate run in the 3.22.1 sweep, in steps (ex24/rep-03, acetone).
SLOWEST_OBSERVED_CONVERGENCE = 227000

#: Longest converging pre-embed membrane_equilibrate stage, in steps: ex17/rep-01's quilt in the
#: 3.22.8 sweep, which finished at 1,455,700 of the 1,500,000 it had been given -- 3% headroom.
#: Its sibling replicas took 704,400 and 600,080, so the budget has to cover the outlier.
SLOWEST_OBSERVED_MEMBRANE_CONVERGENCE = 1455700

#: ex17 has THREE pre-embed stages from TWO budgets: the `patch` protocol is applied to both
#: calibration patches (patchA and patchB), and the `quilt` protocol once.  patchB is not spare
#: capacity -- rep-01's converged at 1,357,350, also past the 1,200,000 ceiling patchA died on --
#: so raising the `patch` budget has to cover two stages that both run long, not one.
SLOWEST_OBSERVED_PATCHB_CONVERGENCE = 1357350

#: DO NOT LOWER THESE CEILINGS ON THE STRENGTH OF A FAST SWEEP.  A probe of ex17 at 3.23.1
#: (2026-09-24, 3 replicas) converged patchA at 196,870 / 252,360 / 236,870 and patchB at
#: 348,220 / 342,350 / 393,220 -- five to seven times below the numbers above.  That is not
#: evidence the budgets are oversized.  It is evidence the SLOW MODE IS RARE, which is the
#: same thing these ceilings exist for.
#:
#: The slow mode is not a property of a seed.  Same seed (27021972), same NAMD build, same
#: config: 3.22.8's rep-01 hit the 1,200,000 ceiling and the probe's rep-01 converged at
#: 196,870, a 6x difference.  Multicore NAMD is not bitwise reproducible, and the adaptive
#: NPgT chunk length is rounded to 10 steps -- usually the rounding absorbs the FP noise
#: (rep-03 reproduced 500/330/490/730 exactly and converged to the same step both times), but
#: rep-01 crossed a boundary at chunk 2 (320 vs 330), inside the first 500 steps, and never
#: rejoined.  So a replica is a draw from a distribution, not a re-run.
#:
#: Tally so far: the tail appeared in 1 of 6 observed ex17 patchA runs.  n=3 cannot resolve a
#: 1-in-6 event either way, so "2,500,000 clears the slow mode" is CONSISTENT WITH the data and
#: NOT ESTABLISHED BY IT.  Settling it needs many more replicas or a way to induce the slow mode
#: directly.  Until then the budget is sized for the tail that was observed, not the median that
#: keeps being re-observed.
SLOW_MODE_OBSERVED_IN_N_OF_M = (1, 6)


def _find(node, name, type_='dict'):
    """The first schema node with this ``name`` (and ``type``), depth first."""
    if isinstance(node, dict):
        if node.get('name') == name and node.get('type') == type_:
            return node
        for v in node.values():
            found = _find(v, name, type_)
            if found:
                return found
    elif isinstance(node, list):
        for v in node:
            found = _find(v, name, type_)
            if found:
                return found
    return None


def _attr_default(task_node, attr):
    for a in task_node.get('attributes', []):
        if a.get('name') == attr:
            return a.get('default')
    raise AssertionError(f'{task_node.get("name")} has no {attr} attribute')


class TestDensityEquilibrateCeiling(unittest.TestCase):

    def test_the_default_ceiling_clears_the_slowest_converging_run_seen(self):
        """Three builds of the 3.22.1 sweep stopped at the old 100000 default with their sibling
        replicas converging in 53000-97000 steps: the margin over a typical run was too thin to
        cover a slow one."""
        schema = yaml.safe_load(open(SCHEMA))
        default = _attr_default(_find(schema, 'density_equilibrate'), 'max_steps')
        self.assertGreaterEqual(default, 200000)


class TestExampleCeilings(unittest.TestCase):
    """Two examples set their own ceilings, and both were hit in the 3.22.1 sweep."""

    def _example(self, number, name):
        path = os.path.join(EXAMPLES, number, 'inputs', name)
        return yaml.safe_load(open(path))

    def _ceilings(self, node, taskname, found=None, path=''):
        """Every ``max_steps`` under ``taskname``, each with the config path it sits at.

        The path is what separates a pre-embed budget from a post-embed one; selecting them by
        sorted value happened to work only while the pre-embed numbers were the larger pair, and
        would have silently picked the wrong stages the moment that stopped being true.
        """
        found = [] if found is None else found
        if isinstance(node, dict):
            for k, v in node.items():
                here = f'{path}/{k}'
                if k == taskname and isinstance(v, dict) and 'max_steps' in v:
                    found.append((here, v['max_steps']))
                self._ceilings(v, taskname, found, here)
        elif isinstance(node, list):
            for i, v in enumerate(node):
                self._ceilings(v, taskname, found, f'{path}[{i}]')
        return found

    @staticmethod
    def _pre_embed(ceilings):
        """The budgets that govern bilayer construction, before the protein is embedded.

        They live under ``make_membrane_system/bilayer/relaxation_protocols``; the post-embed
        stages are top-level tasks.
        """
        return [(p, v) for p, v in ceilings if 'relaxation_protocols' in p]

    def test_acetone_clears_its_slowest_converging_replica(self):
        # acetone decorrelates slowly, so its precision gate needs sampling, not just steps: one
        # replica converged at 227000 and another stopped at the old 250000 ceiling
        cfg = self._example('24', 'subtilisin-acetone.yaml')
        ceilings = self._ceilings(cfg, 'density_equilibrate')
        self.assertEqual(len(ceilings), 1)
        self.assertGreater(ceilings[0][1], 2 * SLOWEST_OBSERVED_CONVERGENCE)

    def test_the_asymmetric_membrane_pre_embed_stages_clear_the_slowest_replica(self):
        """Both pre-embed budgets must have room for ex17's slow replica, not its median one.

        3.22.1: both hit 800000.  3.22.8 raised them to 1,200,000 (patch) and 1,500,000 (quilt) and
        the prediction inverted -- the quilt, thought furthest from settling, converged at
        1,455,700, while the patch, thought nearly fixed, hit its new ceiling.  Across replicas that
        patch spans 236,870 / 372,360 / >1,200,000, so a budget sized on the median is a budget that
        fails one build in three.

        The `patch` budget governs TWO stages, not one: patchA and patchB both run the `patch`
        protocol, and rep-01's patchB converged at 1,357,350 -- itself past the ceiling patchA died
        on.

        The 3.23.1 probe (2026-09-24) did not settle whether 2,500,000 clears the slow mode: all
        three replicas converged fast and none entered it.  See SLOW_MODE_OBSERVED_IN_N_OF_M -- the
        budget stays sized for the tail, and a fast sweep is not a reason to trim it.
        """
        cfg = self._example('17', 'hiv-mpertm3-membrane2.yaml')
        ceilings = self._ceilings(cfg, 'membrane_equilibrate')
        pre_embed = self._pre_embed(ceilings)
        self.assertEqual(len(pre_embed), 2,
                         f'expected the patch and quilt budgets, got {pre_embed}')
        for path, c in pre_embed:
            self.assertGreater(c, 1.5 * SLOWEST_OBSERVED_MEMBRANE_CONVERGENCE, path)
            self.assertGreater(c, 1.5 * SLOWEST_OBSERVED_PATCHB_CONVERGENCE, path)
        self.assertEqual(pre_embed[0][1], pre_embed[1][1],
                         'the patch and quilt budgets are matched so neither is the one that runs '
                         f'out first: {pre_embed}')

    def test_the_pre_embed_budgets_are_not_merely_the_largest_two(self):
        """The selection above must key on where a budget sits, not on how big it is.

        Guards the test itself: if `_pre_embed` fell back to sorting by value, a post-embed budget
        raised above the pre-embed ones would be silently tested in their place.
        """
        cfg = self._example('17', 'hiv-mpertm3-membrane2.yaml')
        ceilings = self._ceilings(cfg, 'membrane_equilibrate')
        self.assertGreater(len(ceilings), 2, 'ex17 should also have post-embed stages')
        inflated = [(p, 10 ** 9 if 'relaxation_protocols' not in p else v) for p, v in ceilings]
        self.assertEqual([v for _, v in self._pre_embed(inflated)],
                         [v for _, v in self._pre_embed(ceilings)])
