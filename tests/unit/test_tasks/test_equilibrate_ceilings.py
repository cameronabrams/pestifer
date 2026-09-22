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

    def _ceilings(self, node, taskname, found=None):
        found = [] if found is None else found
        if isinstance(node, dict):
            for k, v in node.items():
                if k == taskname and isinstance(v, dict) and 'max_steps' in v:
                    found.append(v['max_steps'])
                self._ceilings(v, taskname, found)
        elif isinstance(node, list):
            for v in node:
                self._ceilings(v, taskname, found)
        return found

    def test_acetone_clears_its_slowest_converging_replica(self):
        # acetone decorrelates slowly, so its precision gate needs sampling, not just steps: one
        # replica converged at 227000 and another stopped at the old 250000 ceiling
        cfg = self._example('24', 'subtilisin-acetone.yaml')
        ceilings = self._ceilings(cfg, 'density_equilibrate')
        self.assertEqual(len(ceilings), 1)
        self.assertGreater(ceilings[0], 2 * SLOWEST_OBSERVED_CONVERGENCE)

    def test_the_asymmetric_membrane_pre_embed_stages_clear_the_slowest_replica(self):
        """Both pre-embed stages must have room for ex17's slow replica, not its median one.

        3.22.1: both hit 800000.  3.22.8 raised them to 1,200,000 (patch) and 1,500,000 (quilt) and
        the prediction inverted -- the quilt, thought furthest from settling, converged at
        1,455,700, while the patch, thought nearly fixed, hit its new ceiling.  Across replicas that
        patch spans 236,870 / 372,360 / >1,200,000, so a budget sized on the median is a budget that
        fails one build in three.
        """
        cfg = self._example('17', 'hiv-mpertm3-membrane2.yaml')
        ceilings = sorted(self._ceilings(cfg, 'membrane_equilibrate'))
        pre_embed = ceilings[-2:]        # the two pre-embed stages; the post-embed one converges low
        for c in pre_embed:
            self.assertGreater(c, 1.5 * SLOWEST_OBSERVED_MEMBRANE_CONVERGENCE, ceilings)
        self.assertEqual(pre_embed[0], pre_embed[1],
                         'the patch and quilt budgets are matched so neither is the one that runs '
                         f'out first: {ceilings}')
