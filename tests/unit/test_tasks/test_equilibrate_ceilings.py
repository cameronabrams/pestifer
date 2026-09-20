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

    def test_the_asymmetric_membrane_pre_embed_stages_clear_800000(self):
        # both pre-embed stages (calibration patch, quilt) hit 800000 unconverged in the sweep;
        # the post-embed stage converged there and is left alone
        cfg = self._example('17', 'hiv-mpertm3-membrane2.yaml')
        ceilings = self._ceilings(cfg, 'membrane_equilibrate')
        self.assertEqual(sum(1 for c in ceilings if c > 800000), 2, ceilings)
