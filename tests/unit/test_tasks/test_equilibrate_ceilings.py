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

#: READ A STAGE'S OWN STEPS, NOT THE NUMBER IN ITS LOG LINE.  NAMD's step counter carries across
#: the stages of one make_membrane_system: patchB's first chunk starts at the step patchA ended on,
#: and the quilt's at patchB's (e.g. "patchB: NPgT chunk 1 ... (firsttimestep 2512000)" right after
#: patchA's CEILING at 2,500,000).  ``max_steps`` counts the stage's OWN steps, so "CONVERGED at
#: step N" for patchB or the quilt is not comparable to a ceiling.  Until 2026-11-09 the figures
#: here were read straight off those lines: 3.22.8 rep-01's "patchB converged at 1,357,350, past
#: the ceiling patchA died on" was about 150k steps of its own, and its "quilt at 1,455,700 of
#: 1,500,000, 3% headroom" about 95k.  Neither stage was ever near its budget.
#:
#: Own steps, ex17 on afd96d1d, three seeds (2026-11-09; firsttimestep subtracted from the
#: CONVERGED step, each read from the build log):
#:   patchB  381,040 / 314,520 / 164,520
#:   quilt   497,220 / 489,190 / 670,040
SLOWEST_OBSERVED_PATCHB_CONVERGENCE = 381040
SLOWEST_OBSERVED_QUILT_CONVERGENCE = 670040

#: patchA (the upper, Lo raft patch) is the one stage that does not settle.  Its counter starts
#: fresh (firsttimestep 6000), so its numbers were always its own: it converged at 337k tensionless
#: steps in one seed and hit the 2,500,000 ceiling in the other two, its APL still falling 0.1-0.2%
#: per 100k.  Its `patch` budget is therefore a deliberate CAP (1,200,000), not a budget to converge
#: in -- which is safe only because the leaflet counts come from `bilayer.leaflet_ratio` rather
#: than from this patch's APL.  Stopped at 1.2M the APL read 2.0-2.2% above its 2.1-2.4M value.
#: The earlier advice here -- never lower these ceilings on the strength of a fast sweep -- still
#: holds for a budget meant to be converged in; it does not apply to a cap that the build is
#: designed to hit.


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
        """The quilt budget must clear its slowest replica; the patch cap must clear patchB's.

        Both measured in each stage's OWN steps -- see the note at the top of this module for why
        the log's "CONVERGED at step N" is not that number.  The `patch` protocol runs twice,
        patchA then patchB, so its cap must still leave patchB, which does converge, room to.
        """
        cfg = self._example('17', 'hiv-mpertm3-membrane2.yaml')
        ceilings = self._ceilings(cfg, 'membrane_equilibrate')
        pre_embed = dict((p.split('relaxation_protocols/')[1].split('[')[0], v)
                         for p, v in self._pre_embed(ceilings))
        self.assertEqual(set(pre_embed), {'patch', 'quilt'},
                         f'expected the patch and quilt budgets, got {pre_embed}')
        self.assertGreater(pre_embed['quilt'], 1.5 * SLOWEST_OBSERVED_QUILT_CONVERGENCE)
        self.assertGreater(pre_embed['patch'], 2 * SLOWEST_OBSERVED_PATCHB_CONVERGENCE)

    def test_a_patch_cap_below_patchA_tail_needs_leaflet_ratio(self):
        """patchA does not converge within 2.5M steps, so a patch budget below that is a cap the
        build hits -- and an uncertified patchA may set the leaflet COUNTS only if nothing
        else does.  With `bilayer.leaflet_ratio` unset the tiled route takes the ratio from the
        patches, and a capped patchA would hand it a ratio that depends on where it stopped:
        exactly the irreproducibility leaflet_ratio was added to remove."""
        cfg = self._example('17', 'hiv-mpertm3-membrane2.yaml')
        ceilings = self._ceilings(cfg, 'membrane_equilibrate')
        patch = [v for p, v in self._pre_embed(ceilings) if 'relaxation_protocols/patch' in p]
        self.assertEqual(len(patch), 1, ceilings)
        if patch[0] < 2500000:
            ratios = []

            def walk(node):
                if isinstance(node, dict):
                    if 'relaxation_protocols' in node:
                        ratios.append(node.get('leaflet_ratio'))
                    for v in node.values():
                        walk(v)
                elif isinstance(node, list):
                    for v in node:
                        walk(v)
            walk(cfg)
            self.assertEqual(len(ratios), 1, 'expected one bilayer block')
            self.assertIsInstance(ratios[0], float,
                                  f'patch budget {patch[0]} caps patchA below its observed tail, '
                                  'so bilayer.leaflet_ratio must be set')

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
