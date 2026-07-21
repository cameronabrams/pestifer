import unittest
import logging
import pytest
from pathlib import Path
logger=logging.getLogger(__name__)

pytestmark = pytest.mark.needs_tools
from pidibble.pdbparse import PDBParser
from pidibble.pdbrecord import PDBRecordDict
from pestifer.molecule.residue import ResidueList, ResiduePlaceholder, ResiduePlaceholderList
from pestifer.molecule.atom import Atom, AtomList
from pestifer.objs.seqadv import Seqadv, SeqadvList
from pestifer.core.config import Config
from pestifer.molecule.bioassemb import BioAssembList
from pestifer.util.util import reduce_intlist
from pestifer.objs.resid import ResID
import os
import glob


def _parse_cif(source):
    """Parse ``{source}.cif`` (a cwd symlink set up below) with pidibble as mmCIF."""
    return PDBParser(filepath=f'{source}.cif', input_format='mmCIF').parse().parsed


class TestCIF(unittest.TestCase):

    def setUp(self):
        input_dir = Path(__file__).parents[2] / "inputs"
        zmj = input_dir / '4zmj.cif'
        fae = input_dir / '8fae.cif'
        # copy to cwd
        dest_zmj = Path('4zmj.cif')
        dest_fae = Path('8fae.cif')
        if dest_zmj.exists():
            dest_zmj.unlink()
        if dest_fae.exists():
            dest_fae.unlink()
        os.symlink(zmj.resolve(), dest_zmj)
        os.symlink(fae.resolve(), dest_fae)

    def tearDown(self):
        dest_zmj = Path('4zmj.cif')
        dest_fae = Path('8fae.cif')
        if dest_zmj.exists():
            dest_zmj.unlink()
        if dest_fae.exists():
            dest_fae.unlink()
        logs = Path('.').glob('*.log')
        for log in logs:
            log.unlink()

    def test_load(self):
        p_struct = _parse_cif('8fae')
        self.assertIsInstance(p_struct, PDBRecordDict)

    def test_missing_residue_record(self):
        # pidibble maps mmCIF pdbx_unobs_or_zero_occ_residues onto REMARK.465, with each
        # row carrying both a label and an author identity.
        p_struct = _parse_cif('8fae')
        missing = p_struct['REMARK.465'].tables['MISSING']
        self.assertEqual(len(missing), 22)
        row = missing[0]
        for attr in ['resName', 'chainID', 'seqNum', 'iCode',
                     'auth_chainID', 'auth_resName', 'auth_seqNum']:
            self.assertTrue(hasattr(row, attr), f'missing attr {attr}')

    def test_atoms(self):
        atoms = AtomList.from_cif(_parse_cif('8fae'))
        self.assertEqual(len(atoms), 17693)

    def test_residues(self):
        atoms = AtomList.from_cif(_parse_cif('8fae'))
        self.assertEqual(len(atoms), 17693)
        residues = ResidueList.from_residuegrouped_atomlist(atoms)
        self.assertEqual(len(residues), 2082)
        uCIDs = residues.uniqattrs(['chainID'])['chainID']
        self.assertEqual(len(uCIDs), 82)
        nres = 0
        for c in uCIDs:
            chain = residues.filter(lambda x: x.chainID == c)
            nres += len(chain)
        self.assertEqual(len(residues), nres)

    def test_vmd(self):
        source = '8fae'
        config = Config().configure_new()
        vmd = config.get_scripter('vmd')
        vmd.newscript('testcif')
        vmd.addline(f'mol new {source}.cif')
        atoms = AtomList.from_cif(_parse_cif(source))
        raw_residues = ResidueList.from_residuegrouped_atomlist(atoms)
        raw_residues.apply_segtypes()
        residues = raw_residues.get(lambda x: x.segtype == 'protein')  # 8fae has some glycan residues with resid 0
        uCIDs = residues.uniqattrs(['chainID'])['chainID']
        nres = 0
        for c in uCIDs:
            chain = residues.filter(lambda x: x.chainID == c)
            resids = []
            for x in chain:
                resids.extend([str(y.resid) for y in x.atoms])
            residlist = ' '.join(resids)
            serials = chain.atom_serials(as_type=int)
            vmd_red_list = reduce_intlist(serials)
            vmd.addline(f'set a [atomselect top "serial {vmd_red_list}"]')
            vmd.addline(f'set c [lsort -unique [$a get chain]]')
            vmd.addline(f'$a set chain {c}')
            vmd.addline(f'$a set resid [ list {residlist} ]')
            vmd.addline(f'set cn [lsort -unique [$a get chain]]')
            vmd.addline(f'set resids [$a get resid]')
            vmd.addline(f'puts "CHAIN $cn {c}"')
            vmd.addline(f'set b [atomselect top "chain {c}"]')
            vmd.addline(f'puts "COUNTS [$b num] {len(serials)}"')
            vmd.addline(f'puts "RESIDS && $resids && {residlist}"')
            vmd.addline(f'foreach avmd $resids apyt [list {residlist}] ' + r'{')
            vmd.addline(f'puts "   $avmd $apyt"')
            vmd.addline(r'}')
            nres += len(chain)
        vmd.writescript()
        vmd.runscript()
        with open('testcif.log', 'r') as f:
            output = f.read().split('\n')
        old_logs = glob.glob('%*%')
        for ol in old_logs:
            os.remove(ol)
        self.assertTrue(len(output) > 0)
        for l in output:
            if l.startswith('CHAIN'):
                fields = l.split()
                cvmd = fields[1]
                ccif = fields[2]
                self.assertEqual(ccif, cvmd)
            if l.startswith('COUNTS'):
                fields = l.split()
                cvmd = fields[1]
                ccif = fields[2]
                self.assertEqual(ccif, cvmd)
            if l.startswith('RESIDS'):
                fields = l.split('&&')
                rvmd = fields[1].split()
                rcif = fields[2].split()
                for i, j in zip(rvmd, rcif):
                    self.assertEqual(i, str(ResID.split_ri(j)[0]))  # CIF files do not have insertion codes on their native residue sequence numbers
        Path('testcif.tcl').unlink()

    def test_biomolassemb_cif(self):
        # pidibble normalizes mmCIF biological assemblies onto REMARK.350 records
        BAList = BioAssembList(_parse_cif('8fae'))
        self.assertEqual(len(BAList), 1)

    def test_seqadv_cif(self):
        p = _parse_cif('4zmj')
        Atoms = AtomList.from_cif(p)
        Seqadvs = SeqadvList.from_cif(p)
        EmptyResidues = ResiduePlaceholderList.from_cif(p)
        fromAtoms = ResidueList.from_residuegrouped_atomlist(Atoms)
        fromEmptyResidues = ResidueList.from_ResiduePlaceholderlist(EmptyResidues)
        Residues = fromAtoms + fromEmptyResidues
        Seqadvs.assign_residues(Residues)
        for s in Seqadvs:
            myres = Residues.get(lambda x: x.chainID == s.chainID and x.resid == s.resid)
            logger.debug(f'Seqadv {s.chainID}:{s.resid} typekey {s.typekey} residue {myres.resid if myres else None}')
            if not myres:
                self.assertTrue('engineered' not in s.typekey and 'conflict' not in s.typekey)
            else:
                self.assertTrue('engineered' in s.typekey or 'conflict' not in s.typekey)
