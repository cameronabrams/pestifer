#Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A class for handling molecules
"""
import logging
import os

from pathlib import Path
from pidibble.pdbparse import PDBParser

from .asymmetricunit import AsymmetricUnit
from .bioassemb import BioAssembList, BioAssemb
from .chainidmanager import ChainIDManager
from .segment import Segment

from ..core.objmanager import ObjManager

from ..objs.cleavagesite import CleavageSiteList
from ..objs.graft import Graft

from ..util.cifutil import CIFload
from ..util.stringthings import my_logger

logger=logging.getLogger(__name__)

class Molecule:
    """
    A class for handling molecules, including their asymmetric unit and biological assemblies.
    This class is initialized with a source dictionary that can contain various specifications
    such as PDB or mmCIF identifiers, prebuilt structures, or AlphaFold predictions
    and manages the parsing of these structures into an asymmetric unit and biological assemblies.
    
    Parameters
    ----------
    source : dict
        A dictionary containing the source specifications for the molecule. It can include:

        - ``id``: PDB or mmCIF identifier (e.g., {``id``: ``1ABC``, ``file_format``: ``PDB``})
        - ``prebuilt``: A dictionary with ``psf`` and ``pdb`` keys for prebuilt structures (e.g., {``prebuilt``: {``psf``: ``structure.psf``, ``pdb``: ``structure.pdb``}})
        - ``alphafold``: A dictionary with AlphaFold specifications (e.g., {``alphafold``: {``model``: ``AF-1234``}})
    
    objmanager : ObjManager, optional
        An instance of ObjManager to manage objects within the molecule. If not provided,
        a new ObjManager will be created.
    chainIDmanager : ChainIDManager, optional
        An instance of ChainIDManager to manage chain IDs within the molecule. If not provided,
        a new ChainIDManager will be created.
    reset_counter : bool, optional
        If True, resets the molecule counter to 0. This is useful for testing or reinitialization purposes.
        Default is False.

    Attributes
    ----------
    molid : int
        Unique identifier for the molecule instance, automatically incremented with each new instance.
    objmanager : ObjManager
        An instance of ObjManager that manages objects within the molecule.
    chainIDmanager : ChainIDManager
        An instance of ChainIDManager that manages chain IDs within the molecule.
    sourcespecs : dict
        A dictionary containing the source specifications for the molecule, such as PDB or mmCIF identifiers,
        prebuilt structures, or AlphaFold predictions.
    asymmetric_unit : AsymmetricUnit
        An instance of AsymmetricUnit representing the asymmetric unit of the molecule.
    biological_assemblies : BioAssembList
        An instance of BioAssembList containing the biological assemblies derived from the asymmetric unit.
    parsed_struct : dict
        A dictionary containing the parsed structure of the molecule, which may include atoms, residues, and
        segments parsed from the source specifications.
    rcsb_file_format : str
        The file format of the source specifications, either ``PDB`` or ``mmCIF``.
    """

    """
    Optional attributes for the Molecule class.
    
    Attributes
    ----------
    active_biological_assembly : BioAssemb
        An instance of BioAssemb representing the currently active biological assembly.
        This assembly is derived from the asymmetric unit and can be activated based on user requests.
    """

    _molcounter: int = 0

    def __init__(self, source: dict = {}, objmanager: ObjManager = None, chainIDmanager: ChainIDManager=None, **kwargs):
        psf = None
        reset = kwargs.get('reset_counter', False)
        file_format = source.get('file_format', 'PDB')
        if reset:
            Molecule._molcounter = 0
        else:
            Molecule._molcounter += 1
        if not source:
            logger.debug('Molecule initialized without source.')
            p_struct = None
        else:
            # logger.debug('Molecule initialization')
            if 'source_id' in source and 'source_db' in source:
                ext = '.pdb' if file_format in ['PDB', 'pdb'] else '.cif'
                source_name = source['source_id'] + ext
                source_path = Path(source_name)
                if source_path.exists():
                    if file_format in ['PDB', 'pdb']:
                        p_struct = PDBParser(filepath=source_path).parse().parsed
                    elif file_format in ['mmCIF', 'cif']:
                        logger.debug(f'CIF source {source['source_id']}')
                        p_struct = CIFload(source_path)
                else:
                    if file_format in ['PDB', 'pdb']:
                        p_struct = PDBParser(source_id=source['source_id'], source_db=source['source_db']).parse().parsed
                    elif file_format in ['mmCIF', 'cif']:
                        logger.debug(f'CIF source {source['source_id']}')
                        p_struct = CIFload(source['source_id'], source_db=source['source_db'])
            elif 'prebuilt' in source:
                logger.debug(f'Prebuilt record:')
                my_logger(source["prebuilt"], logger.debug)
                psf = source['prebuilt']['psf']
                pdb = source['prebuilt']['pdb']
                xsc = source['prebuilt'].get('xsc', '')
                logger.debug(f'Using prebuilt psf {psf} and pdb {pdb}')
                p_struct = PDBParser(filepath=pdb).parse().parsed
            else:
                logger.debug(f'Neither "source_id/source_db" or "prebuilt" specified; initializing an empty molecule')
                p_struct = None
        if objmanager is None:
            logger.debug(f'Making an empty ObjManager')
            objmanager = ObjManager()
        if chainIDmanager is None:
            logger.debug(f'Molecule instantiating its own ChainIDManager')
            chainIDmanager = ChainIDManager()
        self.sourcespecs = source
        self.objmanager = objmanager
        self.chainIDmanager = chainIDmanager
        self.rcsb_file_format = file_format
        self.molid = kwargs.get('molid', Molecule._molcounter)
        self.parsed_struct = p_struct
        self.asymmetric_unit = AsymmetricUnit(parsed=p_struct, 
                                              sourcespecs=source, 
                                              objmanager=objmanager, 
                                              chainIDmanager=chainIDmanager, 
                                              psf=psf)
        self.asymmetric_unit.set_parent_molecule(self)
        self.biological_assemblies = BioAssembList(p_struct)
        self.biological_assemblies.set_parent_molecule(self)

    def __repr__(self):
        """
        Return a string representation of the Molecule instance.
        This method provides a concise description of the molecule, including its ID and the number of atoms.
        
        Returns
        -------
        str
            A string representation of the Molecule instance.
        """
        return f'Molecule(molid={self.molid}, num_atoms={self.num_atoms()})'

    def set_coords(self, altcoordsfile):
        """
        Set the coordinates of the asymmetric unit from an alternate coordinates file.
        This method reads the alternate coordinates from a PDB file and updates the asymmetric unit's coordinates.
        
        Parameters
        ----------
        altcoordsfile : str
            The path to the alternate coordinates file in PDB format.
        
        Raises
        ------
        AssertionError
            If the provided file is not in PDB format or if the file does not exist.
        """
        nm, ext = os.path.splitext(altcoordsfile)
        assert ext == '.pdb', f'Alt-coords file must be PDB format'
        altstruct = PDBParser(filepath=altcoordsfile).parse().parsed
        self.asymmetric_unit.set_coords(altstruct)

    def activate_biological_assembly(self, index):
        """
        Activate a biological assembly by its index.
        This method sets the active biological assembly based on the provided index.
        
        Parameters
        ----------
        index : int
            The index of the biological assembly to activate. If the index is 0 or if no biological assemblies are specified,
            the asymmetric unit will be used as the biological assembly.
        
        Returns
        -------
        self : Molecule
            Returns the instance of the Molecule with the active biological assembly set.
            
        Raises
        ------
        AssertionError
            If the specified biological assembly index is invalid.
        """
        if index == 0 or len(self.biological_assemblies) == 0: # we will use the unadulterated A.U. as the B.A.
            self.active_biological_assembly = BioAssemb(self.asymmetric_unit)
            if index != 0:
                logger.warning(f'No biological assemblies specified in input structure.  Using A.U.; ignoring your request of B.A. "{index}"')
        else:
            self.active_biological_assembly = self.biological_assemblies.get(index=index)
            assert self.active_biological_assembly != None, f'No biological assembly "{index:d}" found.'
            logger.info(f'Activating biological assembly {self.active_biological_assembly.name} (idx {index})')
        self.active_biological_assembly.activate(self.asymmetric_unit, self.chainIDmanager)
        return self

    def get_chainmaps(self):
        """
        Get a mapping of chain IDs in the active biological assembly.
        This method returns a dictionary where keys are original chain IDs and values are lists of dictionaries
        containing the transform index and the mapped chain ID.
        
        Returns
        -------
        maps : dict
            A dictionary mapping original chain IDs to lists of dictionaries with transform index and mapped chain ID.
            Each entry in the list corresponds to a transformation applied to the asymmetric unit.
        """
        ba=self.active_biological_assembly
        maps={}
        for transform in ba.transforms:
            for oc,mc in transform.chainIDmap.items():
                if not oc in maps:
                    maps[oc]=[]
                maps[oc].append({'transform':transform.index,'mappedchain':mc})
        return maps

    def loop_counts(self,min_loop_length=1):
        """
        Check if the asymmetric unit contains loops (missing residues) of a specified minimum length.
        This method iterates through the segments of the asymmetric unit and counts the number of loops
        that are in the ``MISSING`` state and have a length greater than or equal to the specified minimum length.

        Parameters
        ----------
        min_loop_length : int
            The minimum length of loops to consider.

        Returns
        -------
        nloops : dict
            A dictionary containing the counts of loops in the asymmetric unit, categorized by segment type.
            The keys are segment types (e.g., ``protein``, ``nucleicacid``) and the values are the counts of loops.
            For example, ``{'protein': 3, 'nucleicacid': 2}`` indicates that there are 3 loops in protein segments
            and 2 loops in nucleic acid segments.
        """
        self.min_loop_length=min_loop_length
        self.has_protein_loops=False
        nloops={}
        nloops['protein']=0
        nloops['nucleicacid']=0
        au=self.asymmetric_unit
        for S in au.segments:
            if S.segtype=='protein':
                for b in S.subsegments:
                    if b.state=='MISSING' and b.num_items()>=min_loop_length:
                        nloops['protein']+=1
            elif S.segtype=='nucleicacid':
                for b in S.subsegments:
                    if b.state=='MISSING' and b.num_items()>=min_loop_length:
                        nloops['nucleicacid']+=1
        if nloops['protein']>0:
            self.has_protein_loops=True
        return nloops

    def nglycans(self):
        """
        Count the number of glycan segments in the asymmetric unit.
        This method iterates through the segments of the asymmetric unit and counts the number of segments
        that are of type ``glycan``.

        Returns
        -------
        nglycans : int
            The number of glycan segments found in the asymmetric unit.
        """
        nglycans=0
        au=self.asymmetric_unit
        for S in au.segments:
            if S.segtype=='glycan':
                nglycans+=1
        return nglycans

    def num_images(self):
        """
        Count the number of images in the active biological assembly.
        This method returns the number of transforms applied to the asymmetric unit in the active biological assembly.
        
        Returns
        -------
        num_images : int
            The number of images (transforms) in the active biological assembly."""
        return len(self.active_biological_assembly.transforms)

    def num_atoms(self):
        """
        Count the number of atoms in the asymmetric unit.
        This method returns the total number of atoms present in the asymmetric unit.
        
        Returns
        -------
        num_atoms : int
            The total number of atoms present in the asymmetric unit.
        """
        return len(self.asymmetric_unit.atoms)
    
    def num_residues(self):
        """
        Count the number of residues in the asymmetric unit.
        This method returns the total number of residues present in the asymmetric unit.

        Returns
        -------
        num_residues : int
            The total number of residues present in the asymmetric unit.
        """
        return len(self.asymmetric_unit.residues)

    def num_segments(self):
        """
        Count the number of segments in the asymmetric unit.
        This method returns the total number of segments present in the asymmetric unit.
        
        Returns
        -------
        num_segments : int
            The total number of segments present in the asymmetric unit.
        """
        return len(self.asymmetric_unit.segments)

    # def write_protein_loop_lines(self,writer,cycles=100,**options):
    #     """
    #     Write Tcl commands to declare loops in the asymmetric unit that are in the ``MISSING``
    #     state and have a length greater than or equal to the specified minimum length.
    #     This method iterates through the segments of the asymmetric unit and generates Tcl commands
    #     to declare the loops.
        
    #     Parameters
    #     ----------
    #     writer : scriptwriter
    #         An instance of a script writer that will be used to write the Tcl commands.
    #     cycles : int, optional
    #         The number of cycles for the loop declaration. Default is 100.
    #     options : dict, optional
    #         Additional options for loop declaration. It can include:

    #         - ``min_length``: The minimum length of loops to consider. Default is 4.
    #         - ``include_c_termini``: If True, includes C-terminal loops in the declaration. If False, skips C-terminal loops. Default is False.
    #     """
    #     ba=self.active_biological_assembly
    #     au=self.asymmetric_unit
    #     min_length=self.min_loop_length
    #     include_c_termini=options.get('include_c_termini',False)
    #     for S in au.segments:
    #         # chainID=S.chainID
    #         if S.segtype in 'protein':
    #             asymm_segname=S.segname
    #             n_subsegs=len(S.subsegments)
    #             for b in S.subsegments:
    #                 is_c_terminus=(S.subsegments.index(b)==(n_subsegs-1))
    #                 is_processible=b.state=='MISSING' and b.num_items()>=min_length
    #                 if is_processible and (not include_c_termini) and is_c_terminus:
    #                     logger.debug(f'A.U. C-terminal loop {b.pstr()} declashing is skipped')
    #                     is_processible=False
    #                 if is_processible:
    #                     reslist=[f'{r.resid.resid}' for r in S.residues[b.bounds[0]:b.bounds[1]+1]]
    #                     tcllist='[list '+' '.join(reslist)+']'
    #                     for transform in ba.transforms:
    #                         cm=transform.chainIDmap
    #                         act_segID=cm.get(asymm_segname,asymm_segname)
    #                         writer.addline(f'declash_loop $mLL {act_segID} {tcllist} {cycles}')

    def write_gaps(self,writer,min_length=4):
        """
        Write Tcl commands to declare gaps in the asymmetric unit that are in the ``MISSING``
        state and have a length greater than or equal to the specified minimum length.
        This method iterates through the segments of the asymmetric unit and generates Tcl commands
        to declare the gaps.
        
        Parameters
        ----------
        writer : scriptwriter
            An instance of a script writer that will be used to write the Tcl commands.
        min_length : int, optional
            The minimum length of gaps to consider. Default is 4.
        """
        ba=self.active_biological_assembly
        au=self.asymmetric_unit
        writer.addline('# fields: segname loop-begin-res loop-end-res connect-to-res')
        for S in au.segments:
            # chainID=S.chainID
            if S.segtype=='protein':
                asymm_segname=S.segname
                for i,b in enumerate(S.subsegments):
                    if b.state=='MISSING':
                        if b.num_items()>=min_length and i<(len(S.subsegments)-1):
                            reslist=[f'{r.resid.resid}' for r in S.residues[b.bounds[0]:b.bounds[1]+1]]
                            bpp=S.subsegments[i+1]
                            nreslist=[f'{r.resid.resid}' for r in S.residues[bpp.bounds[0]:bpp.bounds[1]+1]]
                            assert bpp.state=='RESOLVED'
                            for transform in ba.transforms:
                                cm=transform.chainIDmap
                                act_segID=cm.get(asymm_segname,asymm_segname)
                                writer.addline(f'{act_segID} {reslist[0]} {reslist[-1]} {nreslist[0]}')

    def write_connect_patches(self,writer,min_length=4):
        """
        Write Tcl commands to create connect patches (``LINK``s) for gaps in the asymmetric unit.
        This method iterates through the segments of the asymmetric unit and generates Tcl commands
        to create connect patches for gaps that are in the ``MISSING`` state and have a length
        greater than or equal to the specified minimum length.

        Parameters
        ----------
        writer : scriptwriter
            An instance of a script writer that will be used to write the Tcl commands.
        min_length : int, optional
            The minimum length of gaps to consider. Default is 4.
        """
        ba=self.active_biological_assembly
        au=self.asymmetric_unit
        for S in au.segments:
            # chainID=S.chainID
            if S.segtype=='protein':
                asymm_segname=S.segname
                for i,b in enumerate(S.subsegments):
                    if b.state=='MISSING':
                        if b.num_items()>=min_length and i<(len(S.subsegments)-1):
                            llres=S.residues[b.bounds[1]-1]
                            lres=S.residues[b.bounds[1]]
                            nextb=S.subsegments[i+1]
                            rres=S.residues[nextb.bounds[0]]
                            rrres=S.residues[nextb.bounds[0]+1]
                            for transform in ba.transforms:
                                cm=transform.chainIDmap
                                act_segID=cm.get(asymm_segname,asymm_segname)
                                ll=f'{act_segID}:{llres.resid.resid}'
                                l= f'{act_segID}:{ lres.resid.resid}'
                                r= f'{act_segID}:{ rres.resid.resid}'
                                rr=f'{act_segID}:{rrres.resid.resid}'
                                # undo the C-terminus on the left residue
                                writer.addline(f'patch XCTR {l}')
                                # under the N-terminus on the right residue
                                rname=rres.resname
                                if rname=='GLY':
                                    writer.addline(f'patch XGLP {r}')
                                elif rname=='PRO':
                                    writer.addline(f'patch XPRP {r}')
                                else:
                                    writer.addline(f'patch XNTR {r}')
                                # re-establish the peptide bond linkage
                                writer.addline(f'patch LINK {ll} {l} {r} {rr}')
                                # writer.addline(f'patch HEAL {ll} {l} {r} {rr}')

    def cleave_chains(self, clv_list: CleavageSiteList):
        """
        Cleave segments in the asymmetric unit based on a list of cleavage specifications.
        This method iterates through the list of cleavage specifications, finds the corresponding segments in the asymmetric unit,
        and performs the cleavage operation. It also updates the chain IDs of disulfide bonds and links that are affected by the cleavage.

        Parameters
        ----------
        clv_list : list
            A list of cleavage specifications to apply to the asymmetric unit.
        """
        au = self.asymmetric_unit
        cm = self.chainIDmanager
        topomods = self.objmanager.get('topol', {})
        for clv in clv_list:
            S = au.segments.get(lambda x: x.chainID == clv.chainID)
            if S:
                DchainID = cm.next_unused_chainID()
                logger.debug(f'Cleaving segment {S.segname} at {clv} to make new segment with chainID {DchainID}')
                D = S.cleave(clv, DchainID)
                au.add_segment(D)

                ssbonds = topomods.get('ssbonds', [])
                if ssbonds:
                    for c in ssbonds.filter(lambda x: x.chainID1 == S.segname):
                        for d in D.residues.filter(lambda x: x.resid == c.resid1):
                            c.chainID1 = d.chainID
                    for c in ssbonds.filter(lambda x: x.chainID2 == S.segname):
                        for d in D.residues.filter(lambda x: x.resid == c.resid2):
                            c.chainID2 = d.chainID
                links = topomods.get('links', [])
                if links:
                    logger.debug(f'Examining {len(links)} links for ones needing chainID update from {S.segname} to {DchainID} due to cleavage')
                    for c in links.filter(lambda x: x.chainID1 == S.segname):
                        logger.debug(f'...link {str(c)}')
                        for d in D.residues.filter(lambda x: x.resid == c.resid1):
                            logger.debug(f'...hit! {c.chainID1}->{d.chainID}')
                            c.update_residue(1,chainID=d.chainID)
                    for c in links.filter(lambda x: x.chainID2 == S.segname):
                        logger.debug(f'...link {str(c)}')
                        for d in D.residues.filter(lambda x: x.resid == c.resid2):
                            logger.debug(f'...hit! {c.chainID2}->{d.chainID}')
                            c.update_residue(2,chainID=d.chainID)
                    for c in links:
                        logger.debug(str(c))
            else:
                logger.debug(f'No segment with chainID {clv.chainID} found; no cleavage performed.')

Graft.model_rebuild()
Segment.model_rebuild()
