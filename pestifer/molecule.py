#Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""A class for handling molecules
"""
import logging
import os
from .objmanager import ObjManager
from pidibble.pdbparse import PDBParser
from .util.cifutil import CIFload
from .baseobj import AncestorAwareObj
from .asymmetricunit import AsymmetricUnit
from .bioassemb import BioAssembList,BioAssemb
from .scriptwriters import Psfgen, Filewriter
from .chainidmanager import ChainIDManager
logger=logging.getLogger(__name__)

class Molecule(AncestorAwareObj):
    req_attr=AncestorAwareObj.req_attr+['molid','objmanager','chainIDmanager','sourcespecs','asymmetric_unit','biological_assemblies','parsed_struct','rcsb_file_format']
    opt_attr=AncestorAwareObj.opt_attr+['active_biological_assembly']
    _molcounter=0
    def __init__(self,source={},objmanager=None,chainIDmanager=None,**kwargs):
        psf=None
        reset=kwargs.get('reset_counter',False)
        if reset:
            Molecule._molcounter=0
        else:
            Molecule._molcounter+=1
        if not source:
            logger.debug('Molecule initialized without source.')
            p_struct=None
        else:
            logger.debug('Molecule initialization')
            file_format=source.get('file_format','PDB')
            if source.get('id',{}) or source.get('prebuilt',{}) or source.get('alphafold',{}):
                if source.get('id',{}):
                    logger.debug(f'Molecule initialization from file {source["id"]}.{source["file_format"]}')
                    if file_format=='PDB':
                        p_struct=PDBParser(PDBcode=source['id']).parse().parsed
                    elif file_format=='mmCIF':
                        logger.debug(f'CIF source {source["id"]}')
                        p_struct=CIFload(source['id'])
                elif source.get('prebuilt',{}):
                    logger.debug(f'prebuilt rec [{source["prebuilt"]}]')
                    psf=source['prebuilt']['psf']
                    pdb=source['prebuilt']['pdb']
                    xsc=source['prebuilt'].get('xsc','')
                    logger.debug(f'Using prebuilt psf {psf} and pdb {pdb}')
                    pdb_pseudocode,ext=os.path.splitext(pdb)
                    p_struct=PDBParser(PDBcode=pdb_pseudocode).parse().parsed
                elif source.get('alphafold',{}):
                    ac=source['alphafold']
                    p_struct=PDBParser(alphafold=ac).parse().parsed
            else:
                logger.debug(f'None of "id", "prebuilt", or "alphafold" specified; initializing an empty molecule')
                p_struct=None
        if objmanager==None:
            logger.debug(f'Making an empty ObjManager')
            objmanager=ObjManager()
        if chainIDmanager==None:
            logger.debug(f'Molecule instantiating its own ChainIDManager')
            chainIDmanager=ChainIDManager()
        input_dict={
            'sourcespecs': source,
            'objmanager': objmanager,
            'chainIDmanager':chainIDmanager,
            'rcsb_file_format': file_format,
            'molid': Molecule._molcounter,
            'parsed_struct': p_struct,
            'asymmetric_unit': AsymmetricUnit(parsed=p_struct,sourcespecs=source,objmanager=objmanager,chainIDmanager=chainIDmanager,psf=psf),
            'biological_assemblies': BioAssembList(p_struct)
        }
        super().__init__(input_dict)
        self.asymmetric_unit.claim_descendants(self)
        self.biological_assemblies.claim_descendants(self)

    def set_coords(self,altcoordsfile):
        nm,ext=os.path.splitext(altcoordsfile)
        assert ext=='.pdb',f'Alt-coords file must be PDB format'
        altstruct=PDBParser(PDBcode=nm).parse().parsed
        self.asymmetric_unit.set_coords(altstruct)

    def activate_biological_assembly(self,index):
        if index==0 or len(self.biological_assemblies)==0: # we will use the unadulterated A.U. as the B.A.
            self.active_biological_assembly=BioAssemb(self.asymmetric_unit)
            if index!=0:
                logger.warning(f'No biological assemblies specified in input structure.  Using A.U.; ignoring your request of B.A. "{index}"')
        else:
            self.active_biological_assembly=self.biological_assemblies.get(index=index)
            assert self.activate_biological_assembly!=None,f'No biological assembly "{index:d}" found.'
            logger.info(f'Activating biological assembly {self.active_biological_assembly.name} (idx {index})')
        self.active_biological_assembly.activate(self.asymmetric_unit,self.chainIDmanager)
        return self

    def write_TcL(self,W:Psfgen):
        au=self.asymmetric_unit
        segments=au.segments
        topomods=au.objmanager.get('topol',{})
        ssbonds=topomods.get('ssbonds',[])
        links=topomods.get('links',[])
        ba=self.active_biological_assembly
        for transform in ba.transforms:
            W.banner(f'Transform {transform.index} begins')
            W.banner('The following mappings of A.U. asym ids is used:')
            for k,v in transform.chainIDmap.items():
                W.comment(f'A.U. chain {k}: Image chain {v}')
            W.banner('Segments follow')
            segments.write_TcL(W,transform)
            W.banner('DISU patches follow')
            if ssbonds:
                ssbonds.write_TcL(W,transform)
            W.banner('LINK patches follow')
            if links:
                A=Filewriter().newfile('linkreport.txt')
                links.report(A)
                A.writefile()
                links.write_TcL(W,transform)
            W.banner(f'Transform {transform.index} ends')

    def get_chainmaps(self):
        ba=self.active_biological_assembly
        maps={}
        for transform in ba.transforms:
            for oc,mc in transform.chainIDmap.items():
                if not oc in maps:
                    maps[oc]=[]
                maps[oc].append({'transform':transform.index,'mappedchain':mc})
        return maps

    def has_loops(self,min_loop_length=1):
        nloops=0
        au=self.asymmetric_unit
        for S in au.segments:
            # chainID=S.chainID
            if S.segtype=='protein':
                for b in S.subsegments:
                    if b.state=='MISSING' and b.num_items()>=min_loop_length:
                        nloops+=1
        return nloops

    def nglycans(self):
        nglycans=0
        au=self.asymmetric_unit
        for S in au.segments:
            if S.segtype=='glycan':
                nglycans+=1
        return nglycans

    def num_images(self):
        return len(self.active_biological_assembly.transforms)

    def num_atoms(self):
        return len(self.asymmetric_unit.atoms)
    def num_residues(self):
        return len(self.asymmetric_unit.residues)
    def num_segments(self):
        return len(self.asymmetric_unit.segments)

    def write_loop_lines(self,writer,cycles=100,**options): #min_length=4):
        ba=self.active_biological_assembly
        au=self.asymmetric_unit
        min_length=options.get('min_length',4)
        include_c_termini=options.get('include_c_termini',False)
        for S in au.segments:
            # chainID=S.chainID
            if S.segtype=='protein':
                asymm_segname=S.segname
                n_subsegs=len(S.subsegments)
                for b in S.subsegments:
                    is_c_terminus=(S.subsegments.index(b)==(n_subsegs-1))
                    is_processible=b.state=='MISSING' and b.num_items()>=min_length
                    if is_processible and (not include_c_termini) and is_c_terminus:
                        logger.debug(f'A.U. C-terminal loop {b.pstr()} declashing is skipped')
                        is_processible=False
                    if is_processible:
                        reslist=[f'{r.resseqnum}{r.insertion}' for r in S.residues[b.bounds[0]:b.bounds[1]+1]]
                        tcllist='[list '+' '.join(reslist)+']'
                        for transform in ba.transforms:
                            cm=transform.chainIDmap
                            act_segID=cm.get(asymm_segname,asymm_segname)
                            writer.addline(f'declash_loop $mLL {act_segID} {tcllist} {cycles}')

    def write_gaps(self,writer,min_length=4):
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
                            reslist=[f'{r.resseqnum}{r.insertion}' for r in S.residues[b.bounds[0]:b.bounds[1]+1]]
                            bpp=S.subsegments[i+1]
                            nreslist=[f'{r.resseqnum}{r.insertion}' for r in S.residues[bpp.bounds[0]:bpp.bounds[1]+1]]
                            assert bpp.state=='RESOLVED'
                            for transform in ba.transforms:
                                cm=transform.chainIDmap
                                act_segID=cm.get(asymm_segname,asymm_segname)
                                writer.addline(f'{act_segID} {reslist[0]} {reslist[-1]} {nreslist[0]}')

    def write_connect_patches(self,writer,min_length=4):
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
                                ll=f'{act_segID}:{llres.resseqnum}{llres.insertion}'
                                l= f'{act_segID}:{ lres.resseqnum}{ lres.insertion}'
                                r= f'{act_segID}:{ rres.resseqnum}{ rres.insertion}'
                                rr=f'{act_segID}:{rrres.resseqnum}{rrres.insertion}'
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

    def cleave_chains(self,clv_list):
        au=self.asymmetric_unit
        cm=self.chainIDmanager
        topomods=self.objmanager.get('topol',{})
        for clv in clv_list:
            S=au.segments.get(chainID=clv.chainID)
            if S:
                DchainID=cm.next_unused_chainID()
                D=S.cleave(clv,DchainID)
                au.add_segment(D)

                ssbonds=topomods.get('ssbonds',[])
                if ssbonds:
                    for c in ssbonds.filter(chainID1=S.segname):
                        for d in D.residues.filter(resseqnum=c.resseqnum1,insertion=c.insertion1):
                            c.chainID1=d.chainID
                    for c in ssbonds.filter(chainID2=S.segname):
                        for d in D.residues.filter(resseqnum=c.resseqnum2,insertion=c.insertion2):
                            c.chainID2=d.chainID
                links=topomods.get('links',[])
                if links:
                    logger.debug(f'Examining {len(links)} links for ones needing chainID update from {S.segname} to {DchainID} due to cleavage')
                    for c in links.filter(chainID1=S.segname):
                        logger.debug(f'...link {str(c)}')
                        for d in D.residues.filter(resseqnum=c.resseqnum1,insertion=c.insertion1):
                            logger.debug(f'...hit! {c.chainID1}->{d.chainID}')
                            c.update_residue(1,chainID=d.chainID)
                    for c in links.filter(chainID2=S.segname):
                        logger.debug(f'...link {str(c)}')
                        for d in D.residues.filter(resseqnum=c.resseqnum2,insertion=c.insertion2):
                            logger.debug(f'...hit! {c.chainID2}->{d.chainID}')
                            c.update_residue(2,chainID=d.chainID)
                    for c in links:
                        logger.debug(str(c))
            else:
                logger.debug(f'No segment with chainID {clv.chainID} found; no cleavage performed.')