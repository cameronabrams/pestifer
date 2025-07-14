"""
CHARMM Force Field Residue Database
This module defines the :class:`CHARMMFFResiDatabase` class, which manages a database of CHARMM force field residues.
"""

import logging
from .charmmffcontent import CHARMMFFContent, CHARMMFFStreamID
from .charmmtop import CharmmMasses

logger=logging.getLogger(__name__)

class CHARMMFFResiDatabase:
    """ 
    A class for handling a database of CHARMM force field residues.
    This class loads residues and patches from CHARMM force field topology/stream files, including toplevel files and specific streams.
    It maintains a collection of residues, patches, and masses, and provides methods to access and manipulate this data.
    
    Attributes
    ----------
    charmmff_content : CHARMMFFContent
        An instance of CHARMMFFContent containing the CHARMM force field data.
    streamIDs : list of str
        A list of stream IDs from which residues can be loaded.
    residues : dict
        A dictionary mapping residue names to their corresponding CharmmTopResi objects.
    patches : dict
        A dictionary mapping patch names to their corresponding CharmmTopResi objects.
    masses : CharmmMasses
        An instance of CharmmMasses containing the atom mass records extracted from the CHARMM force field files.
    overrides : dict
        A dictionary containing overrides for specific substreams, such as lipid classification.
    """
    def __init__(self,charmmff_content:CHARMMFFContent,streamIDs=[]):
        logger.debug(f'Initializing CHARMMFFResiDatabase with extra streamIDs {streamIDs}')
        self.charmmff_content=charmmff_content
        self.residues={}
        self.patches={}
        self.masses=CharmmMasses({})
        self.streamIDs=streamIDs
        self.load_charmmresi_from_toplevels()
        logger.debug(f'Loaded {len(self.residues)} residues from toplevels, streams: {self.streamIDs}')
        for streamID in streamIDs:
            logger.debug(f'Loading residues from stream {streamID}')
            self.load_charmmresi_from_stream(streamID)
            if not streamID in self.streamIDs:
                logger.debug(f'Adding stream {streamID} to streams')
                # if the stream is not already in the streams list, add it
                self.streamIDs.append(streamID)
        self.tally_masses()
        self.overrides={
            'substreams':{
                'C6DHPC':'lipid',
                'C7DHPC':'lipid',
                'C8DHPC':'lipid',
                'C9DHPC':'lipid',
                'C10DHPC':'lipid',
                'C11DHPC':'lipid',
                'C12DHPC':'lipid',
                'C13DHPC':'lipid',
                'C14DHPC':'lipid',
                'C15DHPC':'lipid',
                'C16DHPC':'lipid',
                'C17DHPC':'lipid',
                'C18DHPC':'lipid'
            }
        }
        logger.debug(f'Initialized CHARMMFFResiDatabase with extra streamIDs {streamIDs}')

    def load_charmmresi_from_toplevels(self):
        """ 
        Load residues and patches from the toplevel CHARMM force field files.
        This function iterates through the toplevel topology files and extracts residues and patches.
        It updates the `self.residues` and `self.patches` dictionaries with the extracted data.
        """
        for topfile in list(self.charmmff_content.toplevel_top.values())+list(self.charmmff_content.toplevel_toppar.values()):
            new_resis,new_patches=self.load_charmmresi_from_topfile(topfile)
            logger.debug(f'Loaded {len(new_resis)} residues from toplevel {topfile}')
            self.residues.update({x.resname:x for x in new_resis})
            logger.debug(f'Loaded {len(new_patches)} patches from toplevel {topfile}')
            self.patches.update({x.resname:x for x in new_patches})

    def load_charmmresi_from_stream(self,streamID):
        """ 
        Load residues and patches from a specific stream in the CHARMM force field content.
        This function checks if the streamID exists in the CHARMM force field content and loads the residues and patches from the corresponding topology files.

        Parameters
        ----------
        streamID : str
            The ID of the stream from which to load residues and patches.

        Raises
        -------
        Warning
            If the specified streamID is not found in the CHARMM force field content, a warning is logged.
        """
        if streamID not in self.charmmff_content.streams:
            logger.warning(f'load_charmmresi_from_stream: Stream {streamID} not found in CHARMM force field content')
            return
        logger.debug(f'Loading resis from stream {streamID}')
        for topfile in self.charmmff_content.streamfiles[streamID].values():
            this_residues,this_patches=self.load_charmmresi_from_topfile(topfile)
            logger.debug(f'Loaded {len(this_residues)} residues and {len(this_patches)} patches from {topfile} in stream {streamID}')
            self.residues.update({x.resname:x for x in this_residues})
            self.patches.update({x.resname:x for x in this_patches})

    def get_resnames_of_streamID(self,streamID,substreamID=None):
        """ 
        Get the residue names of a specific stream ID.
        This function retrieves the residue names associated with a given stream ID and optional substream ID.

        Parameters
        ----------
        streamID : str
            The ID of the stream from which to retrieve residue names.
        substreamID : str, optional
            The ID of the substream from which to retrieve residue names. If None, all residues in the stream are returned.

        Returns
        -------
        list of str
            A sorted list of residue names associated with the specified stream ID and optional substream ID.

        Raises
        -------
        Warning
            If the specified streamID is not found in the CHARMM force field residue database, a warning is logged and an empty list is returned.
        """
        if streamID not in self.streamIDs:
            logger.warning(f'get_resnames_of_streamID: Stream {streamID} not found in CHARMM force field residue database')
            return []
        resnames=[x.resname for x in self.residues.values() if (x.metadata['streamID']==streamID and ((substreamID is None) or (x.metadata.get('substreamID','')==substreamID)))]
        logger.debug(f'Found {len(resnames)} residues in stream {streamID}')
        resnames.sort()
        return resnames

    def load_charmmresi_from_topfile(self,topfile):
        """ 
        Load residues and patches from a specific topology file.
        This function reads the contents of a topology file, extracts residues and patches, and updates the `self.residues` and `self.patches` dictionaries.


        Parameters
        ----------
        topfile : str
            The name of the topology file to load residues and patches from. This can be a full path or just the basename.

        Returns
        -------
        tuple of list of CharmmTopResi
            A tuple containing two lists: the first list contains the residues extracted from the topology file, and the second list contains the patches.
        """
        logger.debug(f'CHARMMFFResiDatabase::load_charmmresifrom_topfile Loading resis from {topfile}')
        cstr=CHARMMFFStreamID(topfile)
        logger.debug(f'topfile {topfile} CHARMMFFStream: \'{cstr.streamID}\' \'{cstr.substreamID}\'')
        # self.masses.update(self.charmmff_content.masses_from_topfile(topfile))
        all_resis,all_masses=self.charmmff_content.resis_and_masses_from_topfile(topfile,metadata=dict(
                streamID=cstr.streamID,
                substreamID=cstr.substreamID,
                charmmtopfile=topfile
            ))
        self.masses.update(all_masses)
        this_residues=[x for x in all_resis if x.key=='RESI']
        this_patches=[x for x in all_resis if x.key=='PRES']
        for resi in this_residues:
            if resi.resname in self.residues:
                logger.debug(f'Residue {resi.resname} found in {topfile} will overwrite already loaded {resi.resname} from {self.residues[resi.resname].metadata["charmmtopfile"]}')
        return this_residues,this_patches

    def tally_masses(self):
        """ 
        Tally the masses of all residues in the database.
        This function iterates through all residues in the database and sets their masses based on the atom mass records.
        It uses the `set_masses` method of each residue to assign the correct masses.
        This is necessary because the atom mass records are stored throughout the topology files, and we need to ensure that all residues have their masses set correctly after all topology files are read in.
        """
        for resi in self.residues.values():
            resi.set_masses(self.masses)

    def add_stream(self,streamID):
        """ 
        Add a stream to the CHARMM force field residue database.
        This function loads residues and patches from the specified stream ID and updates the `self.streamIDs` list.

        Parameters
        ----------
        streamID : str
            The ID of the stream to add to the CHARMM force field residue database.
        """
        self.load_charmmresi_from_stream(streamID)
        if not streamID in self.streamIDs:
            logger.debug(f'Adding stream {streamID} to streams')
            # if the stream is not already in the streams list, add it
            self.streamIDs.append(streamID)
        self.tally_masses()

    def add_topology(self,topfile,streamIDoverride=None):
        """ 
        Add a topology file to the CHARMM force field residue database.
        This function loads residues and patches from the specified topology file and updates the `self.residues` and `self.patches` dictionaries.
        If a `streamIDoverride` is provided, it will set the `streamID` and `substreamID` metadata for all residues and patches loaded from the topology file.

        Parameters
        ----------
        topfile : str
            The name of the topology file to load residues and patches from. This can be a full path or just the basename.
        streamIDoverride : str, optional
            An optional stream ID to override the default stream ID extracted from the topology file. If provided, it will set the `streamID` and `substreamID` metadata for all residues and patches loaded from the topology file.
        """
        new_resis,new_patches=self.load_charmmresi_from_topfile(topfile)
        if streamIDoverride is not None:
            for resi in new_resis+new_patches:
                resi.metadata['streamID']=streamIDoverride
                resi.metadata['substreamID']=''
        self.residues.update({x.resname:x for x in new_resis})
        self.patches.update({x.resname:x for x in new_patches})
        if streamIDoverride is not None and streamIDoverride not in self.streamIDs:
            self.streamIDs.append(streamIDoverride)
        self.tally_masses()

    def get_resi(self,resname):
        """ 
        Get a residue by its name from the CHARMM force field residue database.
        This function checks if the residue name exists in the database and returns the corresponding CharmmTopResi object.

        Parameters
        ----------
        resname : str
            The name of the residue to retrieve from the database.

        Returns
        -------
        CharmmTopResi or None
            The CharmmTopResi object corresponding to the residue name, or None if the residue is not found in the database.
        """
        if resname in self.residues:
            return self.residues[resname]
        else:
            logger.warning(f'Residue {resname} not found in CHARMM force field residue database')
            return None
    
    def __contains__(self,resname):
        """ 
        Check if a residue name exists in the CHARMM force field residue database.
        This function checks if the specified residue name is present in the ``self.residues`` dictionary.

        Parameters
        ----------
        resname : str
            The name of the residue to check for in the database.
            
        Returns
        -------
        bool
            True if the residue name exists in the database, False otherwise.
        """
        return resname in self.residues

