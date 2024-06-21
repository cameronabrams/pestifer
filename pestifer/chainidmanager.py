# Author: Cameron F. Abrams <cfa22@drexel.edu>.
""" Chain ID manager -- make sure all chains have unique IDs and ID maps are tracked
"""

import logging
logger=logging.getLogger(__name__)

class ChainIDManager:
    """A class for managing chain IDs. 
    
    Methods
    -------
    __init__(format, reserved_maps)
      initializes the repository of chainIDs depending on format
       
    register_chains(chainIDs)
      Registers the list of chainIDs with the manager
    
    unregister_chain(chainID)
      Unregisters the chainID
      
    next_chain()
      Registers the next available chainID and returns it
    
    generate_next_map()
      Given the original chainIDs, returns a dictionary mapping each 
      to a new chainID
    
    thru_map()
      Generates an identity map
      
    cleavage_daughter_chainID(chainID)
      returns a single-entry dictionary mapping the chainID to the
      next available chainID
    """
    def __init__(self,format='PDB',remap={},reserved_maps={}):
        logger.debug(f'New chainIDmanager, format {format}')
        self.format=format
        self.remap=remap
        self.reserved_maps=reserved_maps
        self.ReservedUnused=[]
        U=[chr(i) for i in range(ord('A'),ord('A')+26)]
        if format=='PDB':
            # use all possible 1-byte upper- and lowercase, and single digits
            L=[chr(i) for i in range(ord('a'),ord('a')+26)]
            D=[str(i) for i in range(10)]
            self.Unused=U+L+D
        elif format=='mmCIF':
            # use all possible 1- and 2-byte uppercase
            self.Unused=[]
            for a in ['']+U:
                for b in U:
                    self.Unused.append(b+a)
        for k,v in self.reserved_maps.items():
            assert k in self.Unused,f'Error: reserved map key chainID {k} is not available in the Unused set'
            for mc in v: # only put image chainIDs in ReservedUnused
                assert mc in self.Unused,f'Error: reserved map chainID {mc} is not available in the Unused set'
                self.Unused.remove(mc)
                self.ReservedUnused.append(mc)
                # at this point, the keys in reserved_maps are also elements in Unused
        self.Used=set()

    _init_registration=False
    def register_chains(self,chainIDs):
        """ Registers all chains detected in the list chainIDs with the manager 
        
        Parameters:
        -----------
        chainIDs : list
           List of chainIDs (each is a 1 or 2-byte str) to register

        Returns:
        --------
        reserved_conflicts : list
           List of chainIDs that were attempted to be registered but were already
           named in the reserved maps.  This means that a currently existing chainID
           has been requested to be assigned to a transformed subunit in a 
           biological assembly that is built from BIOMT operations. So the chainIDs
           in this list must be changed to available unused chainIDs by the caller.
        """
        assert not self._init_registration,f'Error: cannot register chains twice'
        self._init_registration=True
        reserved_conflicts=[]
        for c in chainIDs:
            if c in self.ReservedUnused:
                reserved_conflicts.append(c)
            elif c in self.Unused:
                self.Unused.remove(c)
                self.Used.add(c)
            else:
                raise ValueError(f'chainID {c} is not available for registry')
        return reserved_conflicts
    
    def unregister_chain(self,chainID):
        """ Unregister the single chainID with the manager 
        
        Paramter:
        ---------
        chainID : str
           1- or 2-byte chainID to unregister
        """
        self.Used.remove(chainID)
        is_reserved=any([chainID in v for v in self.reserved_maps.values()])
        if is_reserved:
            self.ReservedUnused.append(chainID)
        else:
            self.Unused.append(chainID)

    def next_unused_chainID(self):
        p=self.Unused.pop(0)
        self.Used.add(p)
        return p
    
    def next_reserved_chainID(self,key):
        assert key in self.reserved_maps,f'Key chainID {key} not a key in the prescribed reserve maps'
        mapsto=self.reserved_maps[key]
        assert any([x in self.ReservedUnused for x in mapsto]),f'mapped key {key} has no daughters available'
        for c in mapsto:
            if c in self.ReservedUnused:
                self.ReservedUnused.remove(c)
                self.Used.add(c)
                return c
        return None
        
    def is_already_reserved(self,tst):
        for k,v in self.reserved_maps.items():
            if tst in v:
                return True
        return False

    def generate_next_map(self,chainIDs,active_chains=[]):
        """ Generate a map identifying new chainIDs for each existing chainID
            in the list chainIDs.  If active_chains is not empty, then 
            it must list a subset of the chainIDs considered actually active
        """
        assert len(chainIDs)<=len(self.Unused),f'Not enough available chainIDs'
        myMap={}
        activeChainIDs=chainIDs.copy()
        # logger.debug(f'generating next map from {activeChainIDs} with actives {active_chains}')
        if active_chains:
            inactive_chains=[x for x in activeChainIDs if not x in active_chains]
            for i in inactive_chains:
                activeChainIDs.remove(i)
        assert not (any([c in self.Unused for c in activeChainIDs]) or any([c in self.ReservedUnused for c in activeChainIDs]))
        for c in activeChainIDs:
            if c in self.reserved_maps:
                myMap[c]=self.next_reserved_chainID(c)
            else:
                myMap[c]=self.next_unused_chainID()
        logger.debug(f'generated next chainID map: {myMap}')
        return myMap
    
    def thru_map(self,chainIDs,active_chains=[]):
        activeChainIDs=chainIDs.copy()
        if active_chains:
            inactive_chains=[x for x in activeChainIDs if not x in active_chains]
            for i in inactive_chains:
                activeChainIDs.remove(i)
        logger.debug(f'generating thru_map from {activeChainIDs} with actives {activeChainIDs}')
        return {c:c for c in activeChainIDs}
    
    def cleavage_daughter_chainID(self,chainID):
        assert 1<=len(self.Unused),f'Not enough available chainIDs'
        assert chainID in self.Unused,f'Parent chain {chainID} was never claimed'
        p=self.Unused.pop(0)
        self.Used.add(p)
        return {chainID: p}