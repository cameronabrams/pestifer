# Author: Cameron F. Abrams <cfa22@drexel.edu>.
""" Chain ID manager -- make sure all chains have unique IDs and ID maps are tracked
"""

import logging
logger=logging.getLogger(__name__)

class ChainIDManager:
    """A class for managing chainIDs.  Allows for users to request changes to chainIDs, assign specific chainIDs to transformed subunits, and prevents duplicate usage.
    
    Methods
    -------
    __init__(format, transform_reserves)
      initializes the repository of chainIDs depending on format (PDB or mmCIF)
       
    check(chainID)
      checks to see if proposed chainID is available, and if so, registers it and returns it.  If not, counters proposal with next available chainID, registers that and returns it.
    
    unregister_chain(chainID)
      Unregisters the chainID
      
    next_unused_chainID()
      Registers the next available non-reserved chainID and returns it
    
    next_reserved_chainID(key)
      Registers and returns the next available chainID reserved to key

    generate_next_map()
      Given the original chainIDs, returns a dictionary mapping each 
      to a new chainID
    
    thru_map()
      Generates an identity map
      
    cleavage_daughter_chainID(chainID)
      returns a single-entry dictionary mapping the chainID to the
      next available chainID
    """
    def __init__(self,format='PDB',remap={},transform_reserves={}):
        logger.debug(f'New chainIDmanager, format {format}')
        self.format=format
        self.remap=remap
        self.transform_reserves=transform_reserves
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
        for k,v in self.transform_reserves.items():
            assert k in self.Unused,f'Error: transform-reserved map key chainID {k} is not available in the Unused set'
            for mc in v: # only put image chainIDs in ReservedUnused
                assert mc in self.Unused,f'Error: transform-reserved map chainID {mc} is not available in the Unused set'
                self.Unused.remove(mc)
                self.ReservedUnused.append(mc)
                # at this point, the keys in transform_reserves are also elements in Unused
        self.Used=set()

    def sandbag(self,i_chainIDs):
        """ moves all chainIDs in i_chainIDs to end of Unused so they are not popped too early."""
        logger.debug(f'Sandbagging: unused chains: {self.Unused}')
        assert all([(x in self.Unused) or (x in self.ReservedUnused) for x in i_chainIDs]),f'Cannot sandbag since at least one initial chainID is not found in Unused or ReservedUnusued -- bug!'
        for c in i_chainIDs:
            if c in self.Unused:
                logger.debug(f'Sandbagging unused chainID {c}')
                self.Unused.remove(c)
                self.Unused.append(c)
            else:
                logger.debug(f'Not sandbagging chainID {c} since it is reserved for a transform')

    def check(self,proposed_chainID):
        hold_chainID=proposed_chainID
        if hold_chainID in self.remap.values():
            logger.debug(f'proposed chainID {hold_chainID} is already reserved as a user-map')
            if hold_chainID in self.remap.keys():
                logger.debug('but it is also a key, so...')
                logger.debug(f'proposed chainID {hold_chainID} is user-mapped to chainID {self.remap[hold_chainID]}')
                hold_chainID=self.remap[hold_chainID]    
        elif hold_chainID in self.remap.keys():
            logger.debug(f'proposed chainID {hold_chainID} is user-mapped to chainID {self.remap[hold_chainID]}')
            hold_chainID=self.remap[hold_chainID]
        else:
            logger.debug(f'proposed chainID {hold_chainID} does not collide with any predefined chainID user maps')
        # hold_chainID has not yet been claimed
        #         
        # caller may not propose a chainID that the user has reserved for transforms
        if hold_chainID in self.ReservedUnused:
            logger.debug(f'chainID {hold_chainID} (orig {proposed_chainID}) is reserved for the product of an asymmetric-unit transform')
            hold_chainID=self.next_unused_chainID()
            logger.debug(f'counter-proposed chainID is {hold_chainID}')
        elif hold_chainID in self.Used:
            logger.debug(f'chainID {hold_chainID} (orig {proposed_chainID}) is already used')
            hold_chainID=self.next_unused_chainID()
            logger.debug(f'counter-proposed chainID is {hold_chainID}')
        else:
            logger.debug(f'registering chainID {hold_chainID}')
            self.Unused.remove(hold_chainID)
            self.Used.add(hold_chainID)

        return hold_chainID

    def unregister_chain(self,chainID):
        """ Unregister the single chainID with the manager 
        
        Paramter:
        ---------
        chainID : str
           1- or 2-byte chainID to unregister
        """
        self.Used.remove(chainID)
        is_reserved=any([chainID in v for v in self.transform_reserves.values()])
        if is_reserved:
            self.ReservedUnused.append(chainID)
        else:
            self.Unused.append(chainID)

    def next_unused_chainID(self):
        p=self.Unused.pop(0)
        self.Used.add(p)
        return p
    
    def next_reserved_chainID(self,key):
        assert key in self.transform_reserves,f'Key chainID {key} not a key in the prescribed reserve maps'
        mapsto=self.transform_reserves[key]
        assert any([x in self.ReservedUnused for x in mapsto]),f'mapped key {key} has no daughters available'
        for c in mapsto:
            if c in self.ReservedUnused:
                self.ReservedUnused.remove(c)
                self.Used.add(c)
                return c
        return None
        
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
            if c in self.transform_reserves:
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
    
    # def cleavage_daughter_chainID(self,chainID):
    #     assert 1<=len(self.Unused),f'Not enough available chainIDs'
    #     assert chainID in self.Unused,f'Parent chain {chainID} was never claimed'
    #     return {chainID: self.next_unused_chainID()}