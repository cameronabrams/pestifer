# Author: Cameron F. Abrams <cfa22@drexel.edu>.
""" 
A chain ID manager makes sure all chains have unique IDs and ID maps are tracked.
"""

import logging
logger=logging.getLogger(__name__)

class ChainIDManager:
    """
    A class for managing chainIDs.  Allows for users to request changes to chainIDs, assign specific chainIDs to transformed subunits, and prevents duplicate usage.
    
    Parameters
    ----------
    format : str, optional
        The format of the chain IDs, either 'PDB' or 'mmCIF'.
        Default is 'PDB'.
    remap : dict, optional
        A dictionary mapping user-defined chain IDs to their desired replacements.
        Default is an empty dictionary.
    transform_reserves : dict, optional
        A dictionary mapping chain IDs to lists of reserved chain IDs for transformed subunits.
        This allows for specific chain IDs to be reserved for the products of asymmetric unit transformations.
        Default is an empty dictionary.
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
        """ 
        moves all chainIDs in i_chainIDs to end of Unused so they are not popped too early.

        Parameters
        ----------
        i_chainIDs : list of str
            List of chainIDs to sandbag. These chainIDs will be moved to the end of the Unused list, ensuring they are not used until explicitly requested.
        """
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
        """
        Check if the proposed chainID is available and return a valid chainID.
        If the proposed chainID is already in use or reserved, a new chainID will be proposed.

        Parameters
        ----------
        proposed_chainID : str
            The chainID proposed by the user. This can be a 1- or 2-byte string, depending on the format of the ChainIDManager.
        """
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
        """ 
        Unregister the single chainID with the manager 
        
        Parameters
        ----------
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
        """ 
        Returns the next unused chainID from the Unused list.
        If there are no unused chainIDs, this will raise an error.
        """
        p=self.Unused.pop(0)
        self.Used.add(p)
        return p
    
    def next_reserved_chainID(self,key):
        """
        Returns the next available chainID from the ReservedUnused list that is mapped to the given key.
        If the key is not in the transform_reserves, an error is raised.
        If no reserved chainIDs are available for the key, an error is raised.
        
        Parameters
        ----------
        key : str
            The chainID key for which to find the next available reserved chainID.
        """
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
        """ 
        Generate a map identifying new chainIDs for each existing chainID
        in the list chainIDs.  If active_chains is not empty, then 
        it must list a subset of the chainIDs considered actually active
        
        Parameters
        ----------
        chainIDs : list of str
            List of chainIDs for which to generate a map of next available chainIDs.
        active_chains : list of str, optional
            List of chainIDs that are currently active. If provided, only these chainIDs will be
            considered for the next available chainIDs. If not provided, all chainIDs in chainIDs
            will be considered active.
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
        """
        Generate a map that maps each chainID in chainIDs to itself, but only for the
        chainIDs that are considered active. If active_chains is provided, only those
        chainIDs will be included in the map. If active_chains is empty, all chainIDs in chainIDs will be considered active.

        Parameters
        ----------
        chainIDs : list of str
            List of chainIDs to include in the map.
        active_chains : list of str, optional
            List of chainIDs that are currently active. If provided, only these chainIDs will be
            included in the map. If not provided, all chainIDs in chainIDs will be considered active.
        """
        activeChainIDs=chainIDs.copy()
        if active_chains:
            inactive_chains=[x for x in activeChainIDs if not x in active_chains]
            for i in inactive_chains:
                activeChainIDs.remove(i)
        logger.debug(f'generating thru_map from {activeChainIDs} with actives {activeChainIDs}')
        return {c:c for c in activeChainIDs}
    