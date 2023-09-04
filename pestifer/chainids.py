import logging
logger=logging.getLogger(__name__)

class ChainIDManager:
    def __init__(self,format='PDB'):
        logger.debug(f'New chainIDmanager, format {format}')
        self.format=format
        U=[chr(i) for i in range(ord('A'),ord('A')+26)]
        if format=='PDB':
            L=[chr(i) for i in range(ord('a'),ord('a')+26)]
            D=[str(i) for i in range(10)]
            self.OrderedSupply=U+L+D
        elif format=='mmCIF':
            self.OrderedSupply=[]
            for a in ['']+U:
                for b in U:
                    self.OrderedSupply.append(b+a)
        self.Used=set()

    def register_asymm_chains(self,chainIDs):
        logger.debug(f'Registering asymm chains {chainIDs}')
        for c in chainIDs:
            self.Used.add(self.OrderedSupply.remove(c))

    def next_chain(self):
        p=self.OrderedSupply.pop(0)
        self.Used.add(p)
        return p

    def generate_next_map(self,chainIDs,active_chains=[]):
        assert len(chainIDs)<=len(self.OrderedSupply),f'Not enough available chainIDs'
        myMap={}
        activeChainIDs=chainIDs.copy()
        # logger.debug(f'generating next map from {activeChainIDs} with actives {active_chains}')
        if active_chains:
            inactive_chains=[x for x in activeChainIDs if not x in active_chains]
            for i in inactive_chains:
                activeChainIDs.remove(i)
        for c in activeChainIDs:
            if c in self.OrderedSupply: # should never happen
                self.OrderedSupply.remove(c)
                self.Used.add(c)
        for c in activeChainIDs:
            myMap[c]=self.next_chain()
        logger.debug(f'generated next map: {myMap}')
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
        assert 1<=len(self.OrderedSupply),f'Not enough available chainIDs'
        assert chainID in self.OrderedSupply,f'Parent chain {chainID} was never claimed'
        p=self.OrderedSupply.pop(0)
        self.Used.add(p)
        return {chainID: p}