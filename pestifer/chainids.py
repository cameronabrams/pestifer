import logging
logger=logging.getLogger(__name__)

class ChainIDManager:
    def __init__(self,format='PDB'):
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

    def generate_next_map(self,chainIDs,active_chains=[]):
        assert len(chainIDs)<=len(self.OrderedSupply),f'Not enough available chainIDs'
        myMap={}
        activeChainIDs=chainIDs.copy()
        if active_chains:
            inactive_chains=[x for x in activeChainIDs if not x in active_chains]
            for i in inactive_chains:
                activeChainIDs.remove(i)
        for c in activeChainIDs:
            if c in self.OrderedSupply:
                self.OrderedSupply.remove(c)
                self.Used.add(c)
        for c in activeChainIDs:
            p=self.OrderedSupply.pop(0)
            self.Used.add(p)
            myMap[c]=p
        return myMap
    
    def thru_map(self,chainIDs,active_chains=[]):
        activeChainIDs=chainIDs.copy()
        logger.debug(f'generating thru_map from {activeChainIDs} with actives {active_chains}')
        if active_chains:
            inactive_chains=[x for x in activeChainIDs if not x in active_chains]
            for i in inactive_chains:
                activeChainIDs.remove(i)
        return {c:c for c in activeChainIDs}
    
    def cleavage_daughter_chainID(self,chainID):
        assert 1<=len(self.OrderedSupply),f'Not enough available chainIDs'
        assert chainID in self.OrderedSupply,f'Parent chain {chainID} was never claimed'
        p=self.OrderedSupply.pop(0)
        self.Used.add(p)
        return {chainID: p}