class ChainIDManager:
    def __init__(self):
        U=[chr(i) for i in range(ord('A'),ord('A')+26)]
        L=[chr(i) for i in range(ord('a'),ord('a')+26)]
        D=[str(i) for i in range(10)]
        self.OrderedSupply=U+L+D
        self.Used=set()
    
    def generate_next_map(self,chainIDs):
        assert len(chainIDs)<=len(self.OrderedSupply),f'Not enough available chainIDs'
        myMap={}
        for c in chainIDs:
            if c in self.OrderedSupply:
                self.OrderedSupply.remove(c)
                self.Used.add(c)
        for c in chainIDs:
            p=self.OrderedSupply.pop(0)
            self.Used.add(p)
            myMap[c]=p
        return myMap
    
    def cleavage_daughter_chainID(self,chainID):
        assert 1<=len(self.OrderedSupply),f'Not enough available chainIDs'
        assert chainID in self.OrderedSupply,f'Parent chain {chainID} was never claimed'
        p=self.OrderedSupply.pop(0)
        self.Used.add(p)
        return {chainID: p}