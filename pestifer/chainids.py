class ChainIDManager:
    def __init__(self):
        L=[chr(i) for i in range(ord('a'),ord('a')+26)]
        U=[chr(i) for i in range(ord('A'),ord('A')+26)]
        D=[str(i) for i in range(10)]
        self.OrderedSupply=U+L+D
        self.Used=set()
    
    def generate_map(self,chainIDs):
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