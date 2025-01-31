import matplotlib.pyplot as plt
import matplotlib.cm as cm

cmap=cm.get_cmap('inferno')

fig,ax=plt.subplots(1,1,figsize=(8,6))
ax.set_xlabel("time step")
ax.set_ylabel("Volume (cubic Angstrom)")

for r in range(1,6):
    ts=[]
    v=[]
    for s in range(0,5):
        fp=open('rep{}/run4_stage{}.log'.format(r,s),'r')
        for l in fp:
            tokens=l.split()
            if len(tokens)>0 and tokens[0]=='ENERGY:':
                ts.append(int(tokens[1]))
                v.append(float(tokens[18]))
            elif len(tokens)>0 and tokens[0]=='ETITLE:':
                xl=tokens[1]
                yl=tokens[18]
        fp.close()
    ax.plot(ts,v,label='rep{}'.format(r),color=cmap((r/6)))

ax.legend
plt.savefig('myplot.png')
