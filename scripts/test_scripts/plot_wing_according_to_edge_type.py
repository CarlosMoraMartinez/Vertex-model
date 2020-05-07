


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

basename = 'bud2_moved_0'
e=pd.read_csv(basename + '.edges', sep="\t") #file with edges, pasted from .out file (need to prepare this first)
p=pd.read_csv(basename + '.points', sep="\t") #file with points (like .points but with header (x	y	ind	movable); need to prepare this first)

colors = ['black', 'green', 'yellow', 'orange', 'blue', 'purple','red', 'peru']
for i in range(e.shape[0]):
    ind = e["vertices"][i].split(',')
    paux = p[(p['ind'] == int(ind[0])) | (p['ind'] == int(ind[1]))]
    cc = colors[e["type"][i]] if(e["type"][i] < len(colors)) else colors[-1]       
    fig=plt.plot(paux['x'], paux['y'], color = cc)

plt.show()

mint = np.min(e['tension'][e["type"] != 4])
maxt = np.max(e['tension'][e["type"] != 4]) #borders always have bigger tension
e['tension'][e["type"] == 4] = maxt
normt = (e['tension'] - mint)/(maxt - mint) + 0.1
normt = normt/np.max(normt)


for i in range(e.shape[0]):
    ind = e["vertices"][i].split(',')
    paux = p[(p['ind'] == int(ind[0])) | (p['ind'] == int(ind[1]))]
    cc = colors[e["type"][i]] if(e["type"][i] < len(colors)) else colors[-1]       
    fig=plt.plot(paux['x'], paux['y'], color = 'black', alpha = normt[i] )

plt.show()
