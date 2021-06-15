
import sys

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

"""
This script accepts the name of a single simulation for which .edges, .points files exist in the working directory
and plots:
1) Wing with edges colored according to type
2) Wing with edges shadowed according to tension
3) Angle of edge vs tension, colors according to edge type (Note: tension of border edges is set to max tension of non-border edges)

To plot only a subset of the plots pass a second argument which has letters:
c : plot 1
t : plot 2
p : plot 3

Example:
python plot_wing_according_to_edge_type simulation1 -ctp

It doesn't save the plots. 
"""

basename = sys.argv[1] #'bud2_moved_0'
plot_types = sys.argv[2] if (len(sys.argv) >2) else "-ctp"
e=pd.read_csv(basename + '.edges', sep="\t") #file with edges, pasted from .out file (need to prepare this first)
pf=open(basename + '.points', "r")
p=int(pf.readline()) #discard first line
p=[[float(j) if '.' in j else int(j) for j in i.split("\t")] for i in pf.read().split("\n") if i != ""]
p=pd.DataFrame(p,columns=["x", "y", "ind", "movable"])

#p=pd.read_csv(basename + '.points', sep="\t") #file with points (like .points but with header (x	y	ind	movable); need to prepare this first)

colors = ['black', 'green', 'yellow', 'orange', 'blue', 'purple','red', 'peru', 'peru']

if('c' in plot_types):
    for i in range(e.shape[0]):
        ind = e["vertices"][i].split(',')
        paux = p[(p['ind'] == int(ind[0])) | (p['ind'] == int(ind[1]))]
        cc = colors[e["type"][i]] if(e["type"][i] < len(colors)) else colors[-1]
        fig = plt.plot(paux['x'], paux['y'], color=cc)
    plt.show()

mint = np.min(e['tension'][e["type"] != 4])
maxt = np.max(e['tension'][e["type"] != 4]) #borders always have bigger tension
e['tension'][e["type"] == 4] = maxt
normt = (e['tension'] - mint)/(maxt - mint) + 0.1
normt = normt/np.max(normt)

if('t' in plot_types):
    for i in range(e.shape[0]):
        ind = e["vertices"][i].split(',')
        paux = p[(p['ind'] == int(ind[0])) | (p['ind'] == int(ind[1]))]
        cc = colors[e["type"][i]] if(e["type"][i] < len(colors)) else colors[-1]       
        fig=plt.plot(paux['x'], paux['y'], color = 'black', alpha = normt[i] )
    plt.show()


if('p' in plot_types):
    angles=[]
    colorsthis = []
    for i in range(e.shape[0]):
        ind = e["vertices"][i].split(',')
        paux = p[(p['ind'] == int(ind[0])) | (p['ind'] == int(ind[1]))]
        angles.append(abs(np.arctan2(paux.y.iloc[0] - paux.y.iloc[1], paux.x.iloc[0] - paux.x.iloc[1])*180/np.pi))
        colorsthis.append(colors[e["type"][i]] if(e["type"][i] < len(colors)) else colors[-1])
    plt.scatter(angles, e["tension"], color=colorsthis)
    plt.axvline(x=90)
    plt.xlabel("Edge angle in degrees")
    plt.ylabel("Edge tension (raw value)")
    plt.show()

