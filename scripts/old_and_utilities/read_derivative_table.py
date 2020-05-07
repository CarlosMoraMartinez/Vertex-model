import sys
import os
import re

import pandas as pd
import numpy as np

newdir="derivatives_plots"
try:
    os.mkdir(newdir)
except:
    print("Unable to create directory")

try:
    fname = sys.argv[1]
except:
    fname = "eulertab.csv"

f = open(fname, "r")
raw=[[float(k.split("=")[1]) for k in re.split(re.compile("[:,] "), i)] for i in f]
tab= pd.DataFrame(columns=["ind", "x1","y1","x2","y2","x3","y3","dx","dy"], data=raw)
tab = tab.astype({"ind": int}) 
f.close()    

import matplotlib.pyplot as plt

vars2plot = ["x1","y1","x2","y2","x3","y3","dx","dy"]
for v in vars2plot:
    for i in pd.unique(tab.ind):
        aux = tab.loc[tab.ind == i]
        plt.plot(range(aux.shape[0]), aux[v])
    plt.yscale('symlog', linthreshy=0.015)
    plt.savefig("%s/%s_%s.pdf"%(newdir, fname.split(".")[0], v))
    plt.clf()

