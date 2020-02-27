


import pandas as pd
import matplotlib.pyplot as plt
e=pd.read_csv("ee.csv", sep="\t") #file with edges, pasted from .out file (need to prepare this first)
p=pd.read_csv("pp.csv", sep="\t") #file with points (like .points but with header; need to prepare this first)

colors = ['black', 'green', 'yellow', 'orange', 'blue', 'purple','red', 'peru']
for i in range(e.shape[0]):
    ind = e["vertices"][i].split(',')
    paux = p[(p['ind'] == int(ind[0])) | (p['ind'] == int(ind[1]))]
    cc = colors[e["type"][i]] if(e["type"][i] < len(colors)) else colors[-1]       
    fig=plt.plot(paux['x'], paux['y'], color = cc)

plt.show()
