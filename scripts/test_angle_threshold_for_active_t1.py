
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import copy as cp

spr = 'etournay1_unmoveX3_moved_0.sprtab'
d = pd.read_table(spr)
d2 = cp.copy(d)
d2.x_static = -1*d2.x_static
d2.x_movable = -1*d2.x_movable
d=d.append(d2)
d2 = cp.copy(d)
d2.y_static = -1*d2.y_static
d2.y_movable = -1*d2.y_movable
d=d.append(d2)
d["rown"] = range(d.shape[0])
d=d.set_index("rown")


d=pd.DataFrame({"x_static":np.zeros(400), "y_static":np.zeros(400), "x_movable" : np.linspace(-1,1,20).tolist()*20, "y_movable" : sorted(np.linspace(-1,1,20).tolist()*20)})
maxan = 30
minan = 0
for i in range(d.shape[0]):
    angle = 180*np.arctan2(d.loc[i].x_static - d.loc[i].x_movable, d.loc[i].y_static - d.loc[i].y_movable)/np.pi
    print(i, ": ", angle)
    color = "red"  if(abs(angle)%180 >= minan and abs(angle)%180 <= maxan) else "black"   
    plt.plot([d.loc[i].x_static, d.loc[i].x_movable], [d.loc[i].y_static, d.loc[i].y_movable], color=color)

plt.show()

name = ""
pfile = name+'.points'
efile = name+'.edges'

f=sys.open(pfile)
npoints = int(f.readline())
points = []
for l in f:
  points.append(l.split('\t'))



