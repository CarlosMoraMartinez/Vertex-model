import os
import sys
import re

import pandas as pd
import matplotlib.pyplot as plt

STARTCHAR = '>'
SEPRE = ";"
FIELDRE = "="
NAMERE = ':'
def readOutputFile(fname):
    divisions = []
    t1 = []
    t2 = []
    i=0
    with open(fname) as f:
        for s in f:
            if(not s.startswith(STARTCHAR)):
                continue
            seps = re.split(NAMERE, "".join(s.split()))
            if(seps[0].startswith(">T1")):
                t1.append(dict([val.split(FIELDRE) for val in seps[1].split(SEPRE)]))
            elif(seps[0].startswith(">DIV")):
                divisions.append(dict([val.split(FIELDRE) for val in seps[1].split(SEPRE)]))
            elif(seps[0].startswith(">T2")):
                t2.append(dict([val.split(FIELDRE) for val in seps[1].split(SEPRE)]))
            i+=1
            #if(i>100000):
            #    break
    return (pd.DataFrame(t1).astype(float), pd.DataFrame(divisions).astype(float), pd.DataFrame(t2).astype(float))





fname = sys.argv[1]#'w2F_e18f190.out'#'ensemble17_1-wing2E.out'
t1, d, t2 = readOutputFile(fname)
maxt=int(sys.argv[2])

fig = plt.figure()
ax = plt.subplot(111)


ax.scatter(t1['mov.accepted'], t1['centroid_x'], color = "blue", label = "T1");plt.xlim(0,maxt);plt.ylim(40.0,480.0);
ax.scatter(t2['mov.accepted'], t2['centroid_x'], color = "green", label = "T2");plt.xlim(0,maxt);plt.ylim(40.0,480.0);
ax.scatter(d['movesaccepted'], d['centroid_x'], color = "red", label = "Divisions");plt.xlim(0,maxt);plt.ylim(40.0,480.0);

figtitle=sys.argv[1].strip('.out')

chartBox = ax.get_position()
ax.set_position([chartBox.x0, chartBox.y0 + chartBox.height*0.2, chartBox.width, chartBox.height*0.8])
ax.set_ylabel('P-D axis (x coordinate)')
ax.set_xlabel('time (accepted movements)')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), shadow=True, ncol=3)
plt.title(figtitle)
#plt.show()
plt.savefig(figtitle + '.pdf')


