from operator import index
import sys
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

LM_NAMES = ['1', '2', '3', '4', '12']
colors = ['black', 'green', 'yellow', 'orange', 'blue', 'purple','red', 'peru']

#Edge types
BORDER_TYPE = 4
STRING_TYPE = 7
CUTICLE_TYPE = 8
VEIN_TYPE = 6 
#Cell types
BLADE_TYPE = 0
HINGE_TYPE = 1
BLADE_VEIN_TYPE = 2
HINGE_VEIN_TYPE = 3

Y_PROPORTION = 0.25
EMPTY_CONNECTION = -999

def loadData(basename):
    e=pd.read_csv(basename + '.edges', sep="\t") #file with edges, pasted from .out file (need to prepare this first)
    c=pd.read_csv(basename + '.celltab', sep="\t")
    with open(basename + '.points', 'r') as pf:
        n = int(pf.readline())
        p = [[float(i) if '.' in i else int(i) for i in l.strip('\n').split('\t')] for l in pf]
    p = pd.DataFrame( dict(zip(['x', 'y', 'ind', 'movable'], zip(*p))) ).set_index("ind", drop=False)
    e = e.set_index('ind', drop=False)
    c = c.set_index('ind', drop=False)
    p = p.set_index('ind', drop=False)
    return (p, c, e)

def getMinAndMax(p, border):
    xmin = p.loc[0].x
    ymin = p.loc[0].y
    xmax, ymax = (xmin, ymin)
    for i in border:
        xmin = p.loc[i].x if p.loc[i].x < xmin else xmin
        ymin = p.loc[i].y if p.loc[i].y < ymin else ymin
        xmax = p.loc[i].x if p.loc[i].x > xmax else xmax
        ymax = p.loc[i].y if p.loc[i].y > ymax else ymax
    return (xmin, xmax, ymin, ymax)

def getBorder(v, c, e):
    eborder = e[e['type'] == BORDER_TYPE]
    edgev = [[int(i) for i in x.split(',') if i != ''] + [y] for x, y in zip(eborder['vertices'], eborder['ind'])]
    border = []
    border_e = []
    start = edgev.pop()
    border.append(start[0])
    border.append(start[1])
    border_e.append(start[2])
    current = start[1]
    while(edgev):
        for i in range(len(edgev)):
            if(current == edgev[i][0]):
                e = edgev.pop(i)
                border.append(e[1])
                current = e[1]
                border_e.append(e[2])
                break
            elif(current == edgev[i][1]):
                e = edgev.pop(i)
                border.append(e[0])
                current = e[0]
                border_e.append(e[2])
                break
        else:
            print("Error: border not continuous")
            return []
    border, border_e = map(np.array, [border, border_e])
    xmin, xmax, ymin, ymax = getMinAndMax(v, border)
    closest = -1
    xmin_local = xmax
    roll = 0
    for i, vb in enumerate(border):
        if(v.loc[vb].y - ymin < Y_PROPORTION*(ymax - ymin)):
            if(v.loc[vb].x < xmin_local):
                xmin_local, closest, roll = (v.loc[vb].x, vb, i)
    border = np.roll(border, -1*roll)
    border_e = np.roll(border_e, -1*roll)
    return (border, border_e) #border_e is the edge needed to get to the NEXT vertex (eg. border_e[0] is needed to go from border[0] to border[1])

def plotOrderedBorder(p, border):
    xmin, xmax, ymin, ymax = getMinAndMax(p, border)
    plt.plot(p.loc[border].x, p.loc[border].y)
    plt.scatter(p.loc[border[0]].x, p.loc[border[0]].y, color="black")
    plt.scatter(p.loc[border[10]].x, p.loc[border[10]].y, color="red")
    plt.text(p.loc[border[0]].x, p.loc[border[0]].y, str(border[0]), color= "black")
    plt.text(p.loc[border[10]].x, p.loc[border[10]].y, str(border[10]), color= "red")
    plt.hlines([ymin, ymax], 0, xmax*1.1, color="black")
    plt.vlines([xmin, xmax], 0, ymax*1.1, color="black")
    plt.hlines(ymin + (ymax - ymin)*Y_PROPORTION, 0, xmax*1.1, color="red")
    plt.show()

def getLandmarks1to4(p, c, e, border, border_e):
    veins_v = [[], [], [], []]
    veins_e = [[], [], [], []]
    current_vein = -1
    in_vein = False
    #get landmarks 1-4
    for i, be in enumerate(border_e):
        cells = e.loc[be].cells.split(',')
        cell = int(cells[0]) if int(cells[1]) == EMPTY_CONNECTION else int(cells[1])
        if(c.loc[cell].type == BLADE_VEIN_TYPE):
            if(not in_vein):
                in_vein = True
                current_vein +=1
                veins_e[current_vein].append(be)
                veins_v[current_vein].append(border[i])
                veins_v[current_vein].append(border[i + 1])
            else:
                veins_e[current_vein].append(be)
                veins_v[current_vein].append(border[i + 1])
        else:
            in_vein = False
    landmarks = [list(map(np.mean, zip(*[(p.loc[i].x, p.loc[i].y) for i in vein]))) for vein in veins_v] #should return [[x0, y0], [x1, y1], ...]
    landmarks = sorted(landmarks, key = lambda x : x[1], reverse = True) #
    return landmarks

def getLandmark12_old(p, c, e, border, border_e):
    eveins = e[e['type'] == VEIN_TYPE]
    eveins = eveins.assign(v1= [int(x.split(',')[0]) for x in eveins.vertices], v2= [int(x.split(',')[1]) for x in eveins.vertices] )
    interedges = [xind for x, xind in [[(c.loc[int(i.cells.split(',')[0])].type, c.loc[int(i.cells.split(',')[1])].type), iind] for iind, i in eveins.iterrows()] if(BLADE_VEIN_TYPE in x and HINGE_VEIN_TYPE in x)]
    order_ie = np.argsort([0.5*(p.loc[eveins.loc[i].v1].y  + p.loc[eveins.loc[i].v2].y) for i in interedges])
    order_ie = np.flip(order_ie, axis = 0)
    interedges = [x for _,x in sorted(zip(order_ie, interedges))]

    current_edge = interedges.pop()
    vein_v = set([eveins.loc[current_edge].v1, eveins.loc[current_edge].v2]) #Only want vein with lowest 'y' coordinate (vein 4)
    vein_e = [current_edge]
    while interedges:
        current_edge = interedges.pop()
        if(eveins.loc[current_edge].v1 in vein_v):
            vein_e.append(current_edge)
            vein_v.add(eveins.loc[current_edge].v2)
        elif(eveins.loc[current_edge].v2 in vein_v):
            vein_e.append(current_edge)
            vein_v.add(eveins.loc[current_edge].v1)
    res = list(map(np.mean, zip(*[[p.loc[i].x, p.loc[i].y] for i in vein_v])))
    return res


def getLandmark12(p, c, e, border, border_e):
    #eveins = e[e['type'] == VEIN_TYPE]
    eveins = e.copy()
    eveins = eveins.assign(v1= [int(x.split(',')[0]) for x in eveins.vertices], v2= [int(x.split(',')[1]) for x in eveins.vertices],
     c1= [int(x.split(',')[0]) for x in eveins.cells], c2= [int(x.split(',')[1]) for x in eveins.cells])
    e_celltypes = [[c.loc[a].type, c.loc[b].type] if(a in c.index and b in c.index) else [EMPTY_CONNECTION] for a, b in zip(eveins.c1, eveins.c2)]
    eveins = eveins.assign(interedge = [BLADE_VEIN_TYPE in x and HINGE_VEIN_TYPE in x for x in e_celltypes])
    interedge_tab = eveins.iloc[eveins.interedge.tolist()]
    interedge_tab = interedge_tab.assign(y_position = [0.5*(p.loc[a].y + p.loc[b].y) for a, b in zip(interedge_tab.v1, interedge_tab.v2)])
    order_ie = np.argsort(interedge_tab.y_position.tolist())
    order_ie = np.flip(order_ie, axis = 0)
    interedge_tab = interedge_tab.iloc[order_ie]
    interedges=interedge_tab.ind.tolist()
    #interedges = [x for _,x in sorted(zip(order_ie, interedges))]
    current_edge = interedges.pop()
    vein_v = set([interedge_tab.loc[current_edge].v1, interedge_tab.loc[current_edge].v2]) #Only want vein with lowest 'y' coordinate (vein 4)
    vein_e = [current_edge]
    while interedges:
        current_edge = interedges.pop()
        if(interedge_tab.loc[current_edge].v1 in vein_v):
            vein_e.append(current_edge)
            vein_v.add(interedge_tab.loc[current_edge].v2)
        elif(interedge_tab.loc[current_edge].v2 in vein_v):
            vein_e.append(current_edge)
            vein_v.add(interedge_tab.loc[current_edge].v1)
    res = list(map(np.mean, zip(*[[p.loc[i].x, p.loc[i].y] for i in vein_v])))
    return res


def getLandmarks(p, c, e, border, border_e):
    lms = getLandmarks1to4(p, c, e, border, border_e) + [getLandmark12(p, c, e, border, border_e)]
    return dict(zip(LM_NAMES, lms))


def plotLandmarks(p, c, e, lms, showAllCells=False, show=False, name = 'test'):
    if(not showAllCells): #makes plot faster
        e = e[(e['type']  != HINGE_TYPE) & (e['type'] != BLADE_TYPE)]
    colors = ['black', 'green', 'yellow', 'orange', 'blue', 'purple','red', 'peru']
    for iix, i in e.iterrows():
        indices = [int(x) for x in i.vertices.split(',') if x != '' and x != str(EMPTY_CONNECTION)]
        paux = p.loc[indices]
        cc = colors[i.type] if(i.type < len(colors)) else colors[-1]
        fig=plt.plot(paux['x'], paux['y'], color = cc)
    for n, l in lms.items():
        plt.scatter(l[0], l[1], color = "blue")
        plt.text(l[0], l[1], n, color = "blue")
    if(show):
        plt.show()
    else:
        plt.savefig(name + '_landmarks.pdf')

def plotLandmarks2(p, c, e, lms, showAllCells=2, show=False, name = 'test'):
    if(showAllCells == 1): #makes plot faster
        e = e[(e['type']  != HINGE_TYPE) & (e['type'] != BLADE_TYPE)]
    if(showAllCells < 2): ##Plots edges
        for iix, i in e.iterrows():
            indices = [int(x) for x in i.vertices.split(',') if x != '' and x != str(EMPTY_CONNECTION)]
            paux = p.loc[indices]
            cc = colors[i.type] if(i.type < len(colors)) else colors[-1]
            fig=plt.plot(paux['x'], paux['y'], color = cc)
    elif(showAllCells == 2): #dotplot of cells
        c2 = c.loc[c.centroid_x > 0].loc[c.centroid_y > 0]
        plt.scatter(c2.centroid_x, c2.centroid_y, color = [colors[i] for i in c2.type])
    for n, l in lms.items():
        plt.scatter(l[0], l[1], color = "blue")
        plt.text(l[0], l[1], n, color = "blue")
    if(show):
        plt.show()
    else:
        plt.savefig(name + '_landmarks.pdf')

def processWing(name, plotToFile=2):
    p, c, e = loadData(name)
    border, border_e = getBorder(p, c, e)
    #plotOrderedBorder(p, border)
    lms = getLandmarks(p, c, e, border, border_e)
    if(plotToFile):
        plotLandmarks2(p, c, e, lms, plotToFile, False, name)
    lms.setdefault('name', name)
    plt.clf()
    return lms


def writeLandmarksToCsv(landmarks, basename = 'all'):
    df = pd.DataFrame(landmarks)
    df.to_csv(basename + '_landmarks.csv', sep="\t", index=False)

def main():
    #os.chdir('C:\\Users\\Carlos\\Desktop\\scriptfindlandmarks')
    #name = 'wing2E_moved_4'
    print("WARNING: This script probably does not work well if you are using a cellular cuticle. Remove cuticle vertices or cells beforehand")
    basename = sys.argv[1]
    if(len(sys.argv) > 2):
        plotLms = int(sys.argv[2]) #> 0
        if (plotLms):
            print("Plotting ON")
        else:
            print("Plotting OFF")
    else:
        plotLms = True
        print("Plotting ON")

    files = set([f.split('.')[0] for f in os.listdir() if basename in f])
    landmarks = []
    for name in files:
        try:
            print("Calculating landmarks for ", name)
            landmarks.append(processWing(name, plotLms))
            print("    Calculated landmarks for ", name)
        except:
            print("ERROR calculating landmarks for", name)
    writeLandmarksToCsv(landmarks, basename)

if __name__ == '__main__':
    main()

