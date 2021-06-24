from operator import itemgetter
import sys
import os
import argparse

import numpy as np
import scipy.spatial as spatial
import matplotlib.pyplot as plt 
from matplotlib.widgets import Slider, Button, Lasso, LassoSelector, TextBox, PolygonSelector
import shapely
from shapely.geometry import Polygon, Point
import pandas as pd

BLUE = '#6699cc'
GRAY = '#999999'
RED = '#ff3333'
WHITE = '#ffffff'
BLACK = '#000000'

bladetype = 0
hingetype = 1
veintype = 2
veinhinge = 3
BORDER_TYPE = 4

wingcols = ['blue', 'green', 'black', 'gray']
springcols = ['red', 'yellow', 'purple', 'pink', 'brown']


class wing:
    def __init__(self, name):
        self.name = name
        self.numPoints = 0
        self.pointsList = {}
        self.numCells = 0
        self.plot_cell_types = True
        self.numVerticesCell = []
        self.polygonList = []
        self.cellvertList = []
        self.celltypes = []
        self.numsprings = 0
        self.sprList = []
        self.celltab = None
        self.edgetab = None
        self.fig = None
        self.ax = None
        try:
            self.readData()
        except:
            print("Unable to read data for wing %s"%(self.name))
    def readData(self):
        try:
            self.readPointsFile()
        except:
            print("Error reading .points file")
        try:
            self.readCellsFile()
        except:
            print("Error reading .cells file")
        try:
            self.readCellsTabFile()
        except:
            print("No .celltab file")
        try:
            self.readEdgeTabFile()
        except:
            print("No .edges file")
        self.readSprings()
    def setAx(self, fig, ax, lim=[]):  
        self.fig=fig
        self.ax=ax
        self.setLimits(lim)
    def plot(self):
        from matplotlib.collections import PolyCollection 
        self.setLimits() 
        self.fig, self.ax = plt.subplots()
        plt.xlim(self.limits[0], self.limits[2])
        plt.ylim(self.limits[1], self.limits[3])     
        pc = PolyCollection(self.polygonList, facecolors= [wingcols[self.celltypes[j]] for j in range(self.numCells)], alpha = 0.6)
        pc.set_edgecolors("white")
        self.ax.add_collection(pc)
        for s in self.sprList:
            spcolor = springcols[int(s[2])] if len(s) > 2 else RED
            self.ax.plot([self.pointsList[s[0]][0], self.pointsList[s[1]][0]], [self.pointsList[s[0]][1], self.pointsList[s[1]][1]], color = spcolor)
    def readPointsFile(self):
        pointsFile = open(self.name + ".points", "r")
        self.numPoints = int(pointsFile.readline())
        pointsList = {}
        for row in pointsFile:
            ll = row.split()
            pointsList[ll[2]] = [float(ll[0]), float(ll[1]), int(ll[3])]
        pointsFile.close()
        self.pointsList = pointsList
    def readCellsTabFile(self):
        cells = pd.read_csv(self.name + ".celltab", sep="\t")
        self.celltab = cells
        if(not self.celltypes):
            self.celltypes = list(cells.type) 
        if(not self.polygonList):
            self.numCells = cells.shape[0]
            vertices = list(map(lambda x: [i for i in x.split(',') if i != ''], cells.vertices))
            self.polygonList = [[(self.pointsList[coord][0], self.pointsList[coord][1]) for coord in cell] for cell in vertices]       
    def readCellsFile(self):
        cellsFile = open(self.name + ".cells", "r")
        self.numCells, numVerticesCell = [int(each) for each in cellsFile.readline().split()]
        polygonList = []
        cellvertList = []
        celltypes = []
        for row in cellsFile:
            polygonIndex = [indPoint for indPoint in row.split() if int(indPoint) >= 0]
            celltypes.append(int(polygonIndex.pop())) #CAREFUL: if cells file does not have type, will produce error
            polygonCoords = [[self.pointsList[coord][0], self.pointsList[coord][1]] for coord in polygonIndex]
            polygonList.append(polygonCoords) #Useful for plots
            cellvertList.append(polygonIndex) #Useful for graphs etc
        cellsFile.close()
        self.polygonList, self.cellvertList, self.celltypes = (polygonList, cellvertList, celltypes)
    def readEdgeTabFile(self):
        self.edgetab = pd.read_csv(self.name + ".edges", sep="\t")
    def readSprings(self):
        self.sprList = [] #just in case
        try:
            springsfile = open(self.name + ".spr", "r")
            self.numsprings = int(springsfile.readline())
            self.sprList = [[str(int(j)) for j in line.split('\t')] for line in springsfile ]
            springsfile.close()
        except:
            print("no springs (.spr) file")
    def setLimits(self, limits=[]):
        if(len(limits) == 0):
            self.limits = self.getLimits()
        else: 
            self.limits = limits
    def getLimits(self):
        aux = next(iter(self.pointsList.values()))
        xmin = aux[0]
        ymin = aux[1]
        xmax = aux[0]
        ymax = aux[1]
        for x, y, t in self.pointsList.values():
            if (x < xmin): xmin = x
            if (x > xmax): xmax = x
            if (y < ymin): ymin = y
            if (y > ymax): ymax = y
        rx = abs(xmax - xmin)
        ry = abs(ymax - ymin)
        xmin -= 0.1*rx
        xmax += 0.1*rx
        ymin -= 0.1*ry
        ymax += 0.1*ry
        #print(xmin, ymin, xmax, ymax)
        return (xmin, ymin, xmax, ymax)
    def getBorder(self):
        eborder = self.edgetab[self.edgetab['type'] == BORDER_TYPE]
        edgev = [[int(i) for i in x.split(',') if i] for x in eborder['vertices']]
        edgel = eborder.ind.tolist()
        border = []
        border_e = []
        start = edgev.pop()
        border.append(start[0])
        border.append(start[1])   
        border_e.append(edgel.pop())
        current = start[1]
        while(edgev):
            for i in range(len(edgev)):
                if(current == edgev[i][0]): 
                    e = edgev.pop(i)
                    el = edgel.pop(i)
                    border.append(e[1]) 
                    border_e.append(el)
                    current = e[1]                 
                    break
                elif(current == edgev[i][1]): 
                    e = edgev.pop(i)
                    el = edgel.pop(i)
                    border.append(e[0]) 
                    border_e.append(el)
                    current = e[0]                 
                    break
            else:
                print("Error: border not continuous")
                return []
        return (border, border_e)
    def getMarkPoints(self, points_x, points_y):
        xmin, ymin, xmax, ymax = self.getLimits()
        coord_x = [] 
        for p in points_x:
            coord_x.append(xmin + (xmax - xmin)*p)
        coord_y = []
        for p in points_y:
            coord_y.append(ymin + (ymax - ymin)*p)  
        return (coord_x, coord_y)
    def getClosestCell(self, coord, points, edges, ind=0):
        edgepos = [0.5*(self.pointsList[str(points[i-1])][ind] + self.pointsList[str(points[i])][ind]) for i in range(1, len(points))]
        edgedist = [abs(x - coord) for x in edgepos]
        sorted_edges = [edges[i] for i in np.argsort(edgedist)]
        getcell = lambda x: [d for d in [c for c in self.edgetab.loc[self.edgetab.ind == x].cells.values[0].split(',') if c] if int(d) >= 0][0]
        return list(map(getcell, sorted_edges))
        #closest_edge = edges[np.argmin(edgedist)] 
        #edge_cells = [c for c in self.edgetab.loc[self.edgetab.ind == closest_edge].cells.values[0].split(',') if c]
        #edge_cells = [int(c) for c in edge_cells if int(c) >= 0]
        #return edge_cells[0]
    def getPathExtremesCells(self, points_y=[0.2, 0.05,0.4,0.4,0.6,0.6,0.7,0.8], points_x=[0.1, 0.6, 0.6], thickness = 2):
        px, py = self.getMarkPoints(points_x, points_y)
        py = [(py[i], py[i + 1]) for i in range(0, len(py), 2)]
        border, border_e = self.getBorder()
        proximal_border = [i for i in range(len(border)) if self.pointsList[str(border[i])][0] < px[0]] #px[0] contains x coord at the limit of the "proximal part" (which is arbitrary)
        distal_border = [i for i in range(len(border)) if self.pointsList[str(border[i])][0] > px[1]] #px[1] contains x coord at the limit of the "distal part" (which is arbitrary)
        proximal_points = [border[i] for i in proximal_border]
        proximal_border.pop()
        proximal_edges = [border_e[i] for i in proximal_border]
        distal_points = [border[i] for i in distal_border]
        distal_border.pop()
        distal_edges = [border_e[i] for i in distal_border]     
        cell_pairs = []   
        for vein1, vein2 in py:
            cell_prox = self.getClosestCell(vein1, proximal_points, proximal_edges, 1)
            cell_dist = self.getClosestCell(vein2, distal_points, distal_edges, 1)   
            for layer in range(thickness):       
                cell_pairs.append((cell_prox[layer], cell_dist[layer]))
        import networkx as nx
        elist = [(i, j) for i in range(self.numCells - 1) for j in range(i, self.numCells) if len(set(self.cellvertList[i]) & set(self.cellvertList[j])) == 2]
        G = nx.Graph()
        G.add_edges_from(elist)
        paths = [nx.dijkstra_path(G, int(a), int(b)) for a, b in cell_pairs]
        self.setAsVeins(paths)
        """         
        pairs_merged = [i for p in cell_pairs for i in p] #to avoid removing these from Graph
        for layer in range(thickness):
            paths = [nx.dijkstra_path(G, a, b) for a, b in cell_pairs]
            self.setAsVeins(paths)
            paths_merged = [j for i in paths for j in i if j not in pairs_merged]
            G.remove_nodes_from(paths_merged) """
    def setAsVeins(self, paths):
        for p in paths: 
            for c in p:
                if(self.celltypes[c] < 2):
                    self.celltypes[c] += 2
    def getCellGraph(self):
        import networkx as nx
        elist = [(i, j) for i in range(self.numCells - 1) for j in range(i, self.numCells) if len(set(self.cellvertList[i]) & set(self.cellvertList[j])) == 2]
        G = nx.Graph()
        G.add_edges_from(elist)
        return G


def types_of_a_to_b(aname="", bname="", dead_in_a=[2131, 2132]):
    a = wing(aname)
    b = wing(bname)
    dead = 0
    for i in range(len(a.celltypes)):
        print(i, dead)
        if(i in dead_in_a):
            dead+=1     
            print("**", i)   
        b.celltypes[i + dead] = a.celltypes[i]
    return (a, b)

def rewrite_cells_file(b):
    pols = []
    with open(b.name + '.cells', "r") as f:
        numcells, vpc = [int(i) for i in f.readline().split("\t")]
        for i in range(numcells):
            line = f.readline().split('\t')
            line.pop()
            line.append(str(b.celltypes[i]))
            pols.append(line)
    s = ['\t'.join([str(numcells), str(vpc)])]
    s.extend(['\t'.join([str(j) for j in i]) for i in pols])
    s = '\n'.join(s)
    f = open(b.name + 'ct.cells', "w")
    f.write(s)
    f.close()     

""" def setGraphPaths(w, cell1, cell2):
    import networkx as nx
    elist = [(i, j) for i in range(w.numCells - 1) for j in range(i, w.numCells) if len(set(w.cellvertList[i]) & set(w.cellvertList[j])) == 2]
    G = nx.Graph()
    G.add_edges_from(elist)
    path = nx.dijkstra_path(G,cell1,cell2) """

def make_veins_thinner(aname, iters=1):
    a = wing(aname)
    G = a.getCellGraph()
    for i in range(iters):
        newtypes = dict()
        for c in range(a.numCells):
            if(a.celltypes[c] != veintype and a.celltypes[c] != veinhinge):
                continue
            for nei in G.neighbors(c):
                if(a.celltypes[nei] != veintype and a.celltypes[nei] != veinhinge):
                    newtypes.setdefault(c, a.celltypes[nei])
                    break
        for c, t in newtypes.items():
            a.celltypes[c] = t
    rewrite_cells_file(a)
    return a

parser = argparse.ArgumentParser(description='Assigns cell types from wing A to wing B. This can be useful if we have simulated, for instance, expansion (6-8h APF) and then we have defined cell types on the output wing using the initial_conds/draw_initial_wing.py script. If we want to perform the simulation from 6 to 8h and from 16 to 32h all at once, we need to have the cell types assigned from the beginning. This script accomplishes this. Also, it can "fix" vein issues in output files: fill gaps if veins are broken, make veins thinner if they were originally too thick. Note: this script is very slow for big wings and does not always work well.')
parser.add_argument('-i', '--Inputname', metavar='inputname', type=str, default = "", 
                    help='Identifier of wing to set cell types.')
parser.add_argument('-i2', '--Input2', metavar='input2', type=str, default = "", 
                    help='Identifier of wing from which to copy cell types (only applies if mode="copy"). Also, if mode is not "copy", use this for input.')
parser.add_argument('-d', '--DeadInInput', metavar='start', type=str, default = "", 
                    help='Indices of cells that are present in PUT_CELLTYPES but are dead in CELLTYPES_FROM (only applies if mode="copy"). Format X,X,X...')
parser.add_argument('-m', '--Mode', metavar='mode', type=str, default = 'copy', 
                    help='Mode: One of copy, rect (not implemented), graph, make_thin. \ncopy: copies cell types of wing in argument -i2 to wing in argument -i\ngraph: fills gaps in veins (of -i2) using a graph shortest path algorithm to connect isolated pieces of veins. Works for small gaps.\nmake_thin: makes veins thinner by turning vein cells in contact with hinge/blade into hinge/blade.')

def main():
    #a is the one with cell types. It is a later developmental stage
    #b is the one without cell types. It is an earlier developmental stage (just a bud)
    #Divisions are not important since new cells will be at the end
    #However, T2 must be taken into account. Cells that have died in the process from b to a must be known and passed as argument as "X,X,X"
    #Type must be assigned manually later to new cells and cells that died.
    #Therefore it is recommended to use this script in an interactive way, instead of just executing it
    args = parser.parse_args()
    aname = args.Input2
    bname = args.Inputname
    mode = args.Mode
    if(mode == "copy"):
        dead_in_a = [int(i) for i in sys.argv[3].split(',')]  
        a, b = types_of_a_to_b(aname, bname, dead_in_a)
        a.plot()
        plt.show()
        b.plot()
        plt.show()
        rewrite_cells_file(b)
    elif(mode == "graph"):
        a = wing(aname)
        a.getPathExtremesCells()
        a.plot()
        plt.show()
        b = a
        rewrite_cells_file(b)
    elif(mode == "rect"):
        pass
    elif(mode == "make_thin"):
        b = make_veins_thinner(aname, iters=1)
        b.plot()
        plt.show()
        rewrite_cells_file(b)
    # Example of what should be done manually after plotting wing b:
    # b.celltypes[dead_in_a[0] ] = 1

if __name__ == "__main__":
    main()

