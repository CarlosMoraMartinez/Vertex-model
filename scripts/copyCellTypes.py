from operator import itemgetter
import sys
import os
import argparse

import numpy as np
import scipy.spatial as spatial
from scipy.spatial import Voronoi, voronoi_plot_2d
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
        self.celltypes = []
        self.numsprings = 0
        self.sprList = []
        self.celltab = None
        self.fig = None
        self.ax = None
        try:
            self.readData()
        except:
            print("Unable to read data for wing %s"%(self.name))
    def readData(self):
        self.readPointsFile()
        self.readCellsFile()
        try:
            self.readCellsTabFile()
        except:
            print("No .celltab file")
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
        celltypes = []
        for row in cellsFile:
            polygonIndex = [indPoint for indPoint in row.split() if int(indPoint) >= 0]
            celltypes.append(int(polygonIndex.pop())) #CAREFUL: if cells file does not have type, will produce error
            polygonCoords = [[self.pointsList[coord][0], self.pointsList[coord][1]] for coord in polygonIndex]
            polygonList.append(polygonCoords)
        cellsFile.close()
        self.polygonList, self.celltypes = (polygonList, celltypes)
    def readSprings(self):
        sprList = [] #just in case
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


def main():
    #a is the one with cell types. It is a later developmental stage
    #b is the one without cell types. It is an earlier developmental stage (just a bud)
    #Divisions are not important since new cells will be at the end
    #However, T2 must be taken into account. Cells that have died in the process from b to a must be known and passed as argument as "X,X,X"
    #Type must be assigned manually later to new cells and cells that died.
    #Therefore it is recommended to use this script in an interactive way, instead of just executing it
    aname = sys.argv[1]
    bname = sys.argv[2]
    dead_in_a = [int(i) for i in sys.argv[3].split(',')]  
    a, b = types_of_a_to_b(aname, bname, dead_in_a)
    a.plot()
    plt.show()
    b.plot()
    plt.show()
    # Example of what should be done manually after plotting wing b:
    # b.celltypes[dead_in_a[0] ] = 1
    rewrite_cells_file(b)


