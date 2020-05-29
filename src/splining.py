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

outfile = "landmarks.csv"


class spline:
    landmark_names=[j+i for i in ['1', '2', '3', '4', '12'] for j in ['x', 'y'] ]
    def __init__(self, name):
        self.name = name
        self.wing = wing(self.name)
        self.points=[]
        self.fig, self.ax = plt.subplots()
        self.ax.set_aspect('equal')
        self.ax.set_title(self.name)
        self.wing.setAx(self.fig, self.ax)
        self.wing.plot()
        self.poly = PolygonSelector(self.ax, self.onselect)
    def onselect(self, points):
        self.points = points
    def spline(self):
        mng = plt.get_current_fig_manager() #These two lines are to maximize window
        mng.window.showMaximized()
        plt.show()
        plt.close()
        return dict(zip(["name"] + self.landmark_names, [self.name] + [j for i in self.points for j in i]))



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
        self.fig = None
        self.ax = None
        try:
            self.readData()
        except:
            print("Unable to read data for wing %s"%(self.name))
    def readData(self):
        self.readPointsFile()
        self.readCellsFile()
        self.readSprings()
    def setAx(self, fig, ax, lim=[]):  
        self.fig=fig
        self.ax=ax
        self.setLimits(lim)
    def plot(self):
        from matplotlib.collections import PolyCollection  
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
    def readCellsFile(self):
        cells = pd.read_csv(self.name + ".celltab", sep="\t")
        self.numCells = cells.shape[0]
        self.celltypes = list(cells.type) 
        vertices = list(map(lambda x: [i for i in x.split(',') if i != ''], cells.vertices))
        self.polygonList = [[(self.pointsList[coord][0], self.pointsList[coord][1]) for coord in cell] for cell in vertices]       
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

def getWingSet(names):
    wings = {}
    for n in names:
        wings.setdefault(n, wing(n))
    return wings


def getLimits(wings):
    limits = []
    for w in wings.values():
        limits.append(w.getLimits())
    x0, y0, x1, y1 = [i for i in zip(*limits)]
    return (min(x0), min(y0), max(x1), max(y1))
    

def plotAll(df):
    print("PLOTTING ALL WINGS")
    wings = getWingSet(df.name)
    limits = getLimits(wings)
    df = df.drop("name", axis = 1)
    #maxims = df.max()
    #minims = df.min()
    #xmax = np.max( maxims.iloc[list(map(lambda x: "x" in x, maxims.index))] )
    #ymax = np.max( maxims.iloc[list(map(lambda x: "y" in x, maxims.index))] )
    #xmin = np.min( minims.iloc[list(map(lambda x: "x" in x, minims.index))] )
    #ymin = np.min( minims.iloc[list(map(lambda x: "y" in x, minims.index))] )
    print("COORD LIMITS: ", limits)
    #print(df)
    for n, r in df.iterrows():
        w = wings[n]
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        w.setAx(fig, ax, limits)
        w.plot()
        r = np.array(r).reshape([5, 2])
        plt.scatter(r[:,0], r[:,1], color = "red")
        r2 = np.zeros([6, 2])
        r2[:5,:] = r
        r2[5,:] = r[0,:]
        ax.plot(r2[:,0], r2[:,1], color = "red")
        ax.set_title(n)
        plt.savefig(n + '_landmarks.png', dpi=400)
        plt.close()
        print("plotted ", n)

def main():
    files = os.listdir()
    if(outfile in files):
        try:
            df = pd.read_csv(outfile, index_col=0)
            print("WARNING: Some landmarks already measured: ")
            for l in df.name:
                print("\t", l)
        except:
            print("Detected landmark file but unable to read it. Splining all wings")
    else:
        df = pd.DataFrame(columns = ["name"] + spline.landmark_names)
        df.set_index("name", drop=False)
    wings = [i.replace(".points", "") for i in files if i.endswith(".points")]
    wings = [w for w in wings if not w in list(df.name)]
    #landmarks = []
    for w in wings:
        try:
            spl = spline(w)
            lms = spl.spline()
            #print(" ".join([str(i) for i in lms.values()]))
            lms = pd.Series(lms, name=lms["name"])
            print(lms)
            #landmarks.append(lms)
            df = df.append(lms) #write every time so you don't lose any
            df.to_csv(outfile)
            #print(df)
        except:
            print("ERROR: wing %s not measured."%(w))
    plotAll(df)

if __name__ == '__main__':
    main()



