import numpy as np
import sys
import argparse
import os
from shapely.ops import polygonize
from shapely.geometry import Polygon, MultiPoint, Point
from descartes.patch import PolygonPatch
from math import sqrt
import matplotlib.pyplot as plt

########################################################################################################################
# Plotting operations in general
########################################################################################################################

# Replacing the fiugures module
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


GM = (sqrt(5)-1.0)/2.0
W = 8.0
H = W*GM
SIZE = (W, H)

COLOR_ISVALID = {
    True: BLUE,
    False: RED,
}

def plot_coords(ax, ob, color=GRAY, zorder=1, alpha=1):
    x, y = ob.xy
    ax.plot(x, y, '.', color=color, zorder=zorder, alpha=alpha)

def color_isvalid(ob, valid=BLUE, invalid=RED):
    if ob.is_valid:
        return valid
    else:
        return invalid

def set_limits(ax, x0, xN, y0, yN):
    ax.set_xlim(x0, xN)
    ax.set_xticks(range(x0, xN+1))
    ax.set_ylim(y0, yN)
    ax.set_yticks(range(y0, yN+1))
    ax.set_aspect("equal")

def plot_grid(plot_pos, grid, pointsList, sprList, add_vnums, celltypes, expr, name):
    fig, ax = plt.subplots()
    # In order to plot a MultiPolygon object, I need to iterate over each oplygon
    fig = plt.figure(1, figsize=(5,5), dpi=90)
    ax = fig.add_subplot(plot_pos)
    x, y, t = zip(*pointsList.values())
    #ww = max(x)*1.05
    #hh = max(y)*1.05
    #aa = max([ww, hh])
    ax.set_aspect('equal')
    count = 0
    #ax.scatter(aa, aa, alpha=0.1)
    for k, element in enumerate(grid):
        count += 1
        polygon = Polygon(element)
        plot_coords(ax, polygon.exterior)

        if(len(expr) == 0):
            if(len(celltypes) == 0):
                col = BLUE
            else:
                col = wingcols[celltypes[k]]
            col = color_isvalid(polygon, col)
            transp = 0.5
        else:
            col = wingcols[celltypes[k]]
            transp = expr[k]
       # if(len(celltypes) == 0):
       #     edgec = "#000000"
       # else:
       #     edgec = color_isvalid(polygon, valid=wingcols[celltypes[k]])
        edgec = color_isvalid(polygon, valid=BLACK)
        patch = PolygonPatch(polygon, facecolor=col, edgecolor=edgec, alpha=transp, zorder=2)
        ax.add_patch(patch)
    for s in sprList:
        ax.plot([pointsList[s[0]][0], pointsList[s[1]][0]], [pointsList[s[0]][1], pointsList[s[1]][1]], color = RED)
    #for p in pointsList.keys():
    #    if(pointsList[p][2] == 0):
    #        ax.scatter(pointsList[p][0], pointsList[p][1], color = RED)
       
    if(add_vnums):
        for i in pointsList.keys():
            ax.annotate(i, pointsList[i][0:2], fontsize = 'xx-small')
    plt.savefig(name + '.svg', format='svg', dpi=1200)
    plt.savefig(name + '.png', format='png')
    plt.clf()

##Not tested, probably doesn't works
def plot_grid2(plot_pos, grid, pointsList, sprList, add_vnums, celltypes, expr, name, limits, figg=None):
    from matplotlib.collections import PolyCollection

    if(figg is None):
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
    else:
        fig, ax = figg
        ax.collections.clear()

    plt.xlim(limits[0], limits[2])
    plt.ylim(limits[1], limits[3])
    #col = np.zeros((len(grid), 6, 2))
    pc = PolyCollection(grid, facecolors= [wingcols[celltypes[j]] for j in range(len(grid))], alpha = 0.6)
    pc.set_edgecolors("black")
    ax.add_collection(pc)

    for s in sprList:
        spcolor = springcols[int(s[2])] if len(s) > 2 else RED
        ax.plot([pointsList[s[0]][0], pointsList[s[1]][0]], [pointsList[s[0]][1], pointsList[s[1]][1]], color = spcolor)
    #for p in pointsList.keys():
    #    if(pointsList[p][2] == 0):
    #        ax.scatter(pointsList[p][0], pointsList[p][1], color = RED)       
    if(add_vnums):
        for i in pointsList.keys():
            ax.annotate(i, pointsList[i][0:2])

    plt.savefig(name + '.png')
    plt.savefig(name + '.svg', format='svg', dpi=1200)
    return (fig, ax)
    


def plot_expr(plot_pos, grid, pointsList, sprList, add_vnums, celltypes, expr, name, color_expr):
    for g in color_expr:    
        xx = expr[:,g].tolist()
        #xmin = min(xx)
        xmax = max(xx)
        xx = [i/xmax for i in xx]
        #print("_____ ", xx)
        plot_grid(plot_pos, grid, pointsList, sprList, add_vnums, celltypes, xx, name + 'g' + str(g))


def readPointsFile(name):
    pointsFile = open(name + ".points", "r")
    numPoints = int(pointsFile.readline())

    # This is the list of points
    pointsList = {}#np.zeros((numPoints,2))
    for row in pointsFile:
        ll = row.split()
        pointsList[ll[2]] = [float(ll[0]), float(ll[1]), int(ll[3])]
    pointsFile.close()
    return (numPoints, pointsList)

def readCellsFile(name, plot_cell_types, pointsList):
    cellsFile = open(name + ".cells", "r")
    numCells, numVerticesCell = [int(each) for each in cellsFile.readline().split()]

    polygonList = []
    celltypes = []
    for row in cellsFile:
        polygonIndex = [indPoint for indPoint in row.split() if int(indPoint) >= 0]
        if(plot_cell_types):
            celltypes.append(int(polygonIndex.pop())) #CAREFUL: if cells file does not have type, will produce error
        polygonCoords = [[pointsList[coord][0], pointsList[coord][1]] for coord in polygonIndex]
        polygonList.append(polygonCoords)
    cellsFile.close()
    return (numCells, polygonList, celltypes)

def readSprings(name):
    sprList = [] #just in case
    numsprings = 0
    try:
        springsfile = open(name + ".spr", "r")
        numsprings = int(springsfile.readline())
        #sprList = [[str(int(j)) for j in springsfile.readline().split('\t')] for i in range(numsprings) ]
        sprList = [[str(int(j)) for j in line.split('\t')] for line in springsfile ]
        springsfile.close()
    except:
        print("no springs (.spr) file")
    return (numsprings, sprList)


def readExpr(name):
    xprList = np.loadtxt(name + '.expr') #just in case
    return xprList
    
def getLimits(pointsList):
    aux = next(iter(pointsList.values()))
    xmin = aux[0]
    ymin = aux[1]
    xmax = aux[0]
    ymax = aux[1]
    for x, y, t in pointsList.values():
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
    print(xmin, ymin, xmax, ymax)
    return (xmin, ymin, xmax, ymax)
########################################################################################################################
# Generating the multipolygon object
########################################################################################################################

parser = argparse.ArgumentParser(description='Plot grid arguments.')
parser.add_argument('-i', '--Inputname', metavar='inputname', type=str, default = "hexgrid", 
                    help='Identifier. Used as prefix of input files of cell coordinates. ')
parser.add_argument('-i2', '--Input_expr', metavar='input_expr', type=str, default = "", 
                    help='Identifier. Used as prefix of input files of gene expression. ')
parser.add_argument('-s', '--Start_index', metavar='start', type=int, default = 0, 
                    help='First file to read')
parser.add_argument('-e', '--End_index', metavar='end', type=int, default = 1000, 
                    help='Last file to read')
parser.add_argument('-v', '--write_vertex_number', metavar='writevert', type=bool, default = False, 
                    help='Add number of vertices to plot (recommended only for small grids)')
parser.add_argument('-g', '--genes_to_plot_expression', metavar='gexpr', type=str, default = '', 
                    help='Plot grid with expression of all these genes. Integers separated by commas.')
parser.add_argument('-t', '--plotCellTypes', metavar='plot_cell_types', type=bool, default = True, 
                    help='Color according to cell type or not')


def main():
    args = parser.parse_args()
    plot_cell_types = args.plotCellTypes
    add_vnums = args.write_vertex_number
    color_expr = [int(i) for i in args.genes_to_plot_expression.split(',') if i != '']
    fig = None
    for fnum in range(args.Start_index, args.End_index):
        if(fnum >= 0):
            name = args.Inputname + str(fnum)
        else:
            name = args.Inputname
        if(not os.path.isfile(name + ".points") or not os.path.isfile(name + ".cells")):
            break;
        numPoints, pointsList = readPointsFile(name)
        numCells, polygonList, celltypes = readCellsFile(name, plot_cell_types, pointsList)
        numsprings, sprList = readSprings(name)
        ########################################################################################################################
        # Plotting the final polygons
        ########################################################################################################################
        print("%d files read..."%fnum)
        if(fnum == 0):
            limits = getLimits(pointsList)
        fig = plot_grid2(111, polygonList, pointsList, sprList, add_vnums, celltypes, [] ,name, limits)  
        if(len(color_expr) > 0 and args.Input_expr != ""):
            xprList = readExpr(args.Input_expr + str(fnum))
            plot_expr(111, polygonList, pointsList, sprList, add_vnums, celltypes, xprList, name, color_expr)
        print("%d printed"%fnum)


if __name__ == '__main__':
    main()
