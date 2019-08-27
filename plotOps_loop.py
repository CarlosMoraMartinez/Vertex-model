import numpy as np
import sys
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


bladetype = 0
hingetype = 1
veintype = 2
veinhinge = 3
wingcols = ['blue', 'green', 'black', 'gray']


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

def plot_grid(plot_pos, grid, pointsList, sprList, add_vnums, celltypes):
    # In order to plot a MultiPolygon object, I need to iterate over each oplygon
    fig = plt.figure(1, figsize=(5,5), dpi=90)
    ax = fig.add_subplot(plot_pos)
    count = 0
    for k, element in enumerate(grid):
        count += 1
        polygon = Polygon(element)
        plot_coords(ax, polygon.exterior)
        if(len(celltypes) == 0):
            col = BLUE
        else:
            col = wingcols[celltypes[k]]
        col = color_isvalid(polygon, col)
        patch = PolygonPatch(polygon, facecolor=col, edgecolor=color_isvalid(polygon, valid=BLUE), alpha=0.5, zorder=2)
        ax.add_patch(patch)
    for s in sprList:
        ax.plot([pointsList[s[0]][0], pointsList[s[1]][0]], [pointsList[s[0]][1], pointsList[s[1]][1]], color = RED)
    for p in pointsList.keys():
        if(pointsList[p][2] == 0):
            ax.scatter(pointsList[p][0], pointsList[p][1], color = RED)
       
    if(add_vnums):
        for i in pointsList.keys():
            ax.annotate(i, pointsList[i][0:2])


########################################################################################################################
# Generating the multipolygon object
########################################################################################################################


if(len(sys.argv) <= 5):
    plot_cell_types = True
else:
    plot_cell_types = False

for fnum in range(int(sys.argv[2]), int(sys.argv[3])):
    pointsFile = open(sys.argv[1]+str(fnum)+".points", "r")
    cellsFile = open(sys.argv[1]+str(fnum)+".cells", "r")
    
    numPoints = int(pointsFile.readline())
    numCells, numVerticesCell = [int(each) for each in cellsFile.readline().split()]

# This is the list of points
    pointsList = {}#np.zeros((numPoints,2))
    for row in range(numPoints):
        ll = pointsFile.readline().split()
        pointsList[ll[2]] = [float(ll[0]), float(ll[1]), int(ll[3])]

#print(pointsList)
    polygonList = []
    celltypes = []
    for row in range(numCells):
        polygonIndex = [indPoint for indPoint in cellsFile.readline().split() if int(indPoint) >= 0]
        if(plot_cell_types):
            celltypes.append(int(polygonIndex.pop())) #CAREFUL: if cells file does not have type, will produce error
        polygonCoords = [list(pointsList[coord]) for coord in polygonIndex]
        polygonList.append(polygonCoords)
    #olygonList[row] =
#springlist
    sprList = [] #just in case
    try:
        springsfile = open(sys.argv[1]+str(fnum)+".spr", "r")
        numsprings = int(springsfile.readline())
        sprList = [[str(int(j)) for j in springsfile.readline().split('\t')] for i in range(numsprings) ]
    except:
        print("no springs (.spr) file")
        sprList = [] #just in case

########################################################################################################################
# Plotting the final polygons
########################################################################################################################
    add_vnums = False
    if(len(sys.argv) > 4):
        if(sys.argv[4] == '-n'):
            add_vnums = True
    fig, ax = plt.subplots()
    plot_grid(111, polygonList, pointsList, sprList, add_vnums, celltypes)

#plt.show()
    plt.savefig(sys.argv[1]+str(fnum)+'.png')
    plt.clf()
    print(fnum)


