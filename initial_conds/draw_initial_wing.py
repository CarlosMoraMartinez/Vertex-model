import matplotlib.pyplot as plt 
import numpy as np
import scipy.spatial as spatial
from scipy.spatial import Voronoi, voronoi_plot_2d

import shapely
from shapely.geometry import Polygon

import make_hexagonal_grid as mhg

def makeWingHex():
    argdict = {}
    argdict.setdefault('Size', 3)
    argdict.setdefault('Rows', 100)
    argdict.setdefault('Cols', 200)
    argdict.setdefault('Noise', 0.3)
    argdict.setdefault('Pull', 0)
    argdict.setdefault('Strecht', 1.1)
    argdict.setdefault('Rotate', 0)
    argdict.setdefault('StaticVertices', '')
    argdict.setdefault('Springs', '')
    argdict.setdefault('SpringLength', 0)
    argdict.setdefault('Hinge', -1)
    argdict.setdefault('Veins', '')
    argdict.setdefault('Outname', 'hex')
    # Make hexagonal grid
    hx = mhg.HexGrid(**argdict)
    fig, ax = hx.plotHex2(False)
    #fig, ax = plt.subplots()
    im = plt.imread("wing.png")
    maxx, maxy = hx.getMax()
    ax.imshow(im, extent = [0, maxx, 0, maxy])
    # Draw border by hand 
    num_points, points, lines, fig, ax = getPointsBorder(fig, ax)
    plt.close()
    hx = copyModifiedStructure(points, lines, hx)
    f, ax = hx.plotHex2()
    ax.imshow(im, extent = [0, maxx, 0, maxy])
    plotBorder(points, lines, (f, ax))
    plt.show()
    ax.imshow(im, extent = [0, 100, 0, 100])


def makeWingVoronoi():
    npoints = 2000
    fig, ax = plt.subplots()
    im = plt.imread("wing.png")    
    ax.imshow(im, extent = [0, 100, 0, 100])
    num_points, points, lines, fig, ax = getPointsBorder(fig, ax)
    lim = Polygon(points)
    ranpoints = generatePointsInsideShape(npoints, lim)
    vor = Voronoi(ranpoints)
    f, ax = plt.subplots()
    voronoi_plot_2d(vor, ax=ax)
    plotBorder(points, lines, (f, ax))
    plt.show()

    vor = improveVoronoi(vor, lim)

    #print(vor.regions)
    #plotVorManual(vor, points, lines)
    f, ax = plt.subplots()
    voronoi_plot_2d(vor=vor, ax=ax)
    plotBorder(points, lines, (f, ax))
    plt.show()

def improveVoronoi(vor, lim, iters=10):
    for i in range(iters):
        cents = getCentroids(vor)
        vor = Voronoi(cents)
    return vor

def getCentroids(vor):
    centroids = []
    for r in vor.regions:
        if(-1 in r or len(r) < 3): 
            continue
        pol = []
        #print(r)
        for p in r:
            pol.append((vor.vertices[p, 0], vor.vertices[p, 1]))
        shpol = Polygon(pol)
        #print (shpol)
        centroids.append((shpol.centroid.x, shpol.centroid.y))
    return centroids

def filterVoronoi(vor, lim):
    filtered_regions = []
    for r in vor.regions:
        flag = True
        for p in r:
            if(-1 in vor.vertices[p]):
                flag = False
                break
            x = vor.vertices[p, 0]
            y = vor.vertices[p, 1]
            shpoint = shapely.geometry.Point(x, y)
            if(not lim.contains(shpoint)):
                flag = False
                break
        if(flag):
            filtered_regions.append(r)
    vor.regions = filtered_regions
    return vor
            
def plotVorManual(vor, points, lines):
    f, ax = plt.subplots()
    plotBorder(points, lines, (f, ax))
    for r in vor.regions:
        for p in range(len(r) - 1):
            ax.plot(vor.vertices[r[p]], vor.vertices[r[p+1]] , c = "red")
        #ax.plot(vor.vertices[r[len(r)-1]], vor.vertices[r[0]] , c = "red")
    plt.show()
              


def generatePointsInsideShape(npoints, lim):
    xmin, xmax, ymin, ymax = (min(lim.exterior.xy[0]), max(lim.exterior.xy[0]), min(lim.exterior.xy[1]), max(lim.exterior.xy[1]))
    ranpoints = []
    while(len(ranpoints) < npoints):
        ranx = np.random.uniform(low = xmin, high = xmax)
        rany = np.random.uniform(low = ymin, high = ymax)  
        shpoint = shapely.geometry.Point(ranx, rany)
        if(lim.contains(shpoint)):
            ranpoints.append((ranx, rany))
    return ranpoints
  

def plotBorder(points, lines, fig):
    f, ax = fig
    ax.plot([points[len(points) - 1][0], points[0][0]], [points[len(points) - 1][1], points[0][1]], color = "black")
    for i in range(len(points) - 1):
         ax.plot([points[i][0], points[i+1][0]], [points[i][1], points[i+1][1]], color = "black")

def removeCellsOutside(points, lines, hx):
    hxpol = hx.getShapelyPolygons()
    lim = Polygon(points)
    cells_to_remove = []
    cells_to_fit = []
    for i, cell in enumerate(hxpol):
        if(not lim.contains(cell)):
            if(lim.overlaps(cell)):
                if(lim.intersection(cell).area/cell.area < 0.5):
                    cells_to_remove.append(i)
                else:
                    cells_to_fit.append(i)
            else:
                cells_to_remove.append(i)
    cells_to_remove.sort(reverse = True)
    for i in cells_to_remove:
        hx.removeCell(i)
    return hx

def copyModifiedStructure(points, lines, hx):
    hxpol = hx.getShapelyPolygons()
    lim = Polygon(points)
    cells_to_use = []
    cells_to_fit = []
    for i, cell in enumerate(hxpol):
        if(lim.contains(cell)):
            cells_to_use.append(i)
        elif(lim.overlaps(cell)):
            if(lim.intersection(cell).area/cell.area > 0.2):
                    cells_to_fit.append(i)
    cells_to_use.sort()
    cells_to_fit.sort()

    cells = []
    vertices = []
    vertDict = dict()
    centers= []
    celltypes= []
    springs= []

    vnum = 0
    for cind in cells_to_use:
        c = hx.cells[cind] 
        for i, v in enumerate(c):
            if(v in vertDict):
                c[i] = vertDict[v]
            else:
                vertDict.setdefault(v, vnum)
                c[i] = vnum
                vertices.append( hx.vertices[v] )
                vertices[vnum][2] = vnum
                vnum+=1                
        cells.append(c)
        celltypes.append(hx.celltypes[cind])
        centers.append(hx.centers[cind])
    for s in hx.springs:
        isthere = 0 if s[0] in vertDict else 1 if s[1] in vertDict else -1
        if(isthere == 0):
            newspring = (vertDict[s[0]], vnum)
            vertDict.setdefault(s[1], vnum)
            vertices.append( hx.vertices[s[1]] )
            vertices[vnum][2] = vnum
            vnum+=1                
        elif(isthere == 1):
            newspring = (vnum, vertDict[s[1]])
            vertDict.setdefault(s[0], vnum)
            vertices.append( hx.vertices[s[0]] )
            vertices[vnum][2] = vnum
            vnum+=1                   
    hx.cells = cells
    hx.vertices = vertices
    hx.springs = springs
    hx.celltypes = celltypes
    hx.centers = centers
    return hx       

def getPointsBorder(fig, ax):
    points = []
    lines = []
    num_points = 0
    finished = False
    def onclick(event):
        nonlocal num_points
        nonlocal finished
        if(not finished):
            print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
                   ('double' if event.dblclick else 'single', event.button,
                   event.x, event.y, event.xdata, event.ydata))
            points.append([event.xdata, event.ydata])
            ax.scatter([event.xdata], [event.ydata], color = "blue")
            if(num_points > 0):
                lines.append([num_points - 1, num_points])
                ax.plot([points[num_points - 1][0], points[num_points][0]], [points[num_points - 1][1], points[num_points][1]], color = "black")
            fig.canvas.draw()
            num_points += 1
    def on_key(event):
        nonlocal finished
        nonlocal num_points
        print('you pressed', event.key, event.xdata, event.ydata)
        if(event.key in "mM" and not finished):
            finished = True
            nonlocal cid2
            fig.canvas.mpl_disconnect(cid2)
            lines.append([num_points - 1, 0])
            ax.plot([points[num_points - 1][0], points[0][0]], [points[num_points - 1][1], points[0][1]], color = "black")
            fig.canvas.draw()
        elif(finished):
            nonlocal cid1
            fig.canvas.mpl_disconnect(cid1)
            plt.close()
            return
    cid1 = fig.canvas.mpl_connect('key_press_event', on_key)
    cid2 = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    return (num_points, points, lines, fig, ax)

#makeWingVoronoi()
makeWingHex()
