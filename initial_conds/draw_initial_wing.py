from operator import itemgetter
import sys

import numpy as np
import scipy.spatial as spatial
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt 
from matplotlib.widgets import Slider, Button, Lasso, LassoSelector, TextBox
import shapely
from shapely.geometry import Polygon, Point

import make_hexagonal_grid as mhg

def makeWingHex(wname):
    # Make hexagonal grid
    hx, im = setHexagonSize("wing.png", wname).getGrid()
    fig, ax = plt.subplots()
    fig.suptitle(wname, fontsize=16)
    maxx, maxy = hx.getMax()
    ax.imshow(im, extent = [0, maxx, 0, maxy])

    # Draw border by hand 
    num_points, points, lines, fig, ax = getPointsBorder(fig, ax).getData()
    plt.close()
    hx = hx.copyModifiedStructure(points, lines)
    f, ax, pc = hx.plotHex2()
    ax.imshow(im, extent = [0, maxx, 0, maxy])
    plotBorder(points, lines, (f, ax))

    #Set cell types: hinge and veins
    getCellTypes(hx, pc, f, ax)
    f, ax, pc = hx.plotHex2()
    ax.imshow(im, extent = [0, maxx, 0, maxy])

    #Set borders: static vertices and springs
    getBorders(hx, pc, f, ax)
    hx.writeGrid()


class setHexagonSize:

    def __init__(self, imfile, hxname):

        self.im = plt.imread(imfile)
        self.imratio = self.im.shape[0]/self.im.shape[1]

        self.argdict = {}
        self.argdict.setdefault('Size', 3)
        self.argdict.setdefault('Rows', int(1.1*200*self.imratio))
        self.argdict.setdefault('Cols', 200)
        self.argdict.setdefault('Noise', 0.3)
        self.argdict.setdefault('Pull', 0)
        self.argdict.setdefault('Strecht', 1.1)
        self.argdict.setdefault('Rotate', 0)
        self.argdict.setdefault('StaticVertices', '')
        self.argdict.setdefault('Springs', '')
        self.argdict.setdefault('SpringLength', 0)
        self.argdict.setdefault('Hinge', -1)
        self.argdict.setdefault('Veins', '')
        self.argdict.setdefault('Outname', hxname)
        self.defaultArgs = self.argdict.copy()

        # Make hexagonal grid

        self.fig, self.ax = plt.subplots()
        plt.subplots_adjust(left=0.25, bottom=0.25)
        self.fig.suptitle(self.defaultArgs['Outname'], fontsize=16)
        self.canvas = self.fig.canvas
        self.cid1 = self.canvas.mpl_connect('key_press_event', self.on_key)

        self.makeAndAddToPlot()

        self.axcolor = 'lightgoldenrodyellow'
        self.axrows = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor=self.axcolor)
        self.axnoise = plt.axes([0.25, 0.10, 0.65, 0.03], facecolor=self.axcolor)
        self.axstrecht = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=self.axcolor)
        self.axrotate = plt.axes([0.25, 0.2, 0.65, 0.03], facecolor=self.axcolor)


        self.scols = Slider(self.axrows, 'Num. Cols', 5, 1000, valinit=self.argdict["Rows"], valstep=1)
        self.snoise = Slider(self.axnoise, 'Noise', 0, 0.5, valinit=self.argdict["Noise"], valstep=0.025)
        self.sstrecht = Slider(self.axstrecht, 'Strecht', 0.2, 5, valinit=self.argdict["Strecht"],valstep=0.05)
        self.srotate = Slider(self.axrotate, 'Rotate', -30, 30, valinit=self.argdict["Strecht"],valstep=1)

        self.scols.on_changed(self.update) 
        self.snoise.on_changed(self.update)
        self.sstrecht.on_changed(self.update) 
        self.srotate.on_changed(self.update)

        self.resetax = plt.axes([0.01, 0.025, 0.1, 0.04])
        self.button = Button(self.resetax, 'Reset', color=self.axcolor, hovercolor='0.975')
        self.button.on_clicked(self.reset)


        plt.show()

    def update(self, val):
        self.ax.clear()
        self.argdict['Rows'] = (self.scols.val*self.imratio*self.sstrecht.val).astype(int) 
        self.argdict['Cols']  =  self.scols.val.astype(int)
        self.argdict['Noise'] = self.snoise.val
        self.argdict['Strecht'] = self.sstrecht.val
        self.argdict['Rotate'] =  self.srotate.val
        self.makeAndAddToPlot() 


    def reset(self, event):
        self.argdict = self.defaultArgs.copy()
        self.scols.reset() 
        self.snoise.reset() 
        self.sstrecht.reset() 
        self.srotate.reset() 

    def makeAndAddToPlot(self):
        self.hx = mhg.HexGrid(**self.argdict)
        self.maxx, self.maxy = self.hx.getMax()
        self.ax.imshow(self.im, extent = [0, self.maxx, 0, self.maxy])
        self.hx.plotHex2((self.fig, self.ax), False, 0.3)
        self.canvas.draw_idle()
      
    def on_key(self, event):
        print('you pressed: ', event.key, event.xdata, event.ydata)
        if(event.key in "Mm"):
            self.fig.canvas.mpl_disconnect(self.cid1)
            plt.close()
            return

    def getGrid(self):
        return (self.hx, self.im)




class getBorders:
    def __init__(self, hx, pc, fig, ax):
        self.finished = False
        self.drawing = None
        self.staticPoints = []
        self.springPoints = []
        self.newSpringPoints = []
        self.springPlotRefs = dict()
        self.pointPlotRefs = dict()
        self.springType = 0
        self.hx = hx
        self.collection = pc
        self.fig = fig
        self.canvas = fig.canvas
        self.ax = ax
        self.cid1 = self.canvas.mpl_connect('key_press_event', self.on_key)
        self.cid2 = self.canvas.mpl_connect('button_press_event', self.onpress)

        plt.subplots_adjust(bottom=0.2)
        self.axbox = plt.axes([0.1, 0.05, 0.8, 0.075])

        print("T to set static vertices; Vv to set springs; Jj to remove static vertices or springs")
        plt.show()

    def submitText(self, text):
        self.springType = int(text)
        self.canvas.widgetlock.release(self.text_box)
        print("New springs will be of type: ", self.springType)
        del self.text_box
        self.axbox.clear()
        self.canvas.draw_idle()

    def on_key(self, event):

        print('you pressed: ', event.key, event.xdata, event.ydata)
        if(event.key in "tT" and not self.finished and not self.canvas.widgetlock.locked()):
            self.drawing = "static"
            print("Drawing Static vertices")
        elif(event.key in "Vv" and not self.finished and not self.canvas.widgetlock.locked()):
            self.drawing = "springs"
            print("Drawing Springs")
        elif(event.key in "Jj" and not self.finished and not self.canvas.widgetlock.locked()):
            self.drawing = "remove"
            print("Removing Static and Springs")
        elif (event.key in "Nn" and not self.finished and not self.canvas.widgetlock.locked()):
            self.text_box = TextBox(self.axbox, 'Spring type', initial=str(self.springType))
            self.text_box.on_submit(self.submitText)
            self.canvas.widgetlock(self.text_box)   
            self.canvas.draw_idle()   
            print("Getting new type of Springs")      
        elif(event.key in "Mm" and not self.finished and not self.canvas.widgetlock.locked()):
            self.finished = True
            self.drawing = None
        elif(self.finished):
            self.fig.canvas.mpl_disconnect(self.cid1)
            self.fig.canvas.mpl_disconnect(self.cid2)
            plt.close()
            return
    def onpress(self, event):
        print("Starting Lasso")

        if self.canvas.widgetlock.locked():
            return
        if event.inaxes is None:
            return
        if self.drawing is None:
            return
        # acquire a lock on the widget drawing
        self.lasso = LassoSelector(self.ax, onselect=self.callback)
        self.canvas.widgetlock(self.lasso)
        print("Finishing Lasso")

    def callback(self, verts):
        if(self.drawing == "static"):
            self.setStatic(verts)
        elif(self.drawing == "springs"):
            self.setSprings(verts)
        elif(self.drawing == "remove"):
            self.removeStatic(verts)
        #fig, ax, pc = self.hx.plotHex2(fig=(self.fig, self.ax))
        self.canvas.draw_idle()
        self.canvas.widgetlock.release(self.lasso)
        del self.lasso
        print("T to set static vertices; Vv to set springs; Jj to remove static vertices or springs; Nn to set spring type")

    def setStatic(self, verts):
        for v in verts:
            vstat = self.getClosestVert(v)
            if(vstat not in self.staticPoints):
                self.staticPoints.append(vstat)
        for sp in self.staticPoints:
            self.hx.vertices[sp][3] = 0
        self.updatePlot()
        #self.ax.scatter([self.hx.vertices[sp][0] for sp in self.staticPoints], [self.hx.vertices[sp][1] for sp in self.staticPoints], c="red")

    def setSprings(self, verts):
        vprevious = -1
        spring_positions = []
        springs = []
        for v in verts:
            vstat = self.getClosestVert(v)
            if(vstat == -1):
                continue
            if(vstat != vprevious and vstat not in self.springPoints):
                if(vprevious != -1):
                    x, y = zip(*spring_positions)
                    x = np.mean(x)
                    y = np.mean(y)
                    springs.append((vprevious, x, y, self.springType)) #Vertex in grid, and coordinates of new vertex
                spring_positions = [v] 
                vprevious = vstat
                self.springPoints.append(vstat)
            elif(vstat in self.springPoints):
                spring_positions.append(v)
        for s in springs:
            snum = self.hx.addSpring(*s)
            self.newSpringPoints.append(snum)
        self.updatePlot()

    def updatePlot(self):
        for sl in self.springPlotRefs.values():
            self.ax.lines.remove(sl)
        self.springPlotRefs.clear()
        for i, s in enumerate(self.hx.springs):
            spref = self.ax.plot([self.hx.vertices[s[0]][0], self.hx.vertices[s[1]][0]], [self.hx.vertices[s[0]][1], self.hx.vertices[s[1]][1]], c=mhg.springcols[s[2]])
            self.springPlotRefs.setdefault(i, spref[0])

        for pp in self.pointPlotRefs.values():
            pp.set_visible(False)
        self.pointPlotRefs.clear()
        for sp in self.staticPoints:
            ref = self.ax.scatter(self.hx.vertices[sp][0], self.hx.vertices[sp][1], c="red")
            self.pointPlotRefs.setdefault(sp, ref)
        self.canvas.draw_idle()

    def removeStatic(self, verts):
        springs_to_remove = []
        verts_to_remove = []
        verts_to_movable = []

        #Look for springs to remove and vertices to set as movable
        for v in verts:
            vstat = self.getClosestVert(v, 0)#get the non-movable colsest vertex
            if(vstat == -1):
                continue
            spr = self.hx.getSpringWithVertex(vstat)
            if(spr == -1):
                if(vstat not in verts_to_movable):
                    verts_to_movable.append(vstat)
                    #self.staticPoints.remove(vstat)
            elif(vstat not in verts_to_remove):
                springs_to_remove.append(spr)
                verts_to_remove.append(vstat)

        verts_to_remove.sort(reverse=True)
        springs_to_remove.sort(reverse=True)
        
        #Set vertices to movable
        for v in verts_to_movable:
            self.hx.vertices[v][3] = 1 

        #Remove springs and vertices that are only in springs
        self.hx.removeSpringVertices(springs_to_remove, verts_to_remove)

        ##Now update lists
        self.newSpringPoints.clear()
        self.springPoints.clear() 
        self.staticPoints.clear()

        for s in self.hx.springs:
            for v in s:
                if(self.hx.vertices[v][3] == 1):
                    self.springPoints.append(v)
                else:
                    self.newSpringPoints.append(v)
        for v in self.hx.vertices:
            if(v[3] == 0):
                if(v[2] not in self.newSpringPoints):
                    self.staticPoints.append(v[2])
        self.updatePlot()

    def getClosestVert(self,v, movable=1):   
        dists=[(np.sqrt(np.power(v[0] - v2[0], 2) + np.power(v[1] - v2[1], 2)), v2[2]) for v2 in self.hx.vertices if v2[3] == movable]
        if(len(dists) == 0):
            return -1
        return min(dists,key=itemgetter(0))[1]  


class getCellTypes:

    def __init__(self, hx, pc, fig, ax):
        self.finished = False
        self.celltype = None
        self.hx = hx
        self.collection = pc
        self.fig = fig
        self.canvas = fig.canvas
        self.ax = ax
        self.cid1 = self.canvas.mpl_connect('key_press_event', self.on_key)
        self.cid2 = self.canvas.mpl_connect('button_press_event', self.onpress)
        print("H to set Hinge; B to set Blade; V to set Veins (recommended last)")
        plt.show()

    def on_key(self, event):

        print('you pressed: ', event.key, event.xdata, event.ydata)
        if(event.key in "Hh" and not self.finished and not self.canvas.widgetlock.locked()):
            self.celltype = mhg.hingetype
            print("Drawing Hinge")
        elif(event.key in "Vv" and not self.finished and not self.canvas.widgetlock.locked()):
            self.celltype = mhg.veintype
            print("Drawing Vein")
        elif(event.key in "Bb" and not self.finished and not self.canvas.widgetlock.locked()):
            self.celltype = mhg.bladetype
            print("Drawing Blade")
        elif(event.key in "Mm" and not self.finished and not self.canvas.widgetlock.locked()):
            self.finished = True
            self.celltype = None
        elif(self.finished):
            self.fig.canvas.mpl_disconnect(self.cid1)
            self.fig.canvas.mpl_disconnect(self.cid2)
            plt.close()
            return
    def onpress(self, event):
        print("Starting Lasso")

        if self.canvas.widgetlock.locked():
            return
        if event.inaxes is None:
            return
        if self.celltype is None:
            return
        # acquire a lock on the widget drawing
        self.lasso = Lasso(event.inaxes,
                           (event.xdata, event.ydata),
                           self.callback)
        self.canvas.widgetlock(self.lasso)
        print("Finishing Lasso")
    def callback(self, verts):
        newtypes = []
        try:
            lim = Polygon(verts)
            hxpol = self.hx.getShapelyPolygons()
            for i, cell in enumerate(hxpol):
                if(lim.contains(cell)):
                    newtypes.append((i, self.changeType(self.hx.celltypes[i])))
                    #self.hx.celltypes[i] = self.changeType(self.hx.celltypes[i])
                elif(lim.overlaps(cell)):
                    if(lim.intersection(cell).area/cell.area >= 0.5):
                        newtypes.append((i, self.changeType(self.hx.celltypes[i])))
                        #self.hx.celltypes[i] = self.changeType(self.hx.celltypes[i])
        except:
            print("Error: Error drawing shape")
            self.canvas.widgetlock.release(self.lasso)
            del self.lasso
            return
        for i, t in newtypes:
            self.hx.celltypes[i] = t
        self.canvas.widgetlock.release(self.lasso)
        del self.lasso
        #fig, ax, pc = self.hx.plotHex2(fig=(self.fig, self.ax))
        self.collection.set_facecolors(self.hx.getColors())
        self.canvas.draw_idle()
        print("H to set Hinge; B to set Blade; V to set Veins (recommended last)")

    def changeType(self, t):   
        newt = self.celltype  
        if(newt == mhg.veintype):
            if(t == mhg.hingetype or t == mhg.veinhinge):
                newt = mhg.veinhinge
        return newt



###

###



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

   

class getPointsBorder():
    def __init__(self, fig, ax):
        self.points = []
        self.lines = []
        self.num_points = 0
        self.finished = False
        self.fig = fig
        self.ax = ax

        self.cid1 = self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.cid2 = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        plt.show()

    def onclick(self, event):
        if(not self.finished):
            print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
                   ('double' if event.dblclick else 'single', event.button,
                   event.x, event.y, event.xdata, event.ydata))
            self.points.append([event.xdata, event.ydata])
            self.ax.scatter([event.xdata], [event.ydata], color = "blue")
            if(self.num_points > 0):
                self.lines.append([self.num_points - 1, self.num_points])
                self.ax.plot([self.points[self.num_points - 1][0], self.points[self.num_points][0]], [self.points[self.num_points - 1][1], self.points[self.num_points][1]], color = "black")
            self.fig.canvas.draw()
            self.num_points += 1
    def on_key(self, event):
        print('you pressed', event.key, event.xdata, event.ydata)
        if(event.key in "mM" and not self.finished):
            self.finished = True
            self.fig.canvas.mpl_disconnect(self.cid2)
            self.lines.append([self.num_points - 1, 0])
            self.ax.plot([self.points[self.num_points - 1][0], self.points[0][0]], [self.points[self.num_points - 1][1], self.points[0][1]], color = "black")
            self.fig.canvas.draw()
        elif(self.finished):
            self.fig.canvas.mpl_disconnect(self.cid1)
            plt.close()
            return
    def getData(self):
        return (self.num_points, self.points, self.lines, self.fig, self.ax)


#makeWingVoronoi()
def main():
    makeWingHex(sys.argv[1])

if __name__ == '__main__':
    main()



