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

import make_hexagonal_grid as mhg


inputname = "wing.png"


def makeWingHex(args):
    # Make hexagonal grid
    hx, im = setHexagonSize(args).getGrid()
    fig, ax = plt.subplots()
    fig.suptitle(hx.outname, fontsize=16)
    maxx, maxy = hx.getMax()
    if(len(hx.inputname) > 0):
        maxx*=1.1
        maxy*=1.1
    ax.imshow(im, extent = [0, maxx, 0, maxy])
    hx.writeGrid(add='_a')

    # Draw border by hand 
    if(len(hx.inputname) == 0):
        border = getPointsBorder2(fig, ax, hx, im, maxx, maxy)
        border.disconnect()
        num_points, points, lines, fig, ax, hx = border.getData()
        plt.close()
        f, ax, pc = hx.plotHex2()
        ax.imshow(im, extent = [0, maxx, 0, maxy])
        border.correctBorderManually((f, ax))
        border.writePoints()
        hx.writeGrid(add='_b')
    else:
        plt.close()
    f, ax, pc = hx.plotHex2()
    ax.imshow(im, extent = [0, maxx, 0, maxy])
    #Set cell types: hinge and veins
    getCellTypes(hx, pc, f, ax)
    f, ax, pc = hx.plotHex2()
    ax.imshow(im, extent = [0, maxx, 0, maxy])
    hx.writeGrid(add='_c')

    #Set borders: static vertices and springs
    getBorders(hx, pc, f, ax)
    hx.writeGrid()
    f, ax, pc = hx.plotHex2(save=True)
    plt.show()


class setHexagonSize:

    def __init__(self, args):

        imfile = args.Inputimage 
        inputfile = args.Inputname
        print("IN: ", inputfile, ", len: ", len(inputfile))
        hxname = args.Outname
        
        self.im = plt.imread(imfile)
        self.imratio = self.im.shape[0]/self.im.shape[1]
        if(len(inputfile) > 0):
            self.im[:,:,:] = 255
            self.editing = True
        else:
            self.editing = False

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
        self.argdict.setdefault('Read', inputfile)
 
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

        print("Select parameters. Press 'm' to continue (Inactive for editing)")
        plt.show()

    def update(self, val):
        if(not self.editing):
            self.ax.clear()
            self.argdict['Rows'] = (self.scols.val*self.imratio*self.sstrecht.val).astype(int) 
            self.argdict['Cols']  =  self.scols.val.astype(int)
            self.argdict['Noise'] = self.snoise.val
            self.argdict['Strecht'] = self.sstrecht.val
            self.argdict['Rotate'] =  self.srotate.val
            self.makeAndAddToPlot() 


    def reset(self, event):
        if(not self.editing):
            self.argdict = self.defaultArgs.copy()
            self.scols.reset() 
            self.snoise.reset() 
            self.sstrecht.reset() 
            self.srotate.reset() 

    def makeAndAddToPlot(self):
        self.hx = mhg.HexGrid(**self.argdict)
        self.maxx, self.maxy = self.hx.getMax()
        if(self.editing):
            self.maxx*=1.1
            self.maxy*=1.1
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
        self.stringEdges = []
        self.springPlotRefs = dict()
        self.pointPlotRefs = dict()
        self.stringEdgesPlotRefs = dict()
        self.springType = 0
        self.hx = hx
        self.borderpoints = [self.hx.vertices[i] for i in self.hx.getBorderPoints()]
        self.collection = pc
        self.fig = fig
        self.canvas = fig.canvas
        self.ax = ax
        self.cid1 = self.canvas.mpl_connect('key_press_event', self.on_key)
        self.cid2 = self.canvas.mpl_connect('button_press_event', self.onpress)

        self.updateInnerLists()
        for sp in self.staticPoints:
            ref = self.ax.scatter(self.hx.vertices[sp][0], self.hx.vertices[sp][1], c="red")
            self.pointPlotRefs.setdefault(sp, ref)

        plt.subplots_adjust(bottom=0.2)
        self.axbox = plt.axes([0.1, 0.05, 0.8, 0.075])

        print("T to set static vertices; V to set springs; N to choose spring type; J to remove static vertices or springs; M to continue")
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
            print("Press M to continue")
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
        self.canvas.widgetlock.release(self.lasso)
        del self.lasso
        if(self.drawing == "static"):
            self.setStatic(verts)
        elif(self.drawing == "springs"):
            self.setSprings(verts)
        elif(self.drawing == "remove"):
            self.removeStatic(verts)
        #fig, ax, pc = self.hx.plotHex2(fig=(self.fig, self.ax))
        self.canvas.draw_idle()
        print("T to set static vertices; V to set springs; N to choose spring type; J to remove static vertices or springs; N to set spring type; M to continue")

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
            if(self.springType < 0 and self.stringEdges):
                self.stringEdges.append([self.stringEdges[-1][1], snum])
            elif(self.springType < 0):
                self.stringEdges.append([-1, snum])
        for a, b in self.stringEdges:
            self.hx.addStringEdge(a, b)
        self.updatePlot()

    def updatePlot(self):
        for sl in self.springPlotRefs.values():
            self.ax.lines.remove(sl)
        for sl in self.stringEdgesPlotRefs.values():
            self.ax.lines.remove(sl)
        self.springPlotRefs.clear()
        for i, s in enumerate(self.hx.springs):
            spref = self.ax.plot([self.hx.vertices[s[0]][0], self.hx.vertices[s[1]][0]], [self.hx.vertices[s[0]][1], self.hx.vertices[s[1]][1]], c=mhg.springcols[s[2]])
            self.springPlotRefs.setdefault(i, spref[0])
        for i, s in enumerate(self.hx.stringEdges):
            if(s[0] < 0):
                continue
            spref = self.ax.plot([self.hx.vertices[s[0]][0], self.hx.vertices[s[1]][0]], [self.hx.vertices[s[0]][1], self.hx.vertices[s[1]][1]], c=mhg.springcols[0])
            self.stringEdgesPlotRefs.setdefault(i, spref[0])           
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
        self.updateInnerLists()
        self.updatePlot()

    def updateInnerLists(self):
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

    def getClosestVert(self,v, movable=1, onlyBorder=True):
        if(onlyBorder):
            search=self.borderpoints
        else:
            search=self.hx.vertices   
        dists=[(np.sqrt(np.power(v[0] - v2[0], 2) + np.power(v[1] - v2[1], 2)), v2[2]) for v2 in search if v2[3] == movable]
        if(len(dists) == 0):
            return -1
        return min(dists,key=itemgetter(0))[1]  

class getCellTypes:

    def __init__(self, hx, pc, fig, ax):
        self.finished = False
        self.celltype = None
        self.moving_point = False
        self.point_selected = -1
        self.hx = hx
        self.collection = pc
        self.fig = fig
        self.canvas = fig.canvas
        self.ax = ax
        self.cid1 = self.canvas.mpl_connect('key_press_event', self.on_key)
        self.cid2 = self.canvas.mpl_connect('button_press_event', self.onpress)
        self.dist = lambda x1, x2, y1, y2: ((x1-x2)**2 + (y1-y2)**2)**0.5

        print("N to set Hinge; B to set Blade; V to set Veins (recommended last); Z to move a point; M to continue")
        plt.show()

    def on_key(self, event):

        print('you pressed: ', event.key, event.xdata, event.ydata)
        if(event.key in "Nn" and not self.finished and not self.canvas.widgetlock.locked()):
            self.celltype = mhg.hingetype
            print("Drawing Hinge. Drag mouse to draw.")
        elif(event.key in "Vv" and not self.finished and not self.canvas.widgetlock.locked()):
            self.celltype = mhg.veintype
            print("Drawing Vein. Drag mouse to draw.")
        elif(event.key in "Bb" and not self.finished and not self.canvas.widgetlock.locked()):
            self.celltype = mhg.bladetype
            print("Drawing Blade. Drag mouse to draw.")
        elif(event.key in "Mm" and not self.finished and not self.canvas.widgetlock.locked()):
            self.finished = True
            self.celltype = None
            print("Press M to continue")
        elif(event.key in "Zz" and not self.finished and not self.canvas.widgetlock.locked()):
            if(not self.moving_point): 
                self.moving_point = True
                print("Move point")
            else:
                self.moving_point = False
                print("Finished moving points")               
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
        if self.moving_point:
            self.movePoint(event)
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
        self.canvas.widgetlock.release(self.lasso)
        del self.lasso
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
            #self.canvas.widgetlock.release(self.lasso)
            #del self.lasso
            return
        for i, t in newtypes:
            self.hx.celltypes[i] = t
        #self.canvas.widgetlock.release(self.lasso)
        #del self.lasso
        #fig, ax, pc = self.hx.plotHex2(fig=(self.fig, self.ax))
        self.collection.set_facecolors(self.hx.getColors())
        self.canvas.draw_idle()
        print("N to set Hinge; B to set Blade; V to set Veins (recommended last); Z to move points")

    def changeType(self, t):   
        newt = self.celltype  
        if(newt == mhg.veintype):
            if(t == mhg.hingetype or t == mhg.veinhinge):
                newt = mhg.veinhinge
        return newt

    def movePoint(self, event):
        if(self.point_selected < 0):
            dists = [self.dist(event.xdata, v[0], event.ydata, v[1]) for v in self.hx.vertices] 
            self.point_selected = np.argmin(dists)

        else:
            self.hx.vertices[self.point_selected][0] = event.xdata
            self.hx.vertices[self.point_selected][1] = event.ydata 
            self.ax.collections.clear()
            self.fig, self.ax, self.collection = self.hx.plotHex2((self.fig, self.ax))
            self.point_selected = -1
            self.fig.canvas.draw_idle() 

###

###

class getPointsBorder2():
    def __init__(self, fig, ax, hx, im, maxx, maxy):
        self.points = []
        self.lines = []
        self.num_points = 0
        self.finished = False
        self.fig = fig
        self.ax = ax
        self.hx = hx
        self.im = im
        self.maxx = maxx
        self.maxy = maxy
        self.cid1 = None
        self.cid2 = None
        self.borderpoints = None
        self.borderscatter = None
        self.movePointBlocked = False

        self.poly = PolygonSelector(ax, self.onselect)
        self.dist = lambda x1, x2, y1, y2: ((x1-x2)**2 + (y1-y2)**2)**0.5
        self.dist_to_line = lambda x1, x2, px, y1, y2, py: abs( px*(y2-y1) - py*(x2-x1) - y2*x1 + x2*y1)/(((x2-x1)**2 + (y2-y1)**2)**0.5)
        self.closestInLine_x = lambda a, b, c, x, y: (b*(b*x - a*y) - a*c)/(a**2 + b**2)
        self.closestInLine_y = lambda a, b, c, x, y: (a*(-b*x + a*y) - b*c)/(a**2 + b**2)
        self.point_selected = -1
        self.point_in_border = -1
        self.finished = False
        print("Draw border with mouse: click in each point. Close polygon and close window to continue")
        plt.show()

    def getData(self):
        return (self.num_points, self.points, self.lines, self.fig, self.ax, self.hx)

    def onselect(self, verts):
        self.points = verts
        self.num_points = len(self.points)
        previous = self.num_points - 1
        for i, p in enumerate(verts):
            self.lines.append([previous, i])
            previous = i
        #self.poly.disconnect_events()
    def disconnect(self):
        self.poly.disconnect_events()
        self.hx = self.hx.copyModifiedStructure(self.points, self.lines)
        self.correctBorderAutomatically()
    def plotBorder(self, fig=None):
        if(fig is not None):
            self.fig, self.ax = fig
        self.ax.plot([self.points[len(self.points) - 1][0], self.points[0][0]], [self.points[len(self.points) - 1][1], self.points[0][1]], color = "black")
        for i in range(len(self.points) - 1):
            self.ax.plot([self.points[i][0], self.points[i+1][0]], [self.points[i][1], self.points[i+1][1]], color = "black")

    def correctBorderAutomatically(self):
        self.borderpoints = self.hx.getBorderPoints()
        for p in self.borderpoints:
            distances = np.zeros((len(self.lines) , 9))
            for i, l in enumerate(self.lines):
                p1, p2 = l
                distances[i, 0] = self.dist_to_line(self.points[p1][0], self.points[p2][0], self.hx.vertices[p][0], self.points[p1][1], self.points[p2][1], self.hx.vertices[p][1])
                distances[i, 1] = (self.points[p2][1] - self.points[p1][1])/(self.points[p2][0] - self.points[p1][0])#m
                distances[i, 2] = self.points[p2][1] - distances[i, 1]*self.points[p2][0]#n
                distances[i, 3] = -1*distances[i, 1] #a
                distances[i, 4] = 1 #b
                distances[i, 5] = -1*distances[i, 2] #c
                distances[i, 6] = self.closestInLine_x(distances[i, 3], distances[i, 4], distances[i, 5], self.hx.vertices[p][0], self.hx.vertices[p][1])#newx
                distances[i, 7] = self.closestInLine_y(distances[i, 3], distances[i, 4], distances[i, 5], self.hx.vertices[p][0], self.hx.vertices[p][1])#newy
                distances[i, 8] = self.putPointInside(distances[i, 6], distances[i, 7], self.points[p1][0], self.points[p1][1], self.points[p2][0], self.points[p2][1], self.hx.vertices[p][0], self.hx.vertices[p][1])
                #distances[i, 8] = self.dist(self.hx.vertices[p][0], self.points[p2][0], self.hx.vertices[p][1], self.points[p2][1])
            distances = distances[distances[:,0].argsort()]
            closestline = -1
            for i in range(len(self.lines)):
                if(distances[i, 8] > 0):
                    closestline = i 
                    break

            if(closestline > -1):
                self.hx.vertices[p][0] = distances[closestline, 6] 
                self.hx.vertices[p][1] = distances[closestline, 7]
            #print("Vertex %d: x=%f y=%f; newx:%f, newy:%f, y in eq: %f"%( p, self.hx.vertices[p][0], self.hx.vertices[p][1], newx, newy, m*newx+n))

    def putPointInside(self, x, y, x1, y1, x2, y2, oldx, oldy):
        if((x >= x1 and x <= x2) or (x >= x2 and x <= x1)):
            return 1
        else:
            return -1

    def correctBorderManually(self, fig):
        if(fig is not None):
            self.fig, self.ax = fig
        self.plotBorder(fig)
        self.borderpoints = self.hx.getBorderPoints()

        self.borderscatter = [self.ax.scatter(self.hx.vertices[p][0],  self.hx.vertices[p][1], color = "green") for p in self.borderpoints]
        
        self.ax.scatter([self.hx.vertices[p][0] for p in self.borderpoints],  [self.hx.vertices[p][1] for p in self.borderpoints], color = "green")
        self.cid1 = self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.cid2 = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        print("Select border points to move. Z to stop/resume selection (useful to zoom etc); M to exit")
        plt.show()

    def on_key(self, event):
        print('you pressed', event.key, event.xdata, event.ydata)
        if(event.key in "mM"):
            self.fig.canvas.mpl_disconnect(self.cid2)
            self.fig.canvas.mpl_disconnect(self.cid1)
            self.finished = True
            plt.close()
        elif(event.key in "zZ"):
            if(self.movePointBlocked):
                print("Border Point Selector Unblocked")
                self.movePointBlocked = False
            else:
                print("Border Point Selector Blocked")
                self.movePointBlocked = True
        print("Select border points to move. Z to stop/resume selection (useful to zoom etc); M to exit")
    def onclick(self, event):
        if event.inaxes is None:
            return
        if(self.movePointBlocked):
            return
        if(self.point_selected < 0):
            dists = [self.dist(event.xdata, self.hx.vertices[v][0], event.ydata, self.hx.vertices[v][1]) for v in self.borderpoints] 
            self.point_in_border = np.argmin(dists)
            self.point_selected = self.borderpoints[self.point_in_border]

        else:
            self.hx.vertices[self.point_selected][0] = event.xdata
            self.hx.vertices[self.point_selected][1] = event.ydata 
            self.ax.collections.clear()
            self.hx.plotHex2((self.fig, self.ax))
            self.borderscatter = [self.ax.scatter(self.hx.vertices[p][0],  self.hx.vertices[p][1], color = "green") for p in self.borderpoints]
            self.point_selected = -1
            self.fig.canvas.draw_idle()
        print("Select border points to move. Z to stop/resume selection (useful to zoom etc); M to exit")
    def writePoints(self):
        dirname = ''
        if(not os.path.isdir(self.hx.outname)):
            try:
                os.mkdir(self.hx.outname)
                dirname = self.hx.outname + '/'
            except OSError:
                print ("Creation of the directory %s failed" % self.hx.outname)
        else:
            dirname = self.hx.outname + '/'
        with open(dirname + self.hx.outname + '.border', 'w') as f:
            f.write('x\ty\n')
            for x, y in self.points:
                f.write(str(x) + '\t' + str(y) + '\n')


parser = argparse.ArgumentParser(description='Wing Drawer arguments.')
parser.add_argument('-o', '--Outname', metavar='outname', type=str, default = "hexgrid", 
                                        help='Identifier. Used as prefix of all output files. ')
parser.add_argument('-i', '--Inputname', metavar='inputname', type=str, default = "", 
                                        help='Identifier. Used as prefix to read files. ')
parser.add_argument('-w', '--Inputimage', metavar='inputimage', type=str, default = inputname, 
                                        help='Use this image to draw on it. ')
#makeWingVoronoi()
def main():
    args = parser.parse_args()
    makeWingHex(args)

if __name__ == '__main__':
    main()



