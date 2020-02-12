#import sys
import os
import argparse

import numpy as np
import matplotlib.pyplot as plt


cell_ext = '.cells'
vert_ext = '.points'
cent_ext = '.cent'
spring_ext = '.spr'

size = 2.5
nrow = 10
ncol = 10
noise = 0.1

bladetype = 0
hingetype = 1
veintype = 2
veinhinge = 3
wingcols = ['blue', 'green', 'black', 'gray']
springcols = ['red', 'yellow', 'purple', 'pink', 'brown']

DEFAULT_SPRING_KIND = 0

class HexGrid:
    """
    A class used to represent a squared Hexagonal Grid of nr x nc hexagons.
    Contains all data necessary to run VertexSystem

    AttributesaddStre
    ----------
    s : float
        Size of each hexagon, from center to vertex.
    nr : int
        Number of rows
    nc : int
        Number of columns
    centers : list of tuples of size 4:
        - row index of hexagon
        - col index of hexagon
        - x coordinate
        - y coordinate
    vertices :  list of lists of size 4:
        - x coordinate
        - y coordinate
        - index of vertex
        - movable (1) or not movable (0)        
    cells : list of lists of integers
        Each list has the ordered vertices of a cell
    springs : list of tuples of 2 integers
        Each tuple contains vertices forming the spring
    celltypes : list of integers representing cell types:
        bladetype = 0
        hingetype = 1
        veintype = 2
        veinhinge = 3
    vnum : int
        Number of vertices in grid
    ran : float
        Proportion of size (s) that each vertex is moved at random 
        to make the grid irregular
    pull : float
        Angle between vertices 0 and 6 (vertical hexagon sides). 
        For a regular hexagon the value is 60. 
    strecht : float
        Horizontal scaling of the grid. If greater than 1, grid 
        is stretched horizontally.
        If lower, it is be stretched vertically. 
    rotate : float
        Angles that the whole grid is rotated
    static_vertices : list of integers
        List of vertices that will not move in Vertex Model
    spring_vertices : list of integers
        Vertices that are not part of any cell but touch a spring
    spring_length : float
        Length of springs
    hingelimit
        Cells in columns <= hingelimit will be of cell type hinge
    veinpos : List of integers
        Rows of cells that are of type vein
    outname : str
        Identifier of grid

    Methods
    -------
    getData
    getPoint
    getPoint2
    rotateVertices
    getCells
    addNoise
    plotHex
    writeGrid
    trimStaticPoints
    addStatic
    addCellType
    addStrecht

    """

    width = lambda size: np.sqrt(3) * size
    height = lambda size: 2 * size

    def __init__(self, **kwargs):
        self.__args = kwargs
        self.s = kwargs['Size']
        self.nr = kwargs['Rows']
        self.nc = kwargs['Cols']
        self.ran = kwargs['Noise']
        self.pull = kwargs['Pull']
        self.strecht = kwargs['Strecht']
        self.rotate = kwargs['Rotate']

        self.static_vertices = kwargs['StaticVertices']
        self.spring_vertices = kwargs['Springs']
        self.spring_length = kwargs['SpringLength']

        self.hingelimit = kwargs['Hinge']
        self.veinpos = kwargs['Veins']
        self.outname = kwargs['Outname']
        self.inputname = kwargs['Read']

        self.centers = []
        self.vertices = []
        self.cells = []
        self.springs = []
        self.celltypes = []
        self.vnum = 0
        self.plotcol = None

        if(kwargs['Read'] != ""):
            self.loadData() 
        else:
            self.getCells()

            if(self.ran > 0):
                self.addNoise()
            if(self.strecht > 0):
                self.addStrecht()
            if(self.static_vertices != ''):
                self.addStatic()
            if(self.hingelimit != -1 or self.veinpos != ''):
                self.addCellType()
            else:
                self.celltypes = [bladetype for c in self.cells]

            if(self.rotate > 0):
                self.rotateVertices()

    def __iter__(self):
        for c in self.cells:
            yield c
    def __getitem__(self, ind):
        return self.cells[ind]
    def getData(self):
        return (self.centers, self.vertices, self.cells, self.springs, self.celltypes, self.outname)

    def loadData(self):
        self.vnum, self.vertices = self.readPointsFile()
        #print("verts", self.vertices)
        self.cells, self.celltypes = self.readCellsFile()
        #print("CELLS: ", self.cells)
        self.springs = self.readSprings()
        #print("sprs", self.springs)
        self.centers = self.readCenters()

    def readPointsFile(self):
        pointsFile = open(self.inputname + ".points", "r")
        numPoints = int(pointsFile.readline())

        # This is the list of points
        pointsList = []#np.zeros((numPoints,2))
        for row in range(numPoints):
            ll = pointsFile.readline().split()
            pointsList.append([float(ll[0]), float(ll[1]), int(ll[2]), int(ll[3])])
        pointsFile.close()
        return (numPoints, pointsList)

    def readCellsFile(self):
        cellsFile = open(self.inputname + ".cells", "r")
        numCells, numVerticesCell = [int(each) for each in cellsFile.readline().split()]

        polygonList = []
        celltypes = []
        for row in range(numCells):
            polygonIndex = [int(indPoint) for indPoint in cellsFile.readline().split() if int(indPoint) >= 0]
            celltypes.append(polygonIndex.pop()) #CAREFUL: if cells file does not have type, will produce error
            polygonList.append(polygonIndex)
        cellsFile.close()
        return (polygonList, celltypes)

    def readSprings(self):
        sprList = [] #just in case
        numsprings = 0
        try:
            springsfile = open(self.inputname + ".spr", "r")
            numsprings = int(springsfile.readline())
            sprList = [[int(j) for j in springsfile.readline().split('\t')] for i in range(numsprings) ]
            springsfile.close()
        except:
            print("no springs (.spr) file")
        return sprList
    def readCenters(self):  
        try:    
            centFile = open(self.inputname + ".cent", "r")   
            centers = [[int(c[0]), int(c[1]), float(c[2]), float(c[3])] for c in centFile.readline().split()]
        except: 
            print("no centers (.cent) file")
            centers = []
        return centers
    @staticmethod
    def getPoint(x, y, size, i):
        angle_deg = -1*60*(i + 1) + 30
        angle_rad = np.pi/180 * angle_deg
        return (x + size*np.cos(angle_rad), y + size*np.sin(angle_rad))

    @staticmethod
    def getPoint2(x, y, size, i, pull):
        a2= 60 + (60-pull)/2
        angles = [pull, a2, a2, pull, a2, a2]
        angle_deg = pull/2
        for side in range(i+1):
            angle_deg -= angles[side]
        angle_rad = np.pi/180 * angle_deg
        return (x + size*np.cos(angle_rad), y + size*np.sin(angle_rad))

    def rotateVertices(self):
        xmin, ymin = self.getMin()
        xmax, ymax = self.getMax()
        xorigin = 0.5*(xmin + xmax)
        yorigin = 0.5*(ymin + ymax)
        angle = self.rotate*np.pi/180
        for v in range(len(self.vertices)):
            dx = self.vertices[v][0] - xorigin
            dy = self.vertices[v][1] - yorigin
            xnew = dx*np.cos(angle) - dy*np.sin(angle) + xorigin
            ynew = dy*np.cos(angle) + dx*np.sin(angle) + yorigin
            self.vertices[v][0] = xnew 
            self.vertices[v][1] = ynew
        for c in range(len(self.centers)):
            dx = self.centers[c][2] - xorigin
            dy = self.centers[c][3] - yorigin
            xnew = dx*np.cos(angle) - dy*np.sin(angle) + xorigin
            ynew = dy*np.cos(angle) + dx*np.sin(angle) + yorigin
            newtup = (self.centers[c][0], self.centers[c][1],xnew, ynew)
            self.centers[c] = newtup 

    def getCells(self):
        nrow = self.nr
        ncol = self.nc
        size = self.s
        pull = self.pull
        #max_coord = self.__class__.height(size)*(ncol+1)
        self.centers = [ (j-1, i-1, self.__class__.width(size)*i + self.__class__.width(size)*0.5*(j%2) , self.__class__.height(size)*0.75*j)  for j in range(1, nrow + 1) for i in range(1, ncol + 1)]
        for i, j, x, y in self.centers:
            cell = []
            for v in range(6):
                newv = -1
                if(j > 0):
                    if(v == 2):
                        newv = self.cells[ncol*i + j - 1][0]
                    elif(v == 3):
                        newv = self.cells[ncol*i + j - 1][5]
                if(i > 0):
                    if(i%2 == 1):
                        if (v == 1):
                            newv = self.cells[ncol*(i-1) + j][3]
                        elif(v == 0):
                            newv = self.cells[ncol*(i-1) +j ][4]
                    else:
                        if (v == 2):
                            newv = self.cells[ncol*(i-1) + j][4]
                        elif(v == 1):
                            newv = self.cells[ncol*(i-1) +j ][5]   
                        elif(v == 0 and j < ncol - 1):
                            newv = self.cells[ncol*(i-1) +j  + 1][4]
                if(newv == -1):
                    newv = self.vnum
                    if(pull > 0):
                        a, b = self.getPoint2(x, y, size, v, pull)
                    else:
                        a, b = self.getPoint(x, y, size, v)
                    self.vertices.append([a, b, newv, 1]) # coord_x, coord_y, index, movable
                    self.vnum+=1
                cell.append(newv)
            self.cells.append(cell)

    def addNoise(self):
        size = self.s
        noise = self.ran
        for i in range(len(self.vertices)):
            self.vertices[i][0] = self.vertices[i][0] + np.random.uniform(-1*size*noise, size*noise)
            self.vertices[i][1] = self.vertices[i][1] + np.random.uniform(-1*size*noise, size*noise)

    def plotHex2(self, fig=None, save=False, alpha=0.4):
        from matplotlib.collections import PolyCollection
        if(fig is None):
            fig, ax = plt.subplots()
            fig.suptitle(self.outname, fontsize=16)
        else:
            fig, ax = fig
        ax.set_aspect('equal')
        col = np.zeros((len(self.cells), 6, 2))
        for c, cell in enumerate(self.cells):
            for v, vert in enumerate(cell):
                col[c, v, :] = (self.vertices[vert][0], self.vertices[vert][1])
        print(col.shape)
        pc = PolyCollection(col, facecolors= [wingcols[self.celltypes[j]] for j in range(len(self.cells))], alpha = alpha)
        pc.set_edgecolors("black")
        self.plotcol = ax.add_collection(pc)
        ax.autoscale_view()
        for c in self.springs:
            sscolor = "red" if len(c) < 3 else springcols[c[2]]
            plt.plot( [self.vertices[c[0]][0], self.vertices[c[1]][0]] ,  [self.vertices[c[0]][1], self.vertices[c[1]][1]], c = sscolor)

        if(save):
            self.saveFig()
        return(fig, ax, pc)

    def getPolCollection(self):
        from matplotlib.collections import PolyCollection
        fig, ax = plt.subplots()
        col = np.zeros((len(self.cells), 6, 2))
        for c, cell in enumerate(self.cells):
            for v, vert in enumerate(cell):
                col[c, v, :] = (self.vertices[vert][0], self.vertices[vert][1])
        print(col.shape)
        pc = PolyCollection(col, facecolors= [wingcols[self.celltypes[j]] for j in range(len(self.cells))], alpha = 0 )
        return pc

    def getColors(self):
        return [wingcols[self.celltypes[j]] for j in range(len(self.cells))]
    def plotHex(self, save=True, plot=None):
        if(plot is None):
            f, ax = plt.subplots()
        else:
            f, ax = plot
        #ax.scatter([i[2] for i in self.centers], [i[3] for i in self.centers], c = [wingcols[k] for k in self.celltypes])

        if(len(self.centers) < 50):
            for i in range(len(self.centers)):
                ax.annotate(i, [self.centers[i][2], self.centers[i][3]])

        for j, c in enumerate(self.cells):
            for i in range(len(c)):
                plt.plot([self.vertices[c[i]][0], self.vertices[c[(i+1)%len(c)]][0]] ,  [self.vertices[c[i]][1], self.vertices[c[(i+1)%len(c)]][1]], c = wingcols[self.celltypes[j]])

        for c in self.springs:
            plt.plot( [self.vertices[c[0]][0], self.vertices[c[1]][0]] ,  [self.vertices[c[0]][1], self.vertices[c[1]][1]], c = springcols[c[2]] )

        vplotx = []
        vploty = []
        for v in self.vertices:
            if(v[3] == 0):
                vplotx.append(v[0])
                vploty.append(v[1])
        ax.scatter(vplotx, vploty, c = "red")

        if(save):
            self.saveFig()
        return (f, ax)

    def saveFig(self):
        dirname = ''
        if(not os.path.isdir(self.outname)):
            try:
                os.mkdir(self.outname)
                dirname = self.outname + '/'
            except OSError:
                print ("Creation of the directory %s failed" % self.outname)
        else:
            dirname = self.outname + '/'
        plt.savefig(dirname + self.outname + '.png')
        plt.savefig(dirname + self.outname + '.svg', format='svg', dpi=1200)
        plt.show()


    def writeGrid(self, add=''):
        dirname = ''
        if(not os.path.isdir(self.outname)):
            try:
                os.mkdir(self.outname)
                dirname = self.outname + '/'
            except OSError:
                print ("Creation of the directory %s failed" % self.outname)
        else:
            dirname = self.outname + '/'
        f = open(dirname + self.outname + add + vert_ext,  'w')
        f.write(str(len(self.vertices)))
        f.write('\n')
        for v in self.vertices:
            f.write('\t'.join([str(i) for i in v]))
            f.write('\n')
        f.close()
        f = open(dirname + self.outname + add +  cell_ext,  'w')
        f.write(str(len(self.cells)) + "\t" + str(6))
        f.write('\n')
        for i, c in enumerate(self.cells):
            f.write('\t'.join([str(i) for i in c]))
            f.write('\t-999\t' + str(self.celltypes[i]) + '\n')
        f.close()
        f = open(dirname + self.outname + add +  cent_ext,  'w')
        for c in self.centers:
            f.write('\t'.join([str(i) for i in c]))
            f.write('\n')
        f.close()
        f = open(dirname + self.outname + add + spring_ext,  'w')
        f.write(str(len(self.springs)) + '\n')
        for c in self.springs:
            f.write('\t'.join([str(i) for i in c]))
            f.write('\n')
        f.close()

    @staticmethod
    def trimStaticPoints(spring_vertices):
        ll = [str.split(side, ',') for side in str.split(spring_vertices, ';') ]
        llnum = []
        for l in ll:
            lnum = []
            for el in l:
                if(':' in el):
                    rr = str.split(el, ':')
                    for i in range(int(rr[0]) , int(rr[1])):
                        lnum.append(i)
                elif ('' != el):
                    lnum.append(int(el))
            llnum.append(lnum)
        return llnum
   

    def addStatic(self):
        spring_length = self.spring_length
        nr = self.nr
        nc = self.nc
        left, top, bottom, right = self.trimStaticPoints(self.static_vertices)
        lefts, tops, bottoms, rights = self.trimStaticPoints(self.spring_vertices)
        #static vertices in cells (only in borders)
        for i in left:
            cnum  = self.cells[nc*i]
            self.vertices[cnum[2]][3] = 0
            self.vertices[cnum[3]][3] = 0
        for i in top:
            cnum  = self.cells[nc*(nr-1) + i]
            self.vertices[cnum[3]][3] = 0
            self.vertices[cnum[4]][3] = 0
            self.vertices[cnum[5]][3] = 0
        for i in bottom:
            cnum  = self.cells[i]
            self.vertices[cnum[0]][3] = 0
            self.vertices[cnum[1]][3] = 0
            self.vertices[cnum[2]][3] = 0
        for i in right:
            cnum  = self.cells[nc*i + nc - 1]
            self.vertices[cnum[0]][3] = 0
            self.vertices[cnum[5]][3] = 0
        #vertices joined to springs               
        for i in lefts:
            cnum  = self.cells[nc*i]
            if(i == 0):
                addspringto = [2]
            elif(i == nr-1 and i%2 == 0):
                addspringto = [3]
            elif(i%2 == 1):
                addspringto = [2, 3]
            else:
               addspringto = []
            for j in addspringto:
                v_to_fix = cnum[j]
                x =  self.vertices[v_to_fix][0] - spring_length  
                y =  self.vertices[v_to_fix][1]
                ind = len(self.vertices)  
                self.vertices.append([x, y, ind, 0])
                self.springs.append((v_to_fix, ind, DEFAULT_SPRING_KIND))
        for i in tops:
                cnum  = self.cells[nc*(nr-1) + i]
                v_to_fix = cnum[4]
                x =  self.vertices[v_to_fix][0] 
                y =  self.vertices[v_to_fix][1] + spring_length  
                ind = len(self.vertices)
                self.vertices.append([x, y, ind, 0])
                self.springs.append((v_to_fix, ind, DEFAULT_SPRING_KIND))
        for i in bottoms:
                cnum  = self.cells[i]
                v_to_fix = cnum[1]
                x =  self.vertices[v_to_fix][0] 
                y =  self.vertices[v_to_fix][1] - spring_length  
                ind = len(self.vertices)
                self.vertices.append([x, y, ind, 0])
                self.springs.append((v_to_fix, ind, DEFAULT_SPRING_KIND))
        for i in rights:
                cnum  = self.cells[nc*i + nc - 1]
                if(i == nr-1 and i%2 == 1):
                        addspringto = [5]
                elif(i%2 == 0):
                        addspringto = [0, 5]
                else:
                   addspringto = []
                for j in addspringto:
                        v_to_fix = cnum[j]
                        x =  self.vertices[v_to_fix][0] + spring_length  
                        y =  self.vertices[v_to_fix][1]
                        ind = len(self.vertices)  
                        self.vertices.append([x, y, ind, 0])
                        self.springs.append((v_to_fix, ind, DEFAULT_SPRING_KIND))

    def addSpring(self, v, x, y, kind=DEFAULT_SPRING_KIND):
        ind = len(self.vertices)
        self.vertices.append([x, y, ind, 0])
        self.springs.append([v, ind, kind]) 
        self.vnum+=1
        return ind       

    def getSpringWithVertex(self, v):
        res = -1
        for i, s in enumerate(self.springs):
            if(v in s):
                res = i
                break
        return res

    def addCellType(self):
        nr, nc = (self.nr, self.nc)
        if(self.veinpos != ''):
            self.veinpos = [int(i) for i in str.split(self.veinpos, ',')]
        else:
            self.veinpos = [nr + 1]
        for i, cell in enumerate(self.cells):
            col = i%nc
            row = i//nc
            if(row in self.veinpos):
                if(col > self.hingelimit):
                    self.celltypes.append(veintype)
                else:
                    self.celltypes.append(veinhinge)
            elif(col <= self.hingelimit):
                self.celltypes.append(hingetype)
            else:
                self.celltypes.append(bladetype)

    def addStrecht(self):
        strecht = self.strecht
        for v in range(len(self.vertices)):
           self. vertices[v][0] = self.vertices[v][0]*strecht
        for c in range(len(self.centers)):
            self.centers[c] = (self.centers[c][0], self.centers[c][1], self.centers[c][2]*strecht, self.centers[c][3])

    def removeCell(self, c):
        cell_to_remove = self.cells[c]
        self.cells.pop(c)
        self.centers.pop(c)
        self.celltypes.pop(c)
        v_to_remove = []
        for v in cell_to_remove:
            cont = 0
            for x in self.cells:
                if(v in x):
                    cont+=1
                    break
            if(cont == 0):
                v_to_remove.append(v)
        v_to_remove.sort(reverse=False)
        while(v_to_remove):
            v = v_to_remove.pop()
            spring_vertices = self.removeSpringsWithVert(v)
            #self.remove_spring_vertices(spring_vertices)
            spring_vertices = self.removeVertex(v, spring_vertices)
            if(spring_vertices):
                v_to_remove.extend(spring_vertices)
                v_to_remove.sort(reverse=False)

    def removeVertex(self, v, list_to_update):
        for vert in range(v, len(self.vertices) - 1):
            newvert = self.vertices[vert  + 1][2]
            self.vertices[vert + 1][2] = vert
            self.vertices[vert] = self.vertices[vert + 1]
            for cind, cell in enumerate(self.cells):
                for vicell, vcell in enumerate(cell):
                    if(vcell == newvert):
                        self.cells[cind][vicell] = vert
                        break
            for sind, spring in enumerate(self.springs):
                for vispring, vspring in enumerate(spring):
                    if(vspring == newvert):
                        self.springs[sind] = (vert, spring[1]) if vispring == 0 else (spring[0], vert)
                        break
            for vvind, vv in enumerate(list_to_update):
                if(vv == newvert):
                    list_to_update[vvind] = vert
        self.vertices.pop()  
        self.vnum -= 1
        return list_to_update

    def removeSpringVertices(self, s, v):         
        #v.sort(reverse=True) #Assumes they are sorted
        #s.sort(reverse=True)
        for vert in v:
            self.removeVertex(vert, [])
        for spr in s:
            self.springs.pop(spr)

    
    def removeSpringsWithVert(self, v):
        springs_to_remove = []
        for ind, s in enumerate(self.springs):
            if(v in s):
                springs_to_remove.append(ind)
        other_vertices_to_remove = []
        springs_to_remove.sort(reverse = False)
        for ind in springs_to_remove:
            if(self.springs[ind][0] == v):
                other_vertices_to_remove.append( self.springs[ind][1])
            else:
                 other_vertices_to_remove.append(self.springs[ind][0])
            self.removeSpring(ind)
        return other_vertices_to_remove

    def removeSpring(self, ind):      
        self.springs.pop(ind)  

    def getShapelyPolygons(self):
        from shapely.geometry import Polygon
        pols = []
        for c in self.cells:
            pols.append(Polygon( [(self.vertices[p][0], self.vertices[p][1]) for p in c] ))   
        return pols 

    def getMax(self):
        x = []
        y = []
        for v in self.vertices:
            x.append(v[0])
            y.append(v[1])
        return (max(x), max(y))

    def getMin(self):
        x = []
        y = []
        for v in self.vertices:
            x.append(v[0])
            y.append(v[1])
        return (min(x), min(y))

    def copyModifiedStructure(self, points, lines):
        from shapely.geometry import Polygon
        hxpol = self.getShapelyPolygons()
        lim = Polygon(points)
        cells_to_use = []
        cells_to_fit = []
        for i, cell in enumerate(hxpol):
            if(lim.contains(cell)):
                cells_to_use.append(i)
            elif(lim.overlaps(cell)): 
                arearatio = lim.intersection(cell).area/cell.area
                if(arearatio > 0.2):
                        cells_to_fit.append(i)
                if(arearatio > 0.7):
                        cells_to_use.append(i)
                    #pass
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
            c = self.cells[cind] 
            for i, v in enumerate(c):
                if(v in vertDict):
                    c[i] = vertDict[v]
                else:
                    vertDict.setdefault(v, vnum)
                    c[i] = vnum
                    vertices.append( self.vertices[v] )
                    vertices[vnum][2] = vnum
                    vnum+=1                
            cells.append(c)
            #print([vertices[vv] for vv in c])
            celltypes.append(self.celltypes[cind])
            centers.append(self.centers[cind])
        '''
        for cind in cells_to_fit:
            intersection = lim.intersection(hxpol[cind])
            c = self.cells[cind] 
            newcell = []
            for i, v in enumerate(c):

                pcell = Point(self.vertices[v][0], self.vertices[v][1])
                print("PCELL: ", pcell)
                if(intersection.touches(pcell)): ##Add this point to cell
                    if(v in vertDict):
                        newcell.append(vertDict[v])
                    else:
                        vertDict.setdefault(v, vnum)
                        newcell.append(vnum)
                        vertices.append( self.vertices[v] )
                        vertices[vnum][2] = vnum
                        vnum+=1   
            if(len(newcell) > 3):             
                cells.append(newcell)
                #print([vertices[vv] for vv in newcell])
                celltypes.append(self.celltypes[cind])
                centers.append(self.centers[cind])
            '''        
    
        for s in self.springs:
            isthere = 0 if s[0] in vertDict else 1 if s[1] in vertDict else -1
            if(isthere == 0):
                newspring = (vertDict[s[0]], vnum)
                vertDict.setdefault(s[1], vnum)
                vertices.append( self.vertices[s[1]] )
                vertices[vnum][2] = vnum
                vnum+=1                
            elif(isthere == 1):
                newspring = (vnum, vertDict[s[1]])
                vertDict.setdefault(s[0], vnum)
                vertices.append( self.vertices[s[0]] )
                vertices[vnum][2] = vnum
                vnum+=1    
            springs.append(newspring)     
        self.cells = cells
        self.vertices = vertices
        self.springs = springs
        self.celltypes = celltypes
        self.centers = centers
        #self.adaptToBorders(points, lines)#Doesnt give good results
        return self    

    def adaptToBorders(self, points, lines, precision=100):
        print(lines)
        dist = lambda x1, x2, y1, y2: ((x1-x2)**2 + (y1-y2)**2)**0.5
        for a, b in lines:
            sx = np.linspace(points[a][0], points[b][0], precision)
            sy = np.linspace(points[a][1], points[b][1], precision)
            vert = -1 
            linepoints = []
            for p in range(precision):
                distances = [dist(sx[p], v[0], sy[p], v[1]) for v in self.vertices]
                closest = np.argmin(distances)
                if(closest == vert):
                    linepoints.append(p)
                elif(vert == -1):
                    vert = closest
                    linepoints.append(p)
                else: 
                    self.vertices[vert][0] = np.mean(sx[linepoints])
                    self.vertices[vert][1] = np.mean(sy[linepoints])
                    linepoints = [p]
                    vert = closest
    def getBorderPoints(self):
        touchedcells = np.zeros(len(self.vertices))
        for c in self.cells:
            for v in c:
                touchedcells[v] += 1
        return [i for i in range(len(self.vertices)) if touchedcells[i] > 0 and touchedcells[i] < 3]


def getArgDict(args):
        argdict = {}
        argdict.setdefault('Size', args.Size)
        argdict.setdefault('Rows', args.Rows)
        argdict.setdefault('Cols', args.Cols)
        argdict.setdefault('Noise', args.Noise)
        argdict.setdefault('Pull', args.Pull)
        argdict.setdefault('Strecht', args.Strecht)
        argdict.setdefault('Rotate', args.Rotate)

        argdict.setdefault('StaticVertices', args.StaticVertices)
        argdict.setdefault('Springs', args.Springs)
        argdict.setdefault('SpringLength', args.SpringLength)

        argdict.setdefault('Hinge', args.Hinge)
        argdict.setdefault('Veins', args.Veins)
        argdict.setdefault('Outname', args.Outname + '_s'+str(args.Size) + '_' + str(args.Rows) + 'x' + str(args.Cols) + '_n' + str(args.Noise))
        argdict.setdefault('Read', args.Read)

        print('size: ', args.Size, '; num. rows: ', args.Rows, '; num. cols: ', args.Cols, '; noise: ', args.Noise, '; strecht: ', args.Pull, '; output files: ', argdict['Outname'])

        return argdict


parser = argparse.ArgumentParser(description='Hexagonal grid arguments.')
parser.add_argument('-o', '--Outname', metavar='outname', type=str, default = "hexgrid", 
                                        help='Identifier. Used as prefix of all output files. ')
parser.add_argument('-i', '--Read', metavar='Read', type=str, default = "", 
                                        help='Identifier. Read hexagonal grid. ')
parser.add_argument('-s', '--Size', metavar='size', type=float, default = size,
                                        help=' Hexagon size: distance between center and vertices. ')
parser.add_argument('-r', '--Rows', metavar='rows', type=int, default = nrow,
                                        help='Number rows. ')
parser.add_argument('-c', '--Cols', metavar='cols', type=int, default = ncol,
                                        help='Number of columns.')
parser.add_argument('-n', '--Noise', metavar='noise', type=float, default = noise,
                                        help='Maximum proportion of size that each vertex is moved randomly')
parser.add_argument('-t', '--StaticVertices', metavar='static_vertices', type=str, default = "", 
                                        help='Vertices in margin that are static. Sytax: "left;top;bottom;right", where each position can be a series of numbers or ranges (eg. 0:5) separated by commas. Eg: "0,4:6;;;0:10,12,15"')
parser.add_argument('-p', '--Springs', metavar='springs', type=str, default = "", 
                                        help='Vertices in margin that can move but are attached to springs. Sytax: "left;top;bottom;right", where each position can be a series of numbers or ranges (eg. 0:5) separated by commas. Eg: "0,4:6;;;0:10,12,15"')
parser.add_argument('-l', '--SpringLength', metavar='spring_length', type=float, default = 3.5,
                                        help='Length of springs')
parser.add_argument('-g', '--Hinge', metavar='Hinge_pos', type=int, default = -1,
                                        help='Hexagons in columns in range [0:g] will be considered hinge (value 1).')
parser.add_argument('-v', '--Veins', metavar='Veins_pos', type=str, default = '',
                                        help='Hexagons in rows specified (as integers) will be considered veins ( value 2).')
parser.add_argument('-f', '--Pull', metavar='Change relative size of angles', type=float, default = 0,
                                        help='Value of side angle (default is 60 for a regular hexagon')
parser.add_argument('-k', '--Strecht', metavar='Strecht', type=float, default = 0,
                                        help='Scale horizontally by a factor of k')
parser.add_argument('-j', '--Rotate', metavar='Rotate', type=float, default = 0,
                                        help='Rotate the whole grid specified angles')

def main():
        args = parser.parse_args()
        argdict = getArgDict(args)
        hx = HexGrid(**argdict)
        hx.writeGrid()
        hx.plotHex()

if __name__ == '__main__':
        main()
