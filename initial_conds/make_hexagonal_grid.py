import sys
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

        self.centers = []
        self.vertices = []
        self.cells = []
        self.springs = []
        self.celltypes = []
        self.vnum = 0

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
        angle = self.rotate*np.pi/180
        for v in range(len(self.vertices)):
            xnew = self.vertices[v][0]*np.cos(angle) - self.vertices[v][1]*np.sin(angle)
            ynew = self.vertices[v][1]*np.cos(angle) + self.vertices[v][0]*np.sin(angle)
            self.vertices[v][0] = xnew 
            self.vertices[v][1] = ynew
        for c in range(len(self.centers)):
            xnew = self.centers[c][2]*np.cos(angle) - self.centers[c][3]*np.sin(angle)
            ynew = self.centers[c][3]*np.cos(angle) + self.centers[c][2]*np.sin(angle)
            newtup = (self.centers[c][0], self.centers[c][1],xnew, ynew)
            self.centers[c] = newtup 
    def getCells(self):
        nrow = self.nr
        ncol = self.nc
        size = self.s
        pull = self.pull
        max_coord = self.__class__.height(size)*(ncol+1)
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

    def plotHex(self, save=True, plot=None):
        if(plot is None):
            f, ax = plt.subplots()
        else:
            f, ax = plot
        ax.scatter([i[2] for i in self.centers], [i[3] for i in self.centers], c = [wingcols[k] for k in self.celltypes])

        if(len(self.centers) < 50):
            for i in range(len(self.centers)):
                ax.annotate(i, [self.centers[i][2], self.centers[i][3]])
        for j, c in enumerate(self.cells):
            for i in range(6):
                plt.plot([self.vertices[c[i]][0], self.vertices[c[(i+1)%6]][0]] ,  [self.vertices[c[i]][1], self.vertices[c[(i+1)%6]][1]], c = wingcols[self.celltypes[j]])
        for c in self.springs:
            plt.plot([self.vertices[c[0]][0], self.vertices[c[1]][0]] ,  [self.vertices[c[0]][1], self.vertices[c[1]][1]], c = "red")
        for v in self.vertices:
            if(v[3] == 0):
               ax.scatter(v[0], v[1], c = "red")
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
        plt.show()


    def writeGrid(self):
        dirname = ''
        if(not os.path.isdir(self.outname)):
            try:
                os.mkdir(self.outname)
                dirname = self.outname + '/'
            except OSError:
                print ("Creation of the directory %s failed" % self.outname)
        else:
            dirname = self.outname + '/'
        f = open(dirname + self.outname + vert_ext,  'w')
        f.write(str(len(self.vertices)))
        f.write('\n')
        for v in self.vertices:
            f.write('\t'.join([str(i) for i in v]))
            f.write('\n')
        f.close()
        f = open(dirname + self.outname + cell_ext,  'w')
        f.write(str(len(self.cells)) + "\t" + str(6))
        f.write('\n')
        for i, c in enumerate(self.cells):
            f.write('\t'.join([str(i) for i in c]))
            f.write('\t-999\t' + str(self.celltypes[i]) + '\n')
        f.close()
        f = open(dirname + self.outname + cent_ext,  'w')
        for c in self.centers:
            f.write('\t'.join([str(i) for i in c]))
            f.write('\n')
        f.close()
        f = open(dirname + self.outname + spring_ext,  'w')
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
                self.springs.append((v_to_fix, ind))
        for i in tops:
                cnum  = self.cells[nc*(nr-1) + i]
                v_to_fix = cnum[4]
                x =  self.vertices[v_to_fix][0] 
                y =  self.vertices[v_to_fix][1] + spring_length  
                ind = len(self.vertices)
                self.vertices.append([x, y, ind, 0])
                self.springs.append((v_to_fix, ind))
        for i in bottoms:
                cnum  = self.cells[i]
                v_to_fix = cnum[1]
                x =  self.vertices[v_to_fix][0] 
                y =  self.vertices[v_to_fix][1] - spring_length  
                ind = len(self.vertices)
                self.vertices.append([x, y, ind, 0])
                self.springs.append((v_to_fix, ind))
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
                        self.springs.append((v_to_fix, ind))

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
                pass
                v_to_remove.append(v)
        v_to_remove.sort(reverse=True)
        while(v_to_remove):
            v = v_to_remove.pop()
            spring_vertices = self.removeSpringsWithVert(v)
            self.remove_spring_vertices(spring_vertices)
            v_to_remove = self.removeVertex(v)

    def removeVertex(self, v):
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
        self.vertices.pop()  
        self.vnum -= 1

    def remove_spring_vertices(self, spring_vertices):
        spring_vertices.sort(reverse=True)
        for v in spring_vertices:
            for vert in range(v, len(self.vertices) - 1):
                self.vertices[vert + 1][2] = self.vertices[vert][2]
                self.vertices[vert] = self.vertices[vert + 1]
            self.vertices.pop()
            self.vnum -= 1
               
    def removeSpringsWithVert(self, v):
        springs_to_remove = []
        for ind, s in enumerate(self.springs):
            if(v in s):
                springs_to_remove.append(ind)
        other_vertices_to_remove = []
        springs_to_remove.sort(reverse = True)
        for ind in springs_to_remove:
            other_vertices_to_remove.append( self.springs[ind][1] if(self.springs[ind][0] == v) else self.springs[ind][0] )
            self.removeSpring(ind)
        return other_vertices_to_remove
    def removeSpring(self, ind):      
        self.springs.pop(ind)        

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

        print('size: ', args.Size, '; num. rows: ', args.Rows, '; num. cols: ', args.Cols, '; noise: ', args.Noise, '; strecht: ', args.Pull, '; output files: ', argdict['Outname'])

        return argdict

parser = argparse.ArgumentParser(description='Hexagonal grid arguments.')
parser.add_argument('-o', '--Outname', metavar='outname', type=str, default = "hexgrid", 
                                        help='Identifier. Used as prefix of all output files. ')
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
        hx.removeCell(0)
        hx.writeGrid()
        hx.plotHex()
        hx.removeCell(1)
        hx.writeGrid()
        hx.plotHex()
        hx.removeCell(2)
        hx.writeGrid()
        hx.plotHex()
        hx.removeCell(20)
        hx.writeGrid()
        hx.plotHex()
        hx.removeCell(0)
        hx.writeGrid()
        hx.plotHex()
        hx.removeCell(0)
        hx.writeGrid()
        hx.plotHex()
        hx.removeCell(len(hx.cells)-1)
        hx.writeGrid()
        hx.plotHex()
        hx.removeCell(len(hx.cells)-1)
        hx.writeGrid()
        hx.plotHex()

if __name__ == '__main__':
        main()
