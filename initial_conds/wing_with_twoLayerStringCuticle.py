import os
import argparse
import copy

import numpy as np
import matplotlib.pyplot as plt

import make_hexagonal_grid as mhg

NUM_CELLS_PER_VERTEX = 3
class WingWith2StringLayers:
    dist = staticmethod(lambda b, c:np.sqrt(np.power(b[0] - c[0], 2) + np.power(b[1] - c[1], 2)))
    def __init__(self, arguments): 
        self.original_wing_name = arguments.Inputname
        argdict = {}
        argdict.setdefault('Size', 3)
        argdict.setdefault('Rows', int(1.1*100)) #*self.imratio
        argdict.setdefault('Cols', 210)
        argdict.setdefault('Noise', 0.3)
        argdict.setdefault('Pull', 0)
        argdict.setdefault('Strecht', 1.1)
        argdict.setdefault('Rotate', 0)
        argdict.setdefault('StaticVertices', '')
        argdict.setdefault('Springs', '')
        argdict.setdefault('SpringLength', 0)
        argdict.setdefault('Hinge', -1)
        argdict.setdefault('Veins', '')
        argdict.setdefault('Outname', arguments.Outname)
        argdict.setdefault('Read', self.original_wing_name)
        self.argdict = argdict
        print("Reading base wing")
        self.wing = mhg.HexGrid(**argdict) #Read wing
        #print("Base wing read. Generating grid for cuticle layer.")
        #argdict["Read"] = ""
        #self.cuticle_layer = mhg.HexGrid(**argdict) #Create hexagonal grid for the cuticle
        #Center Grid position
        self.new_vertices = []
        self.new_strings = []
        self.vinstrings = np.unique([i for j in self.strings for i in j]).tolist()
        self.vincells = np.unique([i for j in self.cells for i in j]).tolist()
        self.onlyinstrings = [i for i in self.vinstrings if i not in self.vincells and i > 0]
        self.setCenterPosition()
    @property
    def vertices(self):
        return self.wing.vertices
    @property
    def cells(self):
        return self.wing.cells
    @property
    def celltypes(self):
        return self.wing.celltypes
    @property
    def springs(self):
        return self.wing.springs
    @property
    def spring_vertices(self):
        return self.wing.spring_vertices
    @property
    def strings(self):
        return self.wing.stringEdges
    def setCenterPosition(self):
        x_pos, y_pos, _, _ = [i for i in zip(*[self.vertices[j] for j in self.onlyinstrings])] 
        self.center_x = 0.5*(np.min(x_pos) + np.max(x_pos))
        self.center_y = 0.5*(np.min(y_pos) + np.max(y_pos))
    def makeCuticle(self, size=5, nlayers=1):
        print("a")
        self.wing.cuticle_type = 1 if nlayers == 1 else 2
        connections = [self.strings[i][0] == self.strings[i-1][1] and self.strings[i][1] == self.strings[i+1][0] for i in range(1, len(self.strings) -1) ]
        assert(all(connections), "Error: assumption that strings are ordered is not met. You will have to fix the wing cuticle or change this function")
        #Just ordering them can work
        #Create new vertices in new positions
        assert(all([x[2] == i for i, x in enumerate(self.vertices)]), "Danger!! Number of vertices is not correct")
        prev_layer = self.strings
        next_layer = []
        for layer in range(nlayers):
            print("Starting layer ", layer)
            vind = len(self.vertices)
            previous = -1
            #print("b")
            for i in range(1, len(prev_layer)):
                current_old = prev_layer[i][0]
                self.new_strings.append([current_old, vind])
                #print("c")
                newx, newy = self.getNewPosition(self.vertices[current_old][0], self.vertices[current_old][1], size)
                if(i == 1 or i == (len(prev_layer)-1)):
                    static = 0
                else:
                    static = 1
                self.new_vertices.append([newx, newy, vind, static])
                #print("d")
                if(previous != -1):
                    #print("e")
                    self.new_strings.append([previous, vind])
                    self.new_strings.append([previous, current_old])
                    self.new_strings.append([current_old - 1, vind])
                    next_layer.append([previous, vind])
                    #print("f")
                previous = vind
                vind += 1
                #print("g")
                #Add connections
            self.wing.stringEdges.extend(self.new_strings)
            self.wing.vertices.extend(self.new_vertices)
            self.new_vertices = []
            self.new_strings = []
            prev_layer = next_layer
            next_layer = []
    def scatterNewVertices(self):
        for v in self.new_vertices:
            plt.scatter(v[0], v[1], color="RED")
        for v2 in self.onlyinstrings:
            plt.scatter(self.vertices[v2][0], self.vertices[v2][1], color = "BLUE")
    def getNewPosition(self, x, y, size): 
        dist = np.sqrt(np.power(x-self.center_x, 2) + np.power(y-self.center_y, 2))
        vec_x = x + size*(x-self.center_x)/dist
        vec_y = y + size*(y-self.center_y)/dist
        return (vec_x, vec_y)
    def plotNewWing(self):
        self.wing.plotHex2()
        plt.show()
    def writeNewGrid(self):
        self.wing.writeGrid()
    def getStringLength(self):
        lengths = []
        for a, b in self.strings:
            l = np.sqrt(np.power(self.vertices[a][0] - self.vertices[b][0], 2) + np.power(self.vertices[a][1] - self.vertices[b][1], 2))
            lengths.append(l)
        return lengths
    def removeStringsWithLengthHigherThan(self, maxlen=50):
        lengths = self.getStringLength()
        for i in range(len(lengths) - 1, -1, -1):
            if(lengths[i] >= maxlen):
                x = self.wing.stringEdges.pop(i)
                print("Removed string %d [%d, %d], with length = %f"%(i, x[0], x[1], lengths[i]))


    







inputname = "wing.png"


parser = argparse.ArgumentParser(description='Wing Drawer arguments.')
parser.add_argument('-o', '--Outname', metavar='outname', type=str, default = "hexgrid", 
                                        help='Identifier. Used as prefix of all output files. ')
parser.add_argument('-i', '--Inputname', metavar='inputname', type=str, default = "", 
                                        help='Identifier. Used as prefix to read files. ')
parser.add_argument('-w', '--Inputimage', metavar='inputimage', type=str, default = inputname, 
                                        help='Use this image to draw on it. ')
parser.add_argument('-s', '--CuticleThickness', metavar='thickness', type=int, default = 2, 
                                        help='Number of cell layers of cuticle')
parser.add_argument('-x', '--MinX', metavar='cuticle_start', type=int, default = 2, 
                                        help='Proportion of the wing x axis at which the cuticle starts')                                        
#makeWingVoronoi()
def main():
    args = parser.parse_args()
    #args = {"Inputname":"etournay1_strings10", "Outname":"cuttest4layers"}
    #obj = type("args", (object,), args)
    # ww = WingWith2StringLayers(obj)
    cuticle_width = args.CuticleThickness
    print("WARNING: Make sure that your starting wing does not have vertices that only touch springs or string-like cuticle.")
    ww = WingWith2StringLayers(args)
    ww.makeCuticle(cuticle_width)
    ww.plotNewWing()
    ww.writeNewGrid()
    ww.removeStringsWithLengthHigherThan(50) #Modify length according to your wing and cuticle string length

if __name__ == '__main__':
    main()

