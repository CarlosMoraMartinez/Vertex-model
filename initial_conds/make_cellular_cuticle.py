import os
import argparse
import copy

import numpy as np
import matplotlib.pyplot as plt

import make_hexagonal_grid as mhg

NUM_CELLS_PER_VERTEX = 3
class WingWithCuticle:
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
        print("Base wing read. Generating grid for cuticle layer.")
        argdict["Read"] = ""
        self.cuticle_layer = mhg.HexGrid(**argdict) #Create hexagonal grid for the cuticle
        #Center Grid position
        self.__centerCuticleGridPosition()
        self.borderpoints = []
        self.closest = [] #Closest will hold, for each point in the wing border, which point in the cuticle grid is the closest
        self.grid_cells_used = []
        self.cuticle_cells = []
        self.cuticle_vertices = []
        self.cuticle_vertex_conversion = dict()
        self.cuticle_celltypes = []
        self.cuticle_springs = []
        self.initial_new_ind = 0
        self.new_wing = None
        self.vertex_use_count = []
    def __centerCuticleGridPosition(self):
        x_pos, y_pos, _, _ = [i for i in zip(*self.wing.vertices)] 
        self.center_x = 0.5*(np.min(x_pos) + np.max(x_pos))
        self.center_y = 0.5*(np.min(y_pos) + np.max(y_pos))
        x_pos, y_pos, _, _ = [i for i in zip(*self.cuticle_layer.vertices)]
        old_center_x = 0.5*(np.min(x_pos) + np.max(x_pos))
        old_center_y = 0.5*(np.min(y_pos) + np.max(y_pos))
        diff_x = old_center_x - self.center_x
        diff_y = old_center_y - self.center_y
        for i in range(len(self.cuticle_layer.vertices)):
            self.cuticle_layer.vertices[i][0] -= diff_x
            self.cuticle_layer.vertices[i][1] -= diff_y
    def makeCuticle(self, size):
        import itertools
        print("Making cuticle")
        #Remove string cuticle if exists
        self.initial_new_ind = np.max([i for c in self.wing.cells for i in c]) + 1
        self.wing.stringEdges = []
        del self.wing.vertices[self.initial_new_ind:]
        #Calculate border
        print("Calculating border in original wing")
        #Conditional to get only vertices that can move freely
        self.borderpoints = [self.wing.vertices[i] for i in self.wing.getBorderPoints() if self.wing.vertices[i][3] == 1]
        print("Getting points in the new grid closest to border")
        self.closest = [self.cuticle_layer.vertices[np.argmin([self.dist(i, c) for c in self.cuticle_layer.vertices])][2] for i in self.borderpoints]
        cells_used, cells_added_per_vertex = ([False for i in range(len(self.cuticle_layer.cells))], [0 for i in range(len(self.cuticle_layer.vertices))])
        vertices = [i for i in self.closest]
        #Get different layers in different iterations
        for i in range(size):
            cells_used, cells_added_per_vertex = self.findCellsForAllVertices(vertices, self.cuticle_layer.cells, cells_used, cells_added_per_vertex)
            #get new vertices, i.e., include the ones in new cells
            vertices = [self.cuticle_layer[i] for i in range(len(self.cuticle_layer.cells)) if  cells_used[i]]
            vertices = np.unique(list(itertools.chain.from_iterable(vertices))).tolist()
            #cells_added_per_vertex = cells_added_per_vertex + [0 for i in range(len(vertices) - len(cells_added_per_vertex))]
            print("Added cuticle vertices, layer " + str(i))
            #These lines are unnecessary; just to show which cells in the grid were selected
            for j, selected in enumerate(cells_used):
                if(selected):
                    self.cuticle_layer.celltypes[j] = mhg.cuticletype
            self.cuticle_layer.plotHex2()
            plt.show()
        self.grid_cells_used = [i for i in range(len(cells_used)) if cells_used[i]]
        self.vertex_use_count = cells_added_per_vertex
        new_v_ind = self.initial_new_ind
        for i in self.grid_cells_used:
            oldcell = self.cuticle_layer.cells[i]
            newcell = []
            for vio, v in enumerate(oldcell):
                #if(cells_added_per_vertex[v] == 1 and not v in self.closest):
                #    continue
                newvname = self.cuticle_vertex_conversion.setdefault(v, new_v_ind)
                if newvname == new_v_ind:
                    new_v_ind += 1 
                    self.cuticle_vertices.append(copy.copy(self.cuticle_layer.vertices[v]))
                    self.cuticle_vertices[-1][2] = newvname
                newcell.append(newvname)
            self.cuticle_cells.append(newcell)
        #Set springs and celltypes
        self.cuticle_celltypes = [mhg.cuticletype]*len(self.cuticle_cells)
        self.cuticle_springs = [[self.cuticle_vertex_conversion[self.closest[i]], self.borderpoints[i][2]] for i in range(len(self.closest))]
        #Finally, create new wing from the old one
        self.new_wing = copy.deepcopy(self.wing)
        #self.new_wing.vertices = [self.new_wing.vertices[i] for i in range(self.initial_new_ind)] + self.cuticle_vertices #Deleted strings before
        self.new_wing.vertices = self.new_wing.vertices + self.cuticle_vertices
        self.new_wing.cells = self.new_wing.cells + self.cuticle_cells
        self.new_wing.celltypes = self.new_wing.celltypes + self.cuticle_celltypes
        self.new_wing.vnum = len(self.new_wing.vertices)
        self.new_wing.springs = self.cuticle_springs
        self.plotNewWing()
        print("Printing provisional grid")
        self.writeNewGrid()
        print("Removing border extra vertices (touching only 1 cell)")
        #self.removeExtraVertices()
    def plotNewWing(self):
        self.new_wing.plotHex2()
        plt.show()
    @staticmethod
    def findCellsForAllVertices(vertices, cells, cells_used, cells_added_per_vertex):
        for v_ind, v in enumerate(vertices):
            if(cells_added_per_vertex[v] >= NUM_CELLS_PER_VERTEX): #Should not happen
                continue
            for c_ind in range(len(cells)):
                if(cells_used[c_ind]):
                    continue
                if(vertices[v_ind] in cells[c_ind]):
                    cells_used[c_ind] = True
                    cells_added_per_vertex[v] += 1
        return (cells_used, cells_added_per_vertex)
    def removeExtraVertices(self):        
        """This method removes vertices in the border of the cuticle that only have 2 neighbours."""
        self.new_wing.setVertexCellIndex()
        vcount = [len(v[4]) for v in self.new_wing.vertices]
        #vcount = [0 for i in range(len(self.new_wing.vertices))]
        #for c in self.cuticle_cells:
        #    for v in c:
        #        vcount[v] += 1
        vrem = [i for i, v in enumerate(vcount) if v == 1 and i >= self.initial_new_ind]
        vrem.sort()
        while(vrem):
            id = vrem.pop()
            print("Removing %d"%(id))
            vrem = self.new_wing.removeVertex(id, vrem)
    def writeNewGrid(self):
        self.new_wing.writeGrid()
    

inputname = "wing.png"


parser = argparse.ArgumentParser(description='This script adds a cellular cuticle to a wing initial condition.')
parser.add_argument('-o', '--Outname', metavar='outname', type=str, default = "hexgrid", 
                                        help='Identifier. Used as prefix of all output files. ')
parser.add_argument('-i', '--Inputname', metavar='inputname', type=str, default = "", 
                                        help='Identifier. Used as prefix to read files. ')
parser.add_argument('-w', '--Inputimage', metavar='inputimage', type=str, default = inputname, 
                                        help='Use this image to draw on it. ')
parser.add_argument('-s', '--CuticleThickness', metavar='thickness', type=int, default = 2, 
                                        help='Number of cell layers of cuticle')
parser.add_argument('-x', '--MinX', metavar='cuticle_start', type=int, default = 2, 
                                        help='Proportion of the wing x axis at which the cuticle starts (not implemented)')                                        
#makeWingVoronoi()
def main():
    args = parser.parse_args()
    #args = {"Inputname":"etournay1_strings10", "Outname":"cuttest"}
    #obj = type("args", (object,), args)
    # ww = WingWithCuticle(obj)
    cuticle_width = args.CuticleThickness
    print("WARNING: Make sure that your starting wing does not have vertices that only touch springs or string-like cuticle.")
    ww = WingWithCuticle(args)
    ww.makeCuticle(cuticle_width)
    ww.removeExtraVertices()
    ww.plotNewWing()
    ww.writeNewGrid()

if __name__ == '__main__':
    main()

