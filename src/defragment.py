#import sys
import argparse


EMPTY= '-999'

def readPointsFile(name):
    pointsFile = open(name + ".points", "r")
    numPoints = int(pointsFile.readline())

    # This is the list of points
    pointsList = []#np.zeros((numPoints,2))
    
    for row in pointsFile:
        ll = row.split()
        pointsList.append( [float(ll[0]), float(ll[1]), int(ll[2]),  int(ll[3])])
    pointsFile.close()
    print(numPoints, " points")
    return (numPoints, pointsList)

def readCellsFile(name):
    cellsFile = open(name + ".cells", "r")
    numCells, numVerticesCell = [int(each) for each in cellsFile.readline().split()]

    polygonList = []
    celltypes = []
    for row in cellsFile:
        polygonIndex = [int(indPoint) for indPoint in row.split() if int(indPoint) >= 0]
        celltypes.append(int(polygonIndex.pop())) #CAREFUL: if cells file does not have type, will produce error
        #polygonCoords = [[pointsList[coord][0], pointsList[coord][1]] for coord in polygonIndex]
        polygonList.append(polygonIndex)
    cellsFile.close()
    print(numCells, " cells")
    return (numCells, numVerticesCell, polygonList, celltypes)

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
        return (0, [])
    print(numsprings, " springs")
    return (numsprings, sprList)


def defrag(points, cells, springs):
    for i, p in enumerate(points):
        if(i == p[2]):
            continue
        print("Converting %d to %d"%(p[2], i))
        for c, cell in enumerate(cells):
            for j, n in enumerate(cell):
                if(n == p[2]):
                    cells[c][j] = i
                    print("In cell: %d in cell, %d in vert, converted to %d"%(n, p[2], i))
        for s, spr in springs:
            for j, n in enumerate(spr):
                if(j < 2 and n == p[2]):
                    spr[s][j] = i #j[2] is likely spring type
        p[2] = i
        print("Converted: %d is %d"%(p[2], i))
    return points, cells, springs
        


def writePoints(points, name):
    with open(name + '_df.points', 'w') as f:
        f.write(str(len(points)) + '\n')
        f.write('\n'.join( ['\t'.join([str(x) for x in p]) for p in points] ))  

def writeCells(cells, celltypes, name, lencells=30):
    with open(name + '_df.cells', 'w') as f:
        f.write("%d\t%d\n"%(len(cells), lencells))      
        f.write('\n'.join( ['\t'.join([str(c[i]) if i < len(c) else EMPTY for i in range(lencells)] + [str(t)])  for c, t in zip(cells, celltypes)]))



def writeSprings(spr, name):
    with open(name + '_df.spr', 'w') as f:
        f.write(str(len(spr)) + '\n')
        f.write('\n'.join( ['\t'.join([str(x) for x in p]) for p in spr] ))  

parser = argparse.ArgumentParser(description='This script makes the indices of vertices in a wing consecutive (e.g., 0, 1, 4, 6 -> 0, 1, 2, 3) and changes cells and springs accordingly. It does not change cuticle vertices. Useful when using the result of a simulation as the initial condition for a new simulation, in case there are dead cells or vertices. The output files will have the "_df" suffix. Call this script from inside the directory where the .cells, .points and .spr files are.')
parser.add_argument('-i', '--Inputname', metavar='inputname', type=str, default = "", 
                                        help='Identifier. Used as prefix to read files. ')

def main():
    args = parser.parse_args()
    name = args.Inputname 
    #name = sys.argv[1]
    print("Reading wing " + name)
    numpoints, points = readPointsFile(name)
    numCells, numVerticesCell, cells, celltypes = readCellsFile(name)
    numsprings, spr = readSprings(name)
    print("Files read.")
    points2, cells2, spr2 = defrag(points, cells, spr)
    writePoints(points2, name)
    writeCells(cells2, celltypes, name, numVerticesCell)
    writeSprings(spr2, name)


if __name__ == "__main__":
    main()

