

import sys
import numpy as np
import matplotlib.pyplot as plt


cell_ext = '.cells'
vert_ext = '.points'
cent_ext = '.cent'

size = 2.5
nrow = 20
ncol = 20
noise = 0.1

width = lambda size: np.sqrt(3) * size
height = lambda size: 2 * size


def getPoint(x, y, size, i):
    angle_deg = -1*60*(i + 1) + 30
    angle_rad = np.pi/180 * angle_deg
    return (x + size*np.cos(angle_rad), y + size*np.sin(angle_rad))

#Hexagon centers are separated by width horizontally and by height*3/4 vertically

def getCells(size, nrow, ncol):
    max_coord = height(size)*(ncol+1)
    centers = [ (j-1, i-1, width(size)*i + width(size)*0.5*(j%2) , height(size)*0.75*j)  for j in range(1, nrow + 1) for i in range(1, ncol + 1)]
    vertices = []
    cells = []
    vnum = 0
    for i, j, x, y in centers:
        cell = []
        for v in range(6):
            newv = -1
            if(j > 0):
                if(v == 2):
                    newv = cells[ncol*i + j - 1][0]
                elif(v == 3):
                    newv = cells[ncol*i + j - 1][5]
            if(i > 0):
                if(i%2 == 1):
                    if (v == 1):
                        newv = cells[ncol*(i-1) + j][3]
                    elif(v == 0):
                        newv = cells[ncol*(i-1) +j ][4]
                else:
                    if (v == 2):
                        newv = cells[ncol*(i-1) + j][4]
                    elif(v == 1):
                        newv = cells[ncol*(i-1) +j ][5]   
                    elif(v == 0 and j < ncol - 1):
                        newv = cells[ncol*(i-1) +j  + 1][4]
            if(newv == -1):
                newv = vnum
                a, b = getPoint(x, y, size, v)
                vertices.append([a, b, newv])
                vnum+=1
            cell.append(newv)
        cells.append(cell)
    return (centers, cells, vertices)



 

def addNoise(vertices, size, noise):
    for i in range(len(vertices)):
        vertices[i][0] = vertices[i][0] + np.random.uniform(-1*size*noise, size*noise)
        vertices[i][1] = vertices[i][1] + np.random.uniform(-1*size*noise, size*noise)
    return vertices



def plotHex(centers, vertices, cells, outname = 'hexagonal_grid'):
    f, ax = plt.subplots()
    ax.scatter([i[2] for i in centers], [i[3] for i in centers])
    for i in range(len(centers)):
        ax.annotate(i, [centers[i][2], centers[i][3]])
#    for i in range(len(vertices)):
#        ax.annotate(i, [vertices[i][0], vertices[i][1]])
    for c in cells:
        for i in range(6):
            plt.plot([vertices[c[i]][0], vertices[c[(i+1)%6]][0]] ,  [vertices[c[i]][1], vertices[c[(i+1)%6]][1]])
    plt.savefig( outname + '.png')
    plt.show()



def writeGrid(vertices, cells, centers, outname):

    f = open(outname + vert_ext,  'w')
    f.write(str(len(vertices)))
    f.write('\n')
    for v in vertices:
        f.write('\t'.join([str(i) for i in v]))
        f.write('\n')
    f.close()
    f = open(outname + cell_ext,  'w')
    f.write(str(len(cells)) + "\t" + str(6))
    f.write('\n')
    for c in cells:
        f.write('\t'.join([str(i) for i in c]))
        f.write('\t-999\n')
    f.close()
    f = open(outname + cent_ext,  'w')
    for c in centers:
        f.write('\t'.join([str(i) for i in v]))
        f.write('\n')
    f.close()



def main():
    if(sys.argv[2] != '-1'):
        s = float(sys.argv[2])
    else:
        s = size

    if(sys.argv[3] != '-1'):
        nr = int(sys.argv[3])
    else:
        nr = nrow

    if(sys.argv[4] != '-1'):
        nc = int(sys.argv[4])
    else:
        nc = ncol

    if(sys.argv[5] != '-1'):
        ran = float(sys.argv[5])
    else:
        ran = noise
    outname = sys.argv[1] + '_s'+str(s) + '_' + str(nr) + 'x' + str(nc) + '_n' + str(ran)
    print('size: ', s, '; num. rows: ', nr, '; num. cols: ', nc, '; noise: ', ran, '; output files: ', outname)
    centers, cells, vertices = getCells(s, nr, nc)
    #plotHex(centers, vertices, cells)
    if(ran > 0):
        vertices = addNoise(vertices, s, ran)
    writeGrid(vertices, cells, centers, outname)
    plotHex(centers, vertices, cells, outname)



if __name__ == '__main__':
    main()
