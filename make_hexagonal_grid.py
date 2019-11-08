import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt


cell_ext = '.cells'
vert_ext = '.points'
cent_ext = '.cent'
spring_ext = '.spr'

size = 2.5
nrow = 20
ncol = 20
noise = 0.1

bladetype = 0
hingetype = 1
veintype = 2
veinhinge = 3
wingcols = ['blue', 'green', 'black', 'gray']

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
                vertices.append([a, b, newv, 1]) # coord_x, coord_y, index, movable
                vnum+=1
            cell.append(newv)
        cells.append(cell)
    return (centers, cells, vertices)



 

def addNoise(vertices, size, noise):
    for i in range(len(vertices)):
        vertices[i][0] = vertices[i][0] + np.random.uniform(-1*size*noise, size*noise)
        vertices[i][1] = vertices[i][1] + np.random.uniform(-1*size*noise, size*noise)
    return vertices



def plotHex(centers, vertices, cells, springs, celltypes, outname = 'hexagonal_grid'):
    f, ax = plt.subplots()
    ax.scatter([i[2] for i in centers], [i[3] for i in centers], c = [wingcols[k] for k in celltypes])
    for i in range(len(centers)):
        ax.annotate(i, [centers[i][2], centers[i][3]])
#    for i in range(len(vertices)):
#        ax.annotate(i, [vertices[i][0], vertices[i][1]])
    for j, c in enumerate(cells):
        for i in range(6):
            plt.plot([vertices[c[i]][0], vertices[c[(i+1)%6]][0]] ,  [vertices[c[i]][1], vertices[c[(i+1)%6]][1]], c = wingcols[celltypes[j]])
    for c in springs:
        plt.plot([vertices[c[0]][0], vertices[c[1]][0]] ,  [vertices[c[0]][1], vertices[c[1]][1]], c = "red")
    for v in vertices:
        if(v[3] == 0):
           ax.scatter(v[0], v[1], c = "red")
    plt.savefig( outname + '.png')
    plt.show()



def writeGrid(vertices, cells, centers, springs, celltypes, outname):

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
    for i, c in enumerate(cells):
        f.write('\t'.join([str(i) for i in c]))
        f.write('\t-999\t' + str(celltypes[i]) + '\n')
    f.close()
    f = open(outname + cent_ext,  'w')
    for c in centers:
        f.write('\t'.join([str(i) for i in v]))
        f.write('\n')
    f.close()
    f = open(outname + spring_ext,  'w')
    f.write(str(len(springs)) + '\n')
    for c in springs:
        f.write('\t'.join([str(i) for i in c]))
        f.write('\n')
    f.close()


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
   

def addStatic(cells, vertices, static_vertices, spring_vertices, spring_length, nr, nc):
    left, top, bottom, right = trimStaticPoints(static_vertices)
    lefts, tops, bottoms, rights = trimStaticPoints(spring_vertices)
    #static vertices in cells (only in borders)
    for i in left:
        cnum  = cells[nc*i]
        vertices[cnum[2]][3] = 0
        vertices[cnum[3]][3] = 0
    for i in top:
        cnum  = cells[nc*(nr-1) + i]
        vertices[cnum[3]][3] = 0
        vertices[cnum[4]][3] = 0
        vertices[cnum[5]][3] = 0
    for i in bottom:
        cnum  = cells[i]
        vertices[cnum[0]][3] = 0
        vertices[cnum[1]][3] = 0
        vertices[cnum[2]][3] = 0
    for i in right:
        cnum  = cells[nc*i + nc - 1]
        vertices[cnum[0]][3] = 0
        vertices[cnum[5]][3] = 0
    #vertices joined to springs
    springs = []        
    for i in lefts:
        cnum  = cells[nc*i]
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
            x =  vertices[v_to_fix][0] - spring_length  
            y =  vertices[v_to_fix][1]
            ind = len(vertices)  
            vertices.append([x, y, ind, 0])
            springs.append((v_to_fix, ind))
    for i in tops:
        cnum  = cells[nc*(nr-1) + i]
        v_to_fix = cnum[4]
        x =  vertices[v_to_fix][0] 
        y =  vertices[v_to_fix][1] + spring_length  
        ind = len(vertices)
        vertices.append([x, y, ind, 0])
        springs.append((v_to_fix, ind))
    for i in bottoms:
        cnum  = cells[i]
        v_to_fix = cnum[1]
        x =  vertices[v_to_fix][0] 
        y =  vertices[v_to_fix][1] - spring_length  
        ind = len(vertices)
        vertices.append([x, y, ind, 0])
        springs.append((v_to_fix, ind))
    for i in rights:
        cnum  = cells[nc*i + nc - 1]
        if(i == nr-1 and i%2 == 1):
            addspringto = [5]
        elif(i%2 == 0):
            addspringto = [0, 5]
        else:
           addspringto = []
        for j in addspringto:
            v_to_fix = cnum[j]
            x =  vertices[v_to_fix][0] + spring_length  
            y =  vertices[v_to_fix][1]
            ind = len(vertices)  
            vertices.append([x, y, ind, 0])
            springs.append((v_to_fix, ind))
    return (vertices, springs)



def addCellType(cells, hingelimit, veinpos, nr, nc):
    if(hingelimit is None):
        hingelimit = -1
    if(veinpos != ''):
        veinpos = [int(i) for i in str.split(veinpos, ',')]
    else:
        veimpos = [nr + 1]
    types = []
    for i, cell in enumerate(cells):
        col = i%nc
        row = i//nc
        if(row in veinpos):
            if(col > hingelimit):
                types.append(veintype)
            else:
                types.append(veinhinge)
        elif(col <= hingelimit):
            types.append(hingetype)
        else:
            types.append(bladetype)
    return types



parser = argparse.ArgumentParser(description='Hexagonal grid arguments.')
parser.add_argument('-o', '--Outname', metavar='outname', type=str, default = "hexgrid", 
                    help='Identifier. Used as prefix of all output files. ')
parser.add_argument('-s', '--Size', metavar='size', type=float, default = size,
                    help=' Hexagon size: distance between center and vertices. ')
parser.add_argument('-r', '--Rows', metavar='rows', type=int, default=nrow,
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


parser.add_argument('-g', '--Hinge', metavar='Hinge_pos', type=int, default = None,
                    help='Hexagons in columns in range [0:g] will be considered hinge (value 1).')
parser.add_argument('-v', '--Veins', metavar='Veins_pos', type=str, default = '',
                    help='Hexagons in rows specified (as integers) will be considered veins ( value 2).')
def main():

    args = parser.parse_args()

    s = args.Size
    nr = args.Rows
    nc = args.Cols
    ran = args.Noise

    static_vertices = args.StaticVertices
    spring_vertices = args.Springs
    spring_length = args.SpringLength

    hingelimit = args.Hinge
    veinpos = args.Veins
    outname = args.Outname + '_s'+str(s) + '_' + str(nr) + 'x' + str(nc) + '_n' + str(ran)
    print('size: ', s, '; num. rows: ', nr, '; num. cols: ', nc, '; noise: ', ran, '; output files: ', outname)

    centers, cells, vertices = getCells(s, nr, nc)
  
    if(ran > 0):
        vertices = addNoise(vertices, s, ran)
    if(static_vertices != ''):
        vertices, springs = addStatic(cells, vertices, static_vertices, spring_vertices, spring_length, nr, nc)
    else:
        springs = []
    if(hingelimit is not None or veinpos != ''):
       celltypes = addCellType(cells, hingelimit, veinpos, nr, nc)
    else:
        celltypes = [bladetype for c in cells]

    writeGrid(vertices, cells, centers, springs, celltypes, outname)
    plotHex(centers, vertices, cells, springs, celltypes, outname)



if __name__ == '__main__':
    main()
