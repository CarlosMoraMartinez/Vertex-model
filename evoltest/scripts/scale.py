
import sys
import os

EXT_EDGES = '.edges'
EXT_VERTICES = '.points'
EXT_TARGET = '.border'
EXT_CELL = '.cells'
EXT_CELLTAB = '.celltab'
EXT_SPR = '.spr'

copyfiles = [EXT_EDGES, EXT_CELL, EXT_CELLTAB, EXT_SPR]

def readPoints(fbase):
    with open(fbase + EXT_VERTICES, 'r') as pf:
        n = int(pf.readline())
        p = [[float(i) if '.' in i else int(i) for i in l.strip('\n').split('\t')] for l in pf]
    return (n, p)

def writePoints(fbase, n, p):
    with open(fbase + EXT_VERTICES, 'w') as pf:
        pf.write(str(n) + '\n')
        for pp in p:
            pf.write('\t'.join([str(x) for x in pp]) + '\n')

def scalePoints(n, p, xscale, yscale, xadd, yadd):
    for i in range(n):
        p[i][0] = p[i][0]*xscale + xadd
        p[i][1] = p[i][1]*yscale + yadd
    return p

def copyOtherFiles(i, o):
    for fn in copyfiles:
        if(os.path.exists(i + fn)):
            fi = open(i + fn, 'r')
            fo = open(o + fn, 'w')
            for l in fi:
                fo.write(l)
            fi.close() 
            fo.close()

def main():
    inname = sys.argv[1]
    xscale = float(sys.argv[2])
    yscale = float(sys.argv[3])
    xadd = float(sys.argv[4])
    yadd = float(sys.argv[5])
    outname = sys.argv[6]

    print("starting")
    n, p = readPoints(inname)
    p2 = scalePoints(n, p, xscale, yscale, xadd, yadd)
    writePoints(outname, n, p2)
    print("scaled")
    copyOtherFiles(inname, outname)
    print("Other files copied")


if(__name__ == "__main__"):
    main()


