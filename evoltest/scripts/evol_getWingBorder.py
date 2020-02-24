import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import shapely
from shapely.geometry import Polygon



EXT_EDGES = '.edges'
EXT_VERTICES = '.points'
EXT_TARGET = '.border'
BORDER_TYPE = 4


class WingFitnessEvaluator:
    target_table=None
    polt=None
    def __init__(self):
        pass
    @classmethod
    def setTargetShape(cls, ftarget):
        """
        Assumes that points in input are ordered
        """
        cls.target_table = pd.read_csv(ftarget + EXT_TARGET, sep="\t")
        cls.polt = Polygon(cls.target_table.to_numpy().tolist())
    @classmethod    
    def getFitness(cls, fbase):    
        polw = Polygon(cls.getLinearPath(fbase))
        return cls.errorFunc(cls.polt, polw)
    @classmethod   
    def getTargetAsList(cls):
        return cls.target_table.to_numpy().tolist()
    @classmethod
    def getLinearPath(cls, fbase):
        """
        Input: file name without extension. There must be at least a .edges and a .points file with this name
        Returns a list of consecutive points (as [x, y] lists) representing the contour of the input shape
        """
        e = pd.read_csv(fbase + EXT_EDGES, sep="\t")
        #points file has a different format
        with open(fbase + EXT_VERTICES, 'r') as pf:
            n = int(pf.readline())
            p = [[float(i) if '.' in i else int(i) for i in l.strip('\n').split('\t')] for l in pf]
        p = pd.DataFrame( dict(zip(['x', 'y', 'ind', 'movable'], zip(*p))) ).set_index("ind", drop=False)    
        eborder = e[e['type'] == BORDER_TYPE]
        vertices = [[int(i) for i in x.split(',') if i != ''] for x in eborder['vertices']]
        return cls.getPath(vertices, p)    
    #easyPlot(wingpath, t.to_numpy().tolist())
    @classmethod
    def getPath(cls, edgev, p):
        border = []
        start = edgev.pop()
        border.append(p[['x', 'y']].loc[start[0]])
        border.append(p[['x', 'y']].loc[start[1]])   
        current = start[1]
        while(edgev):
            for i in range(len(edgev)):
                if(current == edgev[i][0]): 
                    e = edgev.pop(i)
                    border.append(p[['x', 'y']].loc[e[1]]) 
                    current = e[1]                 
                    break
                elif(current == edgev[i][1]): 
                    e = edgev.pop(i)
                    border.append(p[['x', 'y']].loc[e[0]]) 
                    current = e[0]                 
                    break
            else:
                print("Error: border not continuous")
                return []
        return border
    @staticmethod
    def easyPlot(path1, path2):
        def plotpath(path, color):
            for i in range(1, len(path)):
                x, y = zip(*[path[i - 1], path[i]])
                plt.plot(x, y, color = color, alpha=0.5)
            x, y = zip(*[path[- 1], path[0]])
            plt.plot(x, y, color = color, alpha=0.5)    
        for p, c in zip([path1, path2], ["red", "blue"]):
            plotpath(p, c)
        plt.show()
    @staticmethod  
    def getPos(v, p):
        return p[['x', 'y']].loc[v].tolist()
    @staticmethod
    def errorFunc(p1, p2):
        intersect = p1.intersection(p2)
        union = p1.union(p2)
        return intersect.area/union.area



def main():
    import sys
    ftarget = sys.argv[1]
    feval = sys.argv[2]
    WingFitnessEvaluator.setTargetShape(ftarget)
    err = WingFitnessEvaluator.getFitness(feval)
    WingFitnessEvaluator.easyPlot(WingFitnessEvaluator.getTargetAsList(), WingFitnessEvaluator.getLinearPath(feval))
    

if(__name__ == "__main__"):
    main()

















