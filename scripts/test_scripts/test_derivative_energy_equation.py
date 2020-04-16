
import numpy as np

pol1=[[1000, 50],[900, 40],[800, 39],[700, 55],[820, 65],[980, 60]] #area=5150
pol2=[[1000, 100*50],[900, 100*40],[800, 100*39],[700, 100*55],[820, 100*65],[980, 100*60]] #area=5150
mx, my = (np.mean([p[0] for p in pol1]), np.mean([p[1] for p in pol1]))
vx, vy = (np.std([p[0] for p in pol1]), np.std([p[1] for p in pol1]))
pol3 = [[100*(p[0] - mx)/vx, 100*(p[1] - my)/vy ] for p in pol1]
pol4 = [[(p[0] - mx), (p[1] - my) ] for p in pol1]


rect = [[0,0],[0,20],[3,20],[3,0]]

a0, tension, contract, h = (5149.0, 1.0, 1.0, 1.0)

def centroid(pol):
    return [i for i in map(np.mean, zip(*pol))]
def area(pol):
    s = len(pol)
    a = 0
    for i in range(s):
        a += pol[i][0]*(pol[(i+1)%s][1] - pol[(s+i-1)%s][1])
    return 0.5*abs(a)

def dist(p1, p2):
    return np.sqrt( (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 )

def perim(pol):
    s = len(pol)
    p = 0
    for i in range(s):
        p += dist(pol[i], pol[(i+1)%s])
    return p


def term1(pol, p=0, coord=0):
    coord2 = 1 if coord==0 else 0
    a = area(pol)
    c = centroid(pol)
    s = len(pol)
    y_after = pol[(p + 1)%s][coord2]
    y_before = pol[(s + p - 1)%s][coord2]
    both = (y_after, y_before)
    if(pol[p][coord] > c[coord]):
        y_after, y_before = (max(both), min(both))
    else:
        y_after, y_before = (min(both), max(both))
    t = (2*a/a0 - 2)*(y_after - y_before)/(2*a0)
    return t


def term2(pol, p=0, coord=0):
    s = len(pol)
    ta = (pol[p][coord] - pol[(p + 1)%s][coord])/dist(pol[p], pol[(p + 1)%s]) #abs
    tb = (pol[p][coord] - pol[(s + p - 1)%s][coord])/dist(pol[p], pol[(s + p - 1)%s]) #abs    
    return tension/a0 * (ta + tb)

def term3(pol, p=0, coord=0):
    s = len(pol)
    per = perim(pol)
    ta = (pol[p][coord] - pol[(p + 1)%s][coord])/dist(pol[p], pol[(p + 1)%s])
    tb = (pol[p][coord] - pol[(s + p - 1)%s][coord])/dist(pol[p], pol[(s + p - 1)%s])
    return (ta + tb) * per*contract/a0

def dEdx(pol, p=0, coord=0):
    return term1(pol, p, coord) + term2(pol, p, coord) + term3(pol, p, coord)


def force(pol, p=0, coord=0):
    return -1*dEdx(pol, p, coord)

def derivativePolygon(pol):
    der = []
    for i in range(len(pol)):
        der.append([])
        for j in range(2):
            d = force(pol, i, j)*h
            der[i].append(d)
    return der

def simulate(pol, steps=1000, get_history=True):
    pol = [p.copy() for p in pol]
    history = [[p.copy() for p in pol]]
    history_der = []
    for t in range(steps):
        der = derivativePolygon(pol)
        #print(der)
        for i in range(len(pol)):
            for j in range(2):
                pol[i][j] += der[i][j]
        if(get_history):
            history.append([p.copy() for p in pol])
            history_der.append([p.copy() for p in der])
    if(get_history):
        der = derivativePolygon(pol)
        history_der.append([p.copy() for p in der])
        return (history, history_der)
    else:
        return pol

def plotPolygon(pol, col="red"):
    s = len(pol)
    for p in range(s):
        plt.plot([pol[p][0], pol[(p+1)%s][0]], [pol[p][1], pol[(p+1)%s][1]], color=col)
        
def plotArrow(point, der, cols=["red"]):
    for i, p in enumerate(pol):
        plt.arrow(p[0], p[1], der[i][0], der[i][1], color = cols[i % len(cols)],  head_width=0.1)

def plotAll(pols, der=None):
    print("First is red, last is black")
    cols = ["red"] +["blue" for i in range(1, len(pols)-1)] + ["black"]
    for p, c in zip(pols, cols):
        plotPolygon(p, c)
    if(der is not None):
        plotArrow(pols[0], der[0], ["red"])
        plotArrow(pols[-1], der[-1], ["black"])
    plt.axis("equal")
    plt.show()

def plotAllArrows(pols, der):
    cols=["red", "blue", "green", "black", "yellow", "pink", "purple", "gray", "orange", "blue", "red"]
    for p, d in zip(pols, der):
        plt.scatter([i[0] for i in p],[i[1] for i in p], color=[cols[i%len(cols)] for i, c in enumerate(p)])
        plotArrow(p, d, cols)
    plt.axis("equal")
    plt.show()

def plotDer(pols):
    points = [p for p in zip(*pols)]
    for i, p in enumerate(points):
        x, y = zip(*p)
        plt.plot(range(len(x)), x, label="x"+str(i) )
        plt.plot(range(len(y)), y, label="y"+str(i) )         
    plt.legend()
    plt.show()

#Eq. area
tension=0;contract=0;a0=area(pol1)*10;h=0.01;his, der=simulate(pol1, 100);plotAll([his[0],his[-1]], der);plotDer(der)

#tension
a0=area(pol3);tension=-100;contract=0;h=10;his, der=simulate(pol3, 1000);plotAll([his[0],his[-1]]);plotDer(der)
a0=area(pol3);tension=1000;contract=0;h=1;his, der=simulate(pol3, 1000);plotAll([his[0],his[-1]]);plotDer(der)

#rectangle to square
a0=area(rect);tension=0;contract=0.1;h=0.1;his, der=simulate(rect, 10000);plotAll([his[0],his[-1]]);plotDer(der)
dist(his[0][0], his[0][1])/dist(his[0][1], his[0][2])
dist(his[-1][0], his[0][1])/dist(his[0][1], his[-1][2])

       

