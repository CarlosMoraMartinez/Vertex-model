import matplotlib.pyplot as plt 
import numpy as np
import scipy.spatial as spatial
from scipy.spatial import Voronoi, voronoi_plot_2d

import make_hexagonal_grid as mhg

def makeWing():
    argdict = {}
    argdict.setdefault('Size', 3)
    argdict.setdefault('Rows', 10)
    argdict.setdefault('Cols', 50)
    argdict.setdefault('Noise', 0.3)
    argdict.setdefault('Pull', 0)
    argdict.setdefault('Strecht', 1.1)
    argdict.setdefault('Rotate', 0)
    argdict.setdefault('StaticVertices', '')
    argdict.setdefault('Springs', '')
    argdict.setdefault('SpringLength', 0)
    argdict.setdefault('Hinge', -1)
    argdict.setdefault('Veins', '')
    argdict.setdefault('Outname', 'hex')
    # Make hexagonal grid
    hx = mhg.HexGrid(**argdict)
    fig, ax = hx.plotHex(False)
    # Draw border by hand 
    num_points, points, lines, fig, ax = getPointsBorder(fig, ax)
    removeCellsOutside

def getPointsBorder(fig, ax):
    points = []
    lines = []
    num_points = 0
    finished = False
    def onclick(event):
        nonlocal num_points
        nonlocal finished
        if(not finished):
            print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
                   ('double' if event.dblclick else 'single', event.button,
                   event.x, event.y, event.xdata, event.ydata))
            points.append([event.xdata, event.ydata])
            ax.scatter([event.xdata], [event.ydata], color = "blue")
            if(num_points > 0):
                lines.append([num_points - 1, num_points])
                ax.plot([points[num_points - 1][0], points[num_points][0]], [points[num_points - 1][1], points[num_points][1]], color = "black")
            fig.canvas.draw()
            num_points += 1
    def on_key(event):
        nonlocal finished
        nonlocal num_points
        print('you pressed', event.key, event.xdata, event.ydata)
        if(event.key in "mM" and not finished):
            finished = True
            nonlocal cid2
            fig.canvas.mpl_disconnect(cid2)
            lines.append([num_points - 1, 0])
            ax.plot([points[num_points - 1][0], points[0][0]], [points[num_points - 1][1], points[0][1]], color = "black")
            fig.canvas.draw()
        elif(finished):
            nonlocal cid1
            fig.canvas.mpl_disconnect(cid1)
    cid1 = fig.canvas.mpl_connect('key_press_event', on_key)
    cid2 = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    return (num_points, points, lines, fig, ax)


makeWing()
