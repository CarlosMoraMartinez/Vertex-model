import sys

import numpy as np
import cv2
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, MultiPoint, Point
from descartes.patch import PolygonPatch

# Necessary steps to show the images
# Arguments: string with the name of the window and the final image object
def plotImage(windowName, imObj):
    cv2.imshow(windowName, imObj)
    cv2.waitKey(0)
    cv2.destroyAllWindows()
    return 0

def createEmptyImage(templateObj):
    emptyImage = np.zeros(templateObj.shape, np.uint8)
    return emptyImage

def plot_grid(plot_pos, grid):

    # In order to plot a MultiPolygon object, I need to iterate over each oplygon
    fig = plt.figure(1, figsize=(5,5), dpi=90)
    ax = fig.add_subplot(plot_pos)

    for element in grid:
            polygon = Polygon(element)
            plot_coords(ax, polygon.exterior)
            patch = PolygonPatch(polygon, facecolor=color_isvalid(polygon), edgecolor=color_isvalid(polygon, valid=BLUE), alpha=0.5, zorder=2)
            ax.add_patch(patch)

# import the image
im = cv2.imread(sys.argv[1])

# Pre process the image with some filters, in order to minimise the details
imgray = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)
ret, thresh = cv2.threshold(imgray, 127, 255, 0)
#plotImage("threshold", thresh)

listContours, h = cv2.findContours(thresh,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)

# Get the centroids for each contour
centroid_x = []
centroid_y = []
area_list = []
cellObjs = []
for eachContour in listContours:

    # Find the moments of the contours, whatever this means
    moment_list = cv2.moments(eachContour)

    # By using the moments, calculate the centroids and the area of each contour
    # Open polygons are open and then return area zero, I exclude these to avoid dividing by 0
    if(moment_list['m00'] != 0.0):
        centroid_x.append(moment_list['m10']/moment_list['m00'])
        centroid_y.append(moment_list['m01']/moment_list['m00'])
        area_list.append(cv2.contourArea(eachContour))

    # Fit the contours to a polygon, by the way, this function return only the vertices!!!!
    contourPerim = cv2.arcLength(eachContour,True)
    polyFittingPrecision = 0.1*cv2.arcLength(eachContour,True)
    polyFittingVertices = cv2.approxPolyDP(eachContour,polyFittingPrecision,True)
    polygonCoords = [eachCoord[0] for eachCoord in polyFittingVertices if len(polyFittingVertices) > 3]
    #for eachCoord in polyFittingVertices:
    #    if len(polyFittingVertices) > 2:
    #        print eachCoord[0]
    #        polygonCoords = eachCoord[0]

    #    if len(polyFittingVertices)> 2:
    #        Polygon(polygonCoords)
    plot_grid(111,polygonCoords)
    #cellObjs.append(polygonCoords)

# Write to file!
#plot_grid(111,cellObjs)
#print len(cellObjs[1])
#for eachPolygon in cellObjs:
#    print eachPolygon
