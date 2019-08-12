

import numpy as np
import matplotlib.pyplot as plt

x1= 4
y1 = 4.5

x2=7
y2=7.2

dist = np.sqrt(np.power(x1-x2, 2) + np.power(y1-y2, 2))

xcen = (x1+x2)*0.5
ycen = (y1+y2)*0.5

angle = np.arctan2(y1-y2, x1-x2);
angle2 = np.arctan2(y2-ycen, x2-xcen); #not used


#newx1 = xcen + np.cos(angle - 0.5*np.pi)*dist*0.75
#newx2 = xcen - np.cos(angle - 0.5*np.pi)*dist*0.75

#newy1 = ycen + np.sin(angle - 0.5*np.pi)*dist*0.75
#newy2 = ycen - np.sin(angle - 0.5*np.pi)*dist*0.75


newx1 = xcen + np.cos(angle + 2*0.5*np.pi)*dist*0.75 -  0*np.sin(angle + 2*0.5*np.pi)*dist*0.75
newx2 = xcen - np.cos(angle + 2*0.5*np.pi)*dist*0.75 -  0*np.sin(angle2 - 2*0.5*np.pi)*dist*0.75

newy1 = ycen + np.sin(angle + 2*0.5*np.pi)*dist*0.75  + 0*np.cos(angle + 2*0.5*np.pi)#*dist*0.75
newy2 = ycen - np.sin(angle + 2*0.5*np.pi)*dist*0.75 + 0*np.cos(angle2 - 2*0.5*np.pi)#*dist*0.75



#newx1 = xcen + np.cos(angle )*dist*0.75 -  0*np.sin(angle + 2*0.5*np.pi)*dist*0.75
#newx2 = xcen + np.cos(angle2 )*dist*0.75 -  0*np.sin(angle2 - 2*0.5*np.pi)*dist*0.75

#newy1 = ycen + np.sin(angle )*dist*0.75  + 0*np.cos(angle + 2*0.5*np.pi)#*dist*0.75
#newy2 = ycen + np.sin(angle2)*dist*0.75 + 0*np.cos(angle2 - 2*0.5*np.pi)#*dist*0.75



newx3= xcen + np.cos(angle  + 2*0.5*np.pi)*dist*0.75 -  0*np.sin(angle + 2*0.5*np.pi)*dist*0.75
newx4 = xcen - np.cos(angle  + 2*0.5*np.pi)*dist*0.75 -  0*np.sin(angle2 - 2*0.5*np.pi)*dist*0.75

newy3 = ycen + np.sin(angle  + 2*0.5*np.pi)*dist*0.75  + 0*np.cos(angle + 2*0.5*np.pi)#*dist*0.75
newy4 = ycen - np.sin(angle + 2*0.5*np.pi)*dist*0.75 + 0*np.cos(angle2 - 2*0.5*np.pi)#*dist*0.75



dist2 = np.sqrt(np.power(newx1-newx2, 2) + np.power(newy1-newy2, 2))

print("angle, angle2, newangle1, newangle2, dist, dist2")
print(angle,angle2, angle + np.pi, angle2 + np.pi, dist, )


xx = [x1, x2, newx1, newx2, xcen, newx3, newx4]
yy = [y1, y2, newy1, newy2, ycen, newy3, newy4]
nn = ['p1', 'p2', 'new p1', 'new p2', 'center',  'new p3', 'new p4']


fig, ax = plt.subplots()
ax.scatter(xx, yy)

for i, tx in enumerate(nn):
	ax.annotate(tx, (xx[i], yy[i]))


plt.plot([x1, x2], [y1, y2])
plt.plot([newx1, newx2], [newy1, newy2])
plt.plot([newx3, newx4], [newy3, newy4])
plt.show()
