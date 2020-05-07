import numpy as np
import matplotlib.pyplot as plt

x1= 4
y1 = 4.5

x2=7
y2=7.2


dist = np.sqrt(np.power(x1-x2, 2) + np.power(y1-y2, 2))

xcen = (x1+x2)*0.5
ycen = (y1+y2)*0.5

angle = np.arctan2(x1-x2, y1-y2);


newx1 = xcen + np.cos(angle - 0.5*np.pi)*dist*0.75 +10
newx2 = xcen - np.cos(angle - 0.5*np.pi)*dist*0.75+10

newy1 = ycen + np.sin(angle - 0.5*np.pi)*dist*0.75+10
newy2 = ycen - np.sin(angle - 0.5*np.pi)*dist*0.75+10

dist2 = np.sqrt(np.power(newx1-newx2, 2) + np.power(newy1-newy2, 2))

print("angle, newangle1, newangle2, dist, dist2")
print(angle, angle - np.pi, angle + np.pi, dist, )


m1 = (y2-y1)/(x2-x1)
b1 = y2 - x2*m1


m2 = (newy2-newy1)/(newx2-newx1)
b2 = newy2 - newx2*m2
print("m1, b1, m2, b2")
print(m1, b1, m2, b2)

cross = (b2-b1)/(m1-m2)

y1cross = b1 + m1*cross
y2cross = b2 + m2*cross


xx = [x1, x2, newx1, newx2, xcen, cross, cross]
yy = [y1, y2, newy1, newy2, ycen, y1cross, y2cross]
nn = ['p1', 'p2', 'new p1', 'new p2', 'center', '     cr1', '        c2']

if(((cross >= x1 and cross <= x2) or (cross >= x2 and cross <= x1)) and ((cross >= newx1 and cross <= newx2) or (cross >= newx2 and cross <= newx1))):
    print("Lines cross")
else:
    print("lines do not cross")

fig, ax = plt.subplots()
ax.scatter(xx, yy)

for i, tx in enumerate(nn):
	ax.annotate(tx, (xx[i], yy[i]))


plt.plot([x1, x2], [y1, y2])
plt.plot([newx1, newx2], [newy1, newy2])
plt.show()
