import numpy as np
import matplotlib.pyplot as plt


a=np.linspace(0, 360, 10000)
b = np.sin(np.pi*a/180)
plt.plot(a, b, color = "orange")
plt.scatter([90, 180, 270, 360], np.sin([0.5*np.pi, np.pi, 1.5*np.pi, 2*np.pi]))
plt.show()





def getStrength(angle, minstr, maxstr, maxangle):
    angle = np.pi*(90 + abs(angle - maxangle))/180
    prop = abs(np.sin(angle))
    res = minstr + prop*(maxstr - minstr)
    return res


def plotExample(minstr, maxstr, maxangle):
    a = np.linspace(maxangle-180, maxangle+180, 10000)
    b = np.array([getStrength(i, minstr, maxstr, maxangle) for i in a])
    plt.plot(a, b, color = "blue")
    plt.hlines([minstr, maxstr], xmin = min(a), xmax=max(a), color = "black", linestyles = "dashed")
    #plt.vlines(x = [-180, -90, 0, 90, 180, 270, 360], ymin=minstr, ymax = maxstr, color = "green", alpha = 0.2)
    plt.vlines(x = [maxangle-180, maxangle, maxangle+180], ymin=minstr, ymax = maxstr, color = "green", alpha = 0.5, label = "Maximum")
    plt.vlines(x = [maxangle-90, maxangle+90], ymin=minstr, ymax = maxstr, color = "red", alpha = 0.5, label = "Minimum")
    plt.legend()
    plt.show()


def plotLines(minstr, maxstr, maxangle, n=100):
    for i in range(n):
        x = np.random.rand(2)
        y = np.random.rand(2)
        angle = 180*np.arctan2(y[1]-y[0], x[1]-x[0])/np.pi
        res = getStrength(angle, minstr, maxstr, maxangle)
        plt.plot(x, y, alpha = res, color = "black")
    plt.show()
	


plotExample(0.1, 7.0, 90)
plotExample(-0.01, 0.02, 80)
plotExample(-0.03, -0.02, -45)
plotExample(0, 25, 0)
plotExample(20, 2, 0)


plotLines(0.2, 1.0, 90)
plotLines(0.2, 1.0, 0)
plotLines(0.2, 1.0, 45)
plotLines(0.0, 0.2, 45)


