import numpy as np
import matplotlib.pyplot as plt


nplots = 10
ff = lambda x, sp : 1 - np.exp(-(x**sp))


maxy = ff(1,1) - 0.05
a = np.linspace(0, 1, 1000)
b = np.linspace(0.1, 5, nplots)

fig = plt.figure()
ax = plt.subplot(111)

for i, x in enumerate(b):
    c = ff(a, x)
    ax.plot(a, c, label="b = " + str(round(x, 2)))
    #plt.text(x=0.1, y=1-maxy*i/nplots, s = str(round(x, 2)))

plt.title('1 -  exp(-x^b)')
chartBox = ax.get_position()
ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.8, chartBox.height])
ax.legend(loc='upper center', bbox_to_anchor=(1.2, 0.8), shadow=True, ncol=1)
plt.show()
