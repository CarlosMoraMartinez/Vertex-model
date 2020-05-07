
import numpy as np
import matplotlib.pyplot as plt

l = []
with open("test_random_1.tsv") as ff:
	for line in ff:
		l.append([float(x) for x in line.split("\t")])

a = np.array([x[0] for x in l])
b = np.array([x[1] for x in l]) 
c = np.array([x[2] for x in l])  

plt.scatter(a, b)
plt.show()
plt.scatter(b, c)
plt.show()
	
