
import numpy as np
import matplotlib.pyplot as plt

t1a = lambda a, a0 : 0.5*(a/a0 - 1)**2
t1b = lambda t, l, a0: t*l/np.sqrt(a0)
t1c = lambda s, l, a0: l**2 * s/(2*a0)

t2a = lambda a, ao : 0.5*(a - a0)**2
t2b = lambda t, l: t*l
t2c = lambda s, l: l**2 * s/2

#These work for step1 (conv ext)
a0=12
t=-0.02
s=0.03

ap = np.linspace(0,24, 1000)
lp = np.linspace(0,10, 1000)
Lp = np.linspace(0,12, 1000)

sol1a = t1a(ap, a0)
sol1b = t1b(t, lp, a0)
sol1c = t1c(s, Lp, a0)

sol2a = t2a(ap, a0)*0.00694
sol2b = t2b(t, lp)*0.2887
sol2c = t2c(s, Lp)*0.0833

plt.plot(sol1a, sol2a);plt.show()
plt.plot(sol1b, sol2b);plt.show()
plt.plot(sol1c, sol2c);plt.show()


np.unique(sol1a/sol2a)
np.unique(sol1b/sol2b)
np.unique(sol1c/sol2c)




#These work for step1 (conv ext)
a0=12
t=-0.02
s=0.03

ap = np.linspace(0,24, 1000)
lp = np.linspace(0,10, 1000)
Lp = np.linspace(0,12, 1000)

sol1a = t1a(ap, a0)
sol1b = t1b(t, lp, a0)
sol1c = t1c(s, Lp, a0)

sol2a = t2a(ap, a0)*0.00694 
sol2b = t2b(t, lp)*0.2887
sol2c = t2c(s, Lp)*0.0833

plt.plot(sol1a, sol2a);plt.show()
plt.plot(sol1b, sol2b);plt.show()
plt.plot(sol1c, sol2c);plt.show()


np.unique(sol1a/sol2a)
np.unique(sol1b/sol2b)
np.unique(sol1c/sol2c)


#These work for step2 (exp)
a0=60
t=-0.02
s=0.03

ap = np.linspace(0,24, 1000)
lp = np.linspace(0,10, 1000)
Lp = np.linspace(0,12, 1000)

sol1a = t1a(ap, a0)
sol1b = t1b(t, lp, a0)
sol1c = t1c(s, Lp, a0)

sol2a = t2a(ap, a0)
sol2b = t2b(t, lp)
sol2c = t2c(s, Lp)

#plt.plot(sol1a, sol2a);plt.show()
#plt.plot(sol1b, sol2b);plt.show()
#plt.plot(sol1c, sol2c);plt.show()


np.unique(sol1a/sol2a)
np.unique(sol1b/sol2b)
np.unique(sol1c/sol2c)
# For a0 = 35: 0.000816, 0.169, 0.02857
# For a0 = 60: 0.000278, 0.129, 0.01667

#These work for step3 (hinge contract) - BLADE
a0=30
t=-0.02
s=0.03

ap = np.linspace(0,24, 1000)
lp = np.linspace(0,10, 1000)
Lp = np.linspace(0,12, 1000)

sol1a = t1a(ap, a0)
sol1b = t1b(t, lp, a0)
sol1c = t1c(s, Lp, a0)

sol2a = t2a(ap, a0)
sol2b = t2b(t, lp)
sol2c = t2c(s, Lp)

#plt.plot(sol1a, sol2a);plt.show()
#plt.plot(sol1b, sol2b);plt.show()
#plt.plot(sol1c, sol2c);plt.show()


np.unique(sol1a/sol2a)
np.unique(sol1b/sol2b)
np.unique(sol1c/sol2c)


#These work for step3 (hinge contract) - HINGE
a0=10
t=-0.02
s=0.03

ap = np.linspace(0,24, 1000)
lp = np.linspace(0,10, 1000)
Lp = np.linspace(0,12, 1000)

sol1a = t1a(ap, a0)
sol1b = t1b(t, lp, a0)
sol1c = t1c(s, Lp, a0)

sol2a = t2a(ap, a0)
sol2b = t2b(t, lp)
sol2c = t2c(s, Lp)

#plt.plot(sol1a, sol2a);plt.show()
#plt.plot(sol1b, sol2b);plt.show()
#plt.plot(sol1c, sol2c);plt.show()


np.unique(sol1a/sol2a)
np.unique(sol1b/sol2b)
np.unique(sol1c/sol2c)

