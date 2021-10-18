# Solution to coulomb potential


import numpy as np
import math as m
import matplotlib.pyplot as plt
from pylab import*

N = 1000  # number of points
x = np.linspace(.1, 10., N)  # x values
h = x[1]-x[0]  # fundamental length being the stepsize

V = np.zeros(N)  # potential values
for i in range(N):
    V[i] = -2/x[i]

Mdd = 1./(h*h)*(np.diag(np.ones(N-1), -1) - 2 *
                np.diag(np.ones(N), 0) + np.diag(np.ones(N-1), 1))
H = -Mdd + np.diag(V)
E, eigenvector1 = np.linalg.eigh(H)
eigenvector = np.transpose(eigenvector1)

# print(eigenvector)


plt.plot(x, V, color="Gray", label="V(x)")
k = 4  # number of eigen states

for i in range(k):
    # plt.plot(x,np.full(N,E[i]),label="$E_{}$={:>8.3f}".format(i,E[i]))
    plt.plot(x, 10*eigenvector[i]+E[i],
             label="$E_{}$={:>8.3f}".format(i, E[i]))


plt.title("Solution to coulomb potential")
plt.legend()
plt.grid()
plt.show()
