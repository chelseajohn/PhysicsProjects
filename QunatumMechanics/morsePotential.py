# Solution to Morse Potential

import numpy as np
import math as m
import matplotlib.pyplot as plt
from pylab import*

N = 1000
x = np.linspace(-0.8, 7., N)  # x values
h = x[1]-x[0]  # fundamental length being the stepsize
alpha = 0.01

V = np.zeros(N)
for i in range(N):
    V[i] = (m.exp(-2*x[i])-2*m.exp(-x[i]))/(alpha**2)


Mdd = 1./(h*h)*(np.diag(np.ones(N-1), -1) - 2 *
                np.diag(np.ones(N), 0) + np.diag(np.ones(N-1), 1))
H = -Mdd + np.diag(V)
E, eigenvector1 = np.linalg.eigh(H)
eigenvector = np.transpose(eigenvector1)

ylim(-11000, -40)
plt.ylabel('$\psi(x)$')
plt.xlabel('$x$')
plt.plot(x, V, color="Gray", label="V(x)")

dissen = E[99]-E[0]
print(dissen)

k = 9  # number of eigen states

for i in range(k):
    plt.plot(x, 800*eigenvector[i]+E[i],
             label="$Ealpha/2 _{}$={:>8.3f}".format(i, alpha*0.5*(-min(V)+E[i])))
    # plt.plot(x,10*eigenvector[i]+E[i],label="$E_{}$={:>8.3f}".format(i,alpha*0.5*(-min(V)+E[i])))


plt.title("Solution to morse potential")
plt.legend()
plt.grid()
plt.show()
