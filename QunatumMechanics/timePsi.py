# Time evolving schrodinger equation for potential well
import numpy as np
import matplotlib.pyplot as plt
from pylab import*
import cmath
import math

N = 1000
pi = math.pi
n = 100
b = -10.0  # peak of gaussian
p0 = 2.0  # momentum
w = 2.0  # width of gaussian
s = 0.0
psi = np.zeros([N, n], dtype=complex)
C = np.zeros(N, dtype=complex)
t1 = np.linspace(0., 40., n)
x = np.linspace(-50, 50, N)
h = x[1]-x[0]
V0 = -50.  # minimum potential
V = np.zeros(N)  # potential values

for i in range(N):
    psi[i, 0] = 100*(cmath.exp(-(x[i]-b)**2/(4*w*w) +
                               (-1j*p0*x[i])))/((2*pi*w*w)**(1/4))
    if x[i] > -.5 and x[i] < .5:  # defining potential well
        V[i] = V0


Mdd = 1./(h*h)*(np.diag(np.ones(N-1), -1)-2 *
                np.diag(np.ones(N), 0) + np.diag(np.ones(N-1), 1))
H = -Mdd + np.diag(V)
E, eigenvector = np.linalg.eigh(H)  # eigen values


for i in range(N):  # calculating Ci's
    C[i] = np.vdot(psi[:, 0], np.conjugate(eigenvector[:, i]))

for j in range(n):
    for i in range(N):
        s = s+C[i]*cmath.exp(-1j*E[i]*t1[j])*eigenvector[:, i]
    psi[:, j] = s
    s = 0.0


p = np.vdot(psi[:, 2]/100, psi[:, 2]/100)*h
print(p)
R = np.vdot(psi[0:500, 20]/100, psi[0:500, 20]/100)*h  # reflection coeffecient
T = np.vdot(psi[500:1000, 20]/100, psi[500:1000, 20]/100) * \
    h  # transmission coeffecient
print(np.abs(R))
print(np.abs(T))

plt.plot(x, np.abs(psi[:, 20]**2))
plt.title("Time evolving schrodinger equation")
plt.xlabel("x")
plt.ylabel("Probability($\psi(x)^2$)")
plt.grid()
plt.plot(x, 10*V, color="Gray", label="V(x)")
plt.legend()
plt.show()
