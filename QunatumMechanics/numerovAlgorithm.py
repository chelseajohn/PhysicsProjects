# Numerov Algorithm
import math as mt
import cmath as cmt
import numpy as np
import scipy
import matplotlib.pyplot as plt


def V(x_l, x_r, n):
    x = np.linspace(x_l, x_r, n+1)
    V = np.zeros(n+1)
    for i in range(0, n+1):
        # if -a<x[i]<a:
            # V[i]=-20
        V[i] = x[i]**2
    # plt.plot(x,V)
    # plt.show()
    return V


def numerov_step(psi1, psi2, k1, k2, k3, h):
    p = 2*(1.0 - 5*(h*h)*(k2)/12.)*psi2
    q = (1.0+1*(h*h)*(k1)/12.)*psi1
    o = 1.0 + 1*(h*h)*(k3)/12.
    t = (p-q)/o
    return t


E = 0  # estimated energy
x0 = -10
xE = 10
n = 1000
x = np.linspace(x0, xE, n+1)
dx = x[1]-x[0]
tol = 0.01
delta = 0.01
a = 1


for i in range(0, 20):
    f = 1.0
    while(f > tol):
        E = E+delta
        dx = x[1]-x[0]
        v = V(x0, xE, n+1)
        # p=np.where(abs(x-a)<0.001) for square
        p = np.where(abs(E-v) < 0.1)
        p = p[0][0]
        k = np.zeros(n+1)
        for i in range(0, len(k)):
            k[i] = (E-v[i])
        psi1 = np.zeros(p)
        psi2 = np.zeros(n+1-p)
        psi1[0] = 0
        psi1[1] = dx*0.1
        psi2[0] = 0
        psi2[1] = dx*0.1
        for j in range(2, p):
            psi1[j] = numerov_step(psi1[j-2], psi1[j-1],
                                   k[j-2], k[j-1], k[j], dx)
        for j in range(2, n-p+1):
            psi2[j] = numerov_step(psi2[j-2], psi2[j-1],
                                   k[n-(j-2)], k[n-(j-1)], k[n-j], dx)
        psi2 = np.flip(psi2, axis=0)
        k1 = psi1[p-1]
        k2 = psi2[0]
        ratio = k1/k2
        psi2 = ratio*psi2
        psi = np.concatenate((psi1, psi2))
        norm = np.sqrt(np.vdot(psi, psi)*dx)
        psi = psi/norm
        d1 = (psi[p-1]-psi[p-2])/dx
        d2 = (psi[p+1]-psi[p])/dx
        f = abs(d2-d1)
    print(E)
    plt.plot(x, psi)
    plt.show()
    E = E+tol*20
