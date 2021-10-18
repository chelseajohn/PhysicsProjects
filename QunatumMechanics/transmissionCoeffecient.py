# Calculating Transmission Coeffecients for potential well

import matplotlib.pyplot as plt
import numpy as np
import cmath
import math
pi = cmath.pi
a = -50
b1 = 50
N = 1000  # number of points

V = np.zeros(N)
x = np.linspace(a, b1, N)

for j in range(N):  # potential well
    if(abs(x[j]) < 0.5):
        V[j] = -50

delta = float(((b1-a)/N))

H = np.zeros([N, N])
D = np.zeros([N, N])

for i in range(N):
    for j in range(N):
        if(i == j):
            D[i][j] = -2
        if(np.abs(i-j) == 1):
            D[i][j] = 1

for i in range(N):
    for j in range(N):
        if(i == j):
            H[i][j] = -1/(delta*delta)*D[i][j] + V[i]
        else:
            H[i][j] = -1/(delta*delta)*D[i][j]

evals, evect = np.linalg.eigh(H)


b = -10  # Peak of gaussian
w = 2  # width of gaussian

p0 = 1  # initial momentum
deltap = 0.2  # momentum stepsize

N2 = 80
p = np.zeros(N2)  # momentum values
T = np.zeros(N2, dtype=complex)  # transmission coeffecient

for j in range(N2):
    phi = np.zeros([N, 100], dtype=complex)

    p0 = p0 + deltap
    p[j] = p0

    maxt = (100/p0)

    for i in range(N):

        phi[i, 0] = (1/(2*pi*w*w)**(1/4)) * \
            math.exp(-(((x[i]-b)**2)/(4*w*w)))*cmath.exp(1j*p0*x[i])

    l = (phi[:, 0])
    c = np.zeros(N, dtype=complex)

    for i in range(N):  # finding Ci's
        c[i] = np.vdot(evect[:, i], l)

    t = np.linspace(0, maxt, 100)

    summ = 0
    for i in range(100):

        for k in range(N):
            summ = summ + c[k]*cmath.exp(-1j*evals[k]*t[i])*evect[:, k]

        phi[:, i] = summ
        summ = 0

    new = phi[:, 25]
    # plt.plot(x,(20*np.abs(phi[:,18]))**2)
    T[j] = np.vdot(new[int(N/2):N], new[int(N/2):N]) * \
        delta  # calculating transmission coeffecient

plt.plot(p, T)

# verification of results using theoratical expression of transmission cofficient

T2 = np.zeros(N2)
V0 = 50
a1 = 0.5
for i in range(N2):
    T2[i] = (16*(p[i]**2 + V0)*p[i]*p[i])/(16*p[i]*p[i]*(p[i]**2 + V0) + 4*V0*V0 *
                                           math.sin(2*a1*math.sqrt(p[i]**2 + V0))*math.sin(2*a1*math.sqrt(p[i]**2 + V0)))


plt.plot(p, T2)
plt.title(" transmission coeffecients vs energy")
plt.xlabel("energy")
plt.ylabel(" transmission coeffecient")
plt.grid()
plt.legend()
plt.show()
