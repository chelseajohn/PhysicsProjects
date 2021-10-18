
# Solution to Finite Square Well
import numpy as np
import matplotlib.pyplot as plt
from pylab import*

N = 4500#number of points
a = 200.0
b = 2.
x = np.linspace(-a/2.,a/2.,N)#x values
h = x[1]-x[0] #fundamental length being the stepsize

V0 =-600.# minimum potential
V=np.zeros(N)#potential values
#defining the potential
for i in range(N):
    if x[i]> -b/2. and x[i]< b/2.:
        V[i]= V0

Mdd = 1./(h*h)*(np.diag(np.ones(N-1),-1) -2* np.diag(np.ones(N),0) + np.diag(np.ones(N-1),1))
H = -Mdd + np.diag(V)#matrix element of the hamiltonian
E,eigenvector1=np.linalg.eigh(H)#finding the eigen values and eigen vectors
eigenvector = np.transpose(eigenvector1) 

#print(eigenvector)

plt.figure(figsize=(10,7))
plt.xlim((-5*b,5*b))
plt.plot(x,V,color="Gray",label="V(x)")#plotting the potential
k=5#number of eigen states
for i in range(k):
    #plt.plot(np.linspace(-b/2.,b/2.,N),np.full(N,E[i]),label="$E_{}$={:>8.3f}".format(i,E[i]))
     plt.plot(x,10*eigenvector[i]+E[i],label="$E_{}$={:>8.3f}".format(i,E[i]))
      
          
        
plt.title("Solutions to the Finite Square Well")
plt.legend()
plt.grid()
plt.show()