

from pylab import *
import pylab
import numpy as np
import random
import matplotlib.pyplot as plt
from tqdm import tqdm
from progressbar import ProgressBar
progress=ProgressBar()

L = 20
N = L**2 # Number of spins
J = 1.0 # Coupling strength


def create_picture(positions,colors): ## Graphical visualisation of spins
    pylab.cla()
    pylab.axis([0, L, 0, L])
    pylab.setp(pylab.gca())
    for pos, col in zip(positions,colors):
        square = pylab.Rectangle((pos[0], pos[1]), 0.8, 0.8, fc = col)
        pylab.gca().add_patch(square)
        
        
def color(i):
    if i == 1: 
        return 'r'# Red: `up'
    else: 
        return 'b' # Blue: `down'
    
##### ######################### Computes neighbors for ith spin #########
############## Identification of opposite edges is taken into account #####
def right(i): # Neighbor to the right
    if (i+1)%L == 0: return i+1-L
    else: return i+1
def left(i): # Neighbor to the left
    if i%L == 0: return i-1+L
    else: return i-1
def up(i): # Neighbor up
    return (i+L)%N
def down(i): # Neighbor down
    return (i-L+N)%N

def wolfclust(mag):
 
    t=random.randint(0, N - 1) ###random spin choice
    cluster=[t]
    frontier=[t]
    while np.size(frontier)!=0:
        j=random.choice(frontier)
        for k in neighbors[j]:
            p=random.uniform(0.0, 1.0)
            if spins[k]==spins[j] and k not in cluster and p<prob:
                cluster.append(k)
                frontier.append(k)
        frontier.remove(j)
    
    for p in cluster:
            spins[p] *=-1
        
    mag=mag+abs(sum(spins))
    magtemp=abs(sum(spins))
    return mag,magtemp
                
            
    
    
    
neighbors = [[right(i),left(i),up(i),down(i)] for i in range(N)] # Neighbor table
##########################################################################
coordinates = [[i%L,i//L] for i in range(N)] ## Coordinates of spins (squares) for graphical output
orientations = [1,-1]
n=1000
spins = []
colors = []
##### initial random spin assignment ##################
for p in range(N):
        spin= random.choice(orientations)
        spins.append(spin)
##################################################


Temp=np.arange(1,4,0.1)
M=np.zeros(len(Temp))
dM=np.zeros(len(Temp))


for i in progress(range(len(Temp))):
    mag=0
    count=0
    T=Temp[i]
    beta = 1/T
    prob=1-np.exp(-2/T)   
    Mtemp=[]
    while(count<n):
        count=count+1
        mag,mtemp=wolfclust(mag)
        Mtemp.append(mtemp)
    dM[i]=(np.std(Mtemp))**2/(N*T)
    M[i]=(float(mag)/(N*count))
    
    
plt.plot(Temp[1:len(Temp)],dM[1:len(Temp)])
plt.show()
plt.plot(Temp,M)
plt.show()
'''
for j in range(N):
    colors.append(color(spins[j]))
create_picture(coordinates,colors)
pylab.show()
'''
    
    