from pylab import *
import pylab
import numpy as np
import random
import math

L=20
N=L**2#number of spins
J=1.0#coupling strength


def create_picture(positions,colors): #graphical visualisation of spins
    pylab.cla()
    pylab.axis([0,L,0,L])
    pylab.setp(pylab.gca())
    for pos,col in zip(positions,colors):
        square=pylab.Rectangle((pos[0],pos[1]),0.8,0.8,fc=col)
        pylab.gca().add_patch(square)

def color(i):
    if i==1: return 'r' #Red:'up'
    else: return 'b' #Blue: 'down'

 #############################computes neighbours for ith spin ###############
#########################identification of opposite edges is taken into account###########   
def right(i):#neighbor to the right
     if(i+1)%L==0:return i+1-L
     else: return i+1
def left(i):#neighbor to the left
    if i%L==0:return i-1+L
    else: return i-1
def up(i):#neighbor up
    return(i+L)%N

def down(i):#neighbor down
    return (i-L+N)%N

neighbors=[[right(i),left(i),up(i),down(i)] for i in range(N)] #neighbor table
coordinates=[[i%L,i//L] for i in range(N)] ##coordinates of spins(squares) for graphical output

orientations=[-1,1]
spins=[]
colors=[]

############initial random spin assignment #########]

for i in range(N):
    spin=random.choice(orientations)
    spins.append(spin)
    
###################################
nsteps=1000*N ##number of steps in Metropolis algorithm
T=1.0     #temperature measured in units of J/K_B
beta=1/T

##########Metropolis algorithm
for step in range(nsteps):
    k=random.randint(0,N-1) ######random spin choice
    delta_E=2.0*spins[k]*sum(spins[j] for j in neighbors[k]) ##change in energy
    if random.uniform(0.0,1.0) < math.exp(-beta*delta_E):#metropolis acceptance probability
        spins[k]*=-1

for i in range(N):
    colors.append(color(spins[i]))

create_picture(coordinates,colors)
pylab.show()





















