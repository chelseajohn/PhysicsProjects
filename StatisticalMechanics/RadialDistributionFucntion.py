from pylab import *
import pylab
import numpy as np
from matplotlib import pyplot as plt
import math as m
from os import path
from progressbar import ProgressBar
progress=ProgressBar()


###############calculating energy 
def energy(vel):
    k=0
    for i in range(Natoms):
        k=k+np.vdot(vel[i],vel[i])
    return 0.5*k

#####################potential energy
def potential(rel):
    r=np.sqrt(np.dot(rel,rel))
    p=4.0*((1/r)**12-(1/r)**6)
    return p

#########virial pressure
def vp(rel):
    r=np.sqrt(np.dot(rel,rel))
    pr=48.0*((1/r)**12-0.5*(1/r)**6)
    return pr
   
#########calculating acceleration ##
def acceleration(rel):
    r=np.sqrt(np.dot(rel,rel))
    a=(48.0/(r**2))*((1/r)**(12)-0.5*((1/r)**6))*rel 
    return a

##########correction to seperation##
def correct_r(rel):
    if(np.abs(rel[0])>0.5*L):
        rel[0]=rel[0]-L*np.sign(rel[0]) 
    if(np.abs(rel[1])>0.5*L):
        rel[1]=rel[1]-L*np.sign(rel[1])   
    return rel

##########correction to position##
def correct_pos(positions):
    for i in range(Natoms):
        if positions[i][0]>L:
            positions[i][0]=positions[i][0]-L
        if positions[i][0]<0:
            positions[i][0]=positions[i][0]+L
        if positions[i][1]>L:
            positions[i][1]=positions[i][1]-L
        if positions[i][1]<0:
            positions[i][1]=positions[i][1]+L
    return positions

####radial function###
def RDF(positions,length):
    distribution_array=np.zeros(length)
    for i in range(Natoms-1):
        for j in range(i+1,Natoms):
            rel=pos[i]-pos[j]
            relative=correct_r(rel)     
            pair_distance=np.sqrt(np.dot(relative,relative))
            ind=int(pair_distance/bin_size)
            if ind <= (length-1):
                distribution_array[ind]+=1
                
    for i in range(length):
        distribution_array[i]= distribution_array[i]/(rho*2*np.pi*(i+1)*bin_size*bin_size)
    return distribution_array           
    
                
    
data=np.loadtxt('data.dat')
L=data[0,4]
r_c=3.0#critical length
bin_size=0.1
length_radial=int(0.5*L/bin_size) #length of RDF array
qw=0
arr=np.zeros([length_radial])
#### Initial Positions of the atoms ####
pos=data[:,0:2]
Natoms=pos.shape[0]
rel_pos=np.zeros([Natoms,2])


####initial velocities#
vel=data[:,2:4]
    
a=np.zeros([Natoms,2])

T=10 #time
deltat=0.01  #delta time
n=int(T/deltat) #time steps
t_relax=1.0
check_time=0.0

kinetic_energy=np.zeros(n)
potential_energy=np.zeros(n)
virial=np.zeros(n)
Total_energy=np.zeros(n)
Inst_temp=np.zeros(n)

rho=Natoms/(L*L)

E_cutoff=4.0*((1.0/r_c)**12-(1/r_c)**6)
k_1=energy(vel)/Natoms


Temp=k_1
v_list=[]

#######verlet algorithm######
for g in  progress(range(n)):   
    pos=pos+vel*deltat*0.5
    pos=correct_pos(pos)
    check_time= g*deltat
    for i in range(Natoms-1):
        for j in range(i+1,Natoms):
            rel=pos[i]-pos[j]
            relative=correct_r(rel)          
            if np.sqrt(np.dot(relative,relative)) < r_c:
                a[i] += acceleration(relative)
                a[j] += -acceleration(relative)
                
    vel=vel + a*deltat
    pos=pos + vel*0.5*deltat 
    if check_time > T*0.5:
        for k in range(Natoms):         
             v_list.append(np.sqrt(np.vdot(vel[k],vel[k])))
    pos=correct_pos(pos)
    
    for i in range(Natoms-1):
        for j in range(i+1,Natoms):
            rel=pos[i]-pos[j]
            relative=correct_r(rel)          
            if np.sqrt(np.dot(relative,relative)) < r_c:
                potential_energy[g] +=potential(relative)-E_cutoff
                virial[g] +=vp(relative)
                
    a=np.zeros([Natoms,2])
    
    kinetic_energy[g]=energy(vel)
    Temp=kinetic_energy[g]/Natoms
    Inst_temp[g]=Temp
    if g%10==0 :
        arr+=RDF(pos,length_radial)
        qw=qw+1
   
arr=arr/qw
Total_energy=kinetic_energy+potential_energy
t=np.linspace(0,T,n)#time to  plot 

distance=np.linspace(0,0.5*L,length_radial)
plt.title("Radial Distribution Function")
plt.xlabel("r")
plt.ylabel("u(r)")
plt.plot(distance,arr)
plt.show()


       

    
    
    
    
    
    
    
    
    
    
    
    
    
