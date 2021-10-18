import pylab
import numpy as np
from matplotlib import pyplot as plt
import math as math
from os import path
from progressbar import ProgressBar
progress=ProgressBar()

# Function to generate an image of the atomic configuration in the box.
def create_picture(positions,t):
    pylab.cla()
    pylab.axis([0, L, 0, L])
    pylab.setp(pylab.gca(), xticks=[0, L], yticks=[0, L])
    for x,y in positions:
        atom = pylab.Circle((x, y), Ratom, fc='r')
        plt.gca().add_patch(atom)
    plt.savefig(path.join('lennardimages/',"image-{0}.png".format(t)))   
   
###############calculating energy 
def energy(vel):
    k=0
    for i in range(Natoms):
        k=k+np.vdot(vel[i],vel[i])
    return 0.5*k

def potential(rel):
    r=np.sqrt(np.dot(rel,rel))
    p=4*((1/r)**12-(1/r)**6)
    return p

#####################
    
    
#########calculating acceleration ##
def acceleration(rel):
    r=np.sqrt(np.dot(rel,rel))
    a=(48/(r**2))*((1/r)**(12)-0.5*((1/r)**6))*rel 
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


m=10
Natoms=m*m #Number of atoms
dist_sep=1.2
L=(m+1)*dist_sep #length of box
r_c=5#critical length
Ratom=0.3#Radius of atom

#### Initial Positions of the atoms ####
pos=np.zeros([Natoms,2])
rel_pos=np.zeros([Natoms,2])
'''
i=-1
while(True):
    i=i+1
    pos[i,:]=np.random.uniform(Ratom,1-Ratom,[1,2])
    
    for j in range(0,i):
        if(np.sqrt(np.dot(pos[i,:]-pos[j,:],pos[i,:]-pos[j,:]))<2*Ratom):
            i=-1
            break
    if(i==Natoms-1):
        break
'''
##periodic arrangememt#####

count=0
columns=0
rows=0
while(True):
    x=dist_sep + dist_sep*columns
    y=dist_sep + dist_sep*rows
    pos[count][0]=x
    pos[count][1]=y
    columns += 1
    count=count+1
    if (L-x)<dist_sep:
        rows=rows+1
        columns=0
        count=count- 1   
    if(count==Natoms):
        break
 
create_picture(np.array(pos),0)
plt.show()

##initial velocity and accelerations of atoms##
vel=np.zeros([Natoms,2])
mean=0.0
sigma=1.0
vel[:,0]=np.random.normal(mean,sigma,Natoms)
vel[:,1]=np.random.normal(mean,sigma,Natoms)
cm_momentum=sum(vel)/Natoms
vel=vel-cm_momentum
a=np.zeros([Natoms,2])

T=10 #time
deltat=0.01  #delta time
n=int(T/deltat) #time steps
t=np.linspace(0,T,n)
kinetic_energy=np.zeros(n)
potential_energy=np.zeros(n)
Total_energy=np.zeros(n)
E_cutoff=4*((1/r_c)**12-(1/r_c)**6)
k_1=energy(vel)/Natoms
T_target=0.2
Temp=k_1
v_list=[]
#######verlet algorithm######
for g in progress(range(n)):   
    pos=pos+vel*deltat*0.5
    pos=correct_pos(pos)
    for i in range(Natoms-1):
        for j in range(i+1,Natoms):
            rel=pos[i]-pos[j]
            relative=correct_r(rel)          
            if np.sqrt(np.dot(relative,relative)) < r_c:
                a[i] += acceleration(relative)
                a[j] += -acceleration(relative)
                #potential_energy[g] +=potential(relative)-E_cutoff
    vel=vel + a*deltat
    pos=pos + vel*0.5*deltat 
    for k in range(Natoms):         
         v_list.append(np.sqrt(np.vdot(vel[k],vel[k])))
    pos=correct_pos(pos)
    for i in range(Natoms-1):
        for j in range(i+1,Natoms):
            rel=pos[i]-pos[j]
            relative=correct_r(rel)          
            if np.sqrt(np.dot(relative,relative)) < r_c:
                potential_energy[g] +=potential(relative)-E_cutoff
    a=np.zeros([Natoms,2])
    kinetic_energy[g]=energy(vel)
    k_mean=kinetic_energy[g]/Natoms
    Temp=k_mean
    vel=vel*np.sqrt(T_target/Temp)  
    #if(g%10==0):
        #create_picture(np.array(pos),int(g/10))

Total_energy=kinetic_energy+potential_energy

plt.plot(t,kinetic_energy,label="kinetic energy")
plt.plot(t,potential_energy,label="potential energy")
plt.plot(t,Total_energy,label="total energy")
plt.legend()
plt.show()
v=np.linspace(0,max(v_list),100)
maxwell=v*np.exp(-v*v/(2*T_target))/T_target
plt.plot(v,maxwell)   
plt.hist(np.array(v_list),numbins=10,normed=1)
plt.show()

   















