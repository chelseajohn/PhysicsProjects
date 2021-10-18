
import pylab
import numpy as np
from matplotlib import pyplot as plt
import math as m
from os import path
from progressbar import ProgressBar
progress=ProgressBar()
from tqdm import tqdm

# Function to generate an image of the atomic configuration in the box.
def create_picture(positions,t):
    pylab.cla()
    pylab.axis([0, L, 0, L])
    pylab.setp(pylab.gca(), xticks=[0, L], yticks=[0, L])
    for x,y in positions:
        atom = pylab.Circle((x, y), Ratom, fc='r')
        plt.gca().add_patch(atom)
    plt.savefig(path.join('virialimages/',"image-{0}.png".format(t)))   
   
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

u=3
P=np.linspace(0.1,1.5,u)
rho_avg=np.zeros(u)
no=5
dist_sep=1.0
Natoms=no*no #Number of atoms
L=(no+1)*dist_sep #length of box
r_c=3.0#critical length
Ratom=0.3#Radius of atom

#### Initial Positions of the atoms ####
pos=np.zeros([Natoms,2])
rel_pos=np.zeros([Natoms,2])
'''
i=-1
while(True):
    i=i+1
    pos[i,:]=np.random.uniform(Ratom,1-Ratom,[1,2])8
    
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

for co in progress(range(u)):

    
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
    Inst_press=np.zeros(n)
    rho=np.zeros(n)
    
    E_cutoff=4.0*((1.0/r_c)**12-(1/r_c)**6)
    k_1=energy(vel)/Natoms
    
    T_target=1.0/6.0
    Temp=k_1
    v_list=[]
    
    #######verlet algorithm######
    for g in tqdm(range(n)):   
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
        
        pressure=Natoms*(Temp/L**2)+virial[g]/(2*L*L)
        Inst_press[g]=pressure
        P_target=P[co]
        rho[g]=Natoms/L**2
        
        lamb=np.sqrt(1+((T_target/Temp)-1)*(deltat/t_relax))
        mu=np.sqrt(1+(pressure-P_target)*(deltat/t_relax))
        
        vel=vel*lamb
        pos=pos*mu
        L=L*mu
        
        
    rho_avg[co]=float(sum(rho))/n
    #create_picture(np.array(pos),co)
    #plt.show()
    Total_energy=kinetic_energy+potential_energy
    t=np.linspace(0,T,n)#time to  plot 
    
plt.plot(rho_avg,P,label="P vs Density")
plt.legend()
plt.show()
       
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
