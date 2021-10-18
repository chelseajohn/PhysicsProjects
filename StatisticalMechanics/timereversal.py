
############## EVENT-DRIVEN COLLISIONS ####################################
###########################################################################

#### This Program simulates the motion of four atoms in a 2D box ######

import pylab
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation


#####  Function to compute time for wall collision ######
def wall_time(coord, velcomp, rad):  # pos1 and pos2 are positions of atoms 1 and 2, vel1 and vel2 are their 

    if velcomp > 0.0:
        del_t = (1.0 - rad - coord) / velcomp
    elif velcomp < 0.0:
        del_t = (coord - rad) / abs(velcomp)
    else:
        del_t = float('inf')
    return del_t


# Function to calculate time it takes for a pair of atoms to collide. 
# velocities. rad is the radius of the atoms. 
def pair_time(pos1, vel1, pos2, vel2, rad):  
    rel_pos = pos2 - pos1
    rel_vel = vel2 - vel1
    rel_dist_squar = np.dot(rel_pos,rel_pos)
    rel_speed_squar = np.dot(rel_vel,rel_vel)
    scal_prod = np.dot(rel_pos,rel_vel)
    a = scal_prod ** 2 - rel_speed_squar * ( rel_dist_squar - 4.0 * rad **2)
    if a > 0.0 and scal_prod < 0.0: ## Conditions for collision.
        del_t = - (scal_prod + np.sqrt(a)) / rel_speed_squar ## Collision time.
    else:
        del_t = float('inf')
    return del_t


# Function to generate an image of the atomic configuration in the box.
def create_picture(positions):
    pylab.cla()
    pylab.axis([0, L, 0, L])
    pylab.setp(pylab.gca(), xticks=[0, L], yticks=[0, L])
    for x,y in positions:
        atom = pylab.Circle((x, y), Ratom, fc='r')
        pylab.gca().add_patch(atom)


##########################################################################


######################### Initialization #################################
L = 1. ## Box edge length
Natoms = 10 #Number of atoms
density = 0.2 # Fraction of area of box occupied by the atoms
Ratom = np.sqrt(density/(Natoms*np.pi))   ## Radius of an atom.


#### Initial Positions and velocities of the atoms ####
pos=np.zeros([Natoms,2])
vel=np.zeros([Natoms,2])
i=-1
while(True):
    i=i+1
    pos[i,:]=np.random.uniform(Ratom,1-Ratom,[1,2])
    velx=np.random.uniform(0,np.sqrt(2))
    vely=np.sqrt(2-velx*velx)
    vel[i,:]=[velx,vely]
    
    for j in range(0,i):
        if(np.sqrt(np.dot(pos[i,:]-pos[j,:],pos[i,:]-pos[j,:]))<2*Ratom):
            i=-1
            break
    if(i==Natoms-1):
        break
#################################################################################


## List indexing all pairs of atoms.
pairs = [[i,j] for i in range(Natoms) for j in range(i+1, Natoms)] 

## Positions and velocities as numpy arrays.
positions = np.array(pos)
velocities = np.array(vel)

##### Creates a picture of the atoms in the box ########
create_picture(positions)
pylab.show()
print(positions,"Initial config")



wall_times = [wall_time(positions[i][j], velocities[i][j], Ratom) for i in range(Natoms) for j in range(2)] 
# Pair collision times
pair_times = [pair_time(positions[i], velocities[i], positions[j], velocities[j], Ratom) for i,j in pairs] 
# The next collision event is the minimum of wall and pair collision times.
time = min(wall_times + pair_times)  

t = 0.0 # Initial time.
n_events =100# Number of collision events.


############### Event Loop #############################################    
for event in range(n_events):
    # Wall collision times for all atoms and their velocity components.
    wall_times = [wall_time(positions[i][j], velocities[i][j], Ratom) for i in range(Natoms) for j in range(2)] 
    # Pair collision times
    pair_times = [pair_time(positions[i], velocities[i], positions[j], velocities[j], Ratom) for i,j in pairs] 
    # The next collision event is the minimum of wall and pair collision times.
    next_event = min(wall_times + pair_times)  
    t += next_event 
    
    for i in range(Natoms):
        positions[i] += velocities[i]*next_event # Evolve positions to collision event
    if min(wall_times) < min(pair_times): # Check if next event is a collision with a wall
        wall_index = wall_times.index(next_event)
        particle, component = wall_index/2, wall_index%2
        velocities[particle][component] *= -1.0 ## Velocity component normal to wall changes sign
    else:
        pair_index = pair_times.index(next_event)
        particle_1, particle_2 = pairs[pair_index] # Indices of atoms participating in collision.
        rel_pos = positions[particle_2] - positions[particle_1]
        rel_vel = velocities[particle_2] - velocities[particle_1]
        distance = np.sqrt(np.dot(rel_pos,rel_pos))
        unit_perp = rel_pos/distance
        scal_prod = np.dot(rel_vel,unit_perp)
        velocities[particle_1] += scal_prod*unit_perp # Change in velocities of atoms colliding with each other
        velocities[particle_2] -= scal_prod*unit_perp
    

##################################################################################

##### Creates a picture of the atoms in the box ########
create_picture(positions)
pylab.show()
print(positions,"forward nth event")

'''
velocities *= -1.0

t = 0.0 # Initial time.

############### Event Loop ##############################################
for event in range(n_events):
    # Wall collision times for all atoms and their velocity components.
    wall_times = [wall_time(positions[i][j], velocities[i][j], Ratom) for i in range(Natoms) for j in range(2)] 
    # Pair collision times
    pair_times = [pair_time(positions[i], velocities[i], positions[j], velocities[j], Ratom) for i,j in pairs] 
    # The next collision event is the minimum of wall and pair collision times.
    next_event = min(wall_times + pair_times)  
    t += next_event 
    
    for i in range(Natoms):
        positions[i] += velocities[i]*next_event # Evolve positions to collision event
    if min(wall_times) < min(pair_times): # Check if next event is a collision with a wall
        wall_index = wall_times.index(next_event)
        particle, component = wall_index/2, wall_index%2
        velocities[particle][component] *= -1.0 ## Velocity component normal to wall changes sign
    else:
        pair_index = pair_times.index(next_event)
        particle_1, particle_2 = pairs[pair_index] # Indices of atoms participating in collision.
        rel_pos = positions[particle_2] - positions[particle_1]
        rel_vel = velocities[particle_2] - velocities[particle_1]
        distance = np.sqrt(np.dot(rel_pos,rel_pos))
        unit_perp = rel_pos/distance
        scal_prod = np.dot(rel_vel,unit_perp)
        velocities[particle_1] += scal_prod*unit_perp # Change in velocities of atoms colliding with each other
        velocities[particle_2] -= scal_prod*unit_perp

##################################################################################

positions += velocities*time # Evolve positions to collision event

##### Creates a picture of the atoms in the box ########
create_picture(positions)
pylab.show()

print(positions,"reversed time")
'''