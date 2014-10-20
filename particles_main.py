import math
import random
from visual import *
import time, sys
import numpy as np


'''
TODO:
1. Turn everything to np array
2. Parallelization
3. Optimize loops
4. add shape matching algorithm for rigid body(müller's paper)
5. elasticity simulation(müller's paper)
'''

def integrate(dt, pos, vel, index, gravity, globalDamping, boundaryDamping, radius_s):

    
    vel[index] = vel[index]+(gravity*dt)
    vel[index] = vel[index]*globalDamping
    pos[index]=pos[index]+vel[index]*dt
    

    if pos[index].x>1.0-radius_s:
        pos[index].x=1.0-radius_s
        vel[index].x=vel[index].x*boundaryDamping
    if pos[index].x<-1.0+radius_s:
        pos[index].x=-1.0+radius_s
        vel[index].x=vel[index].x*boundaryDamping

    if pos[index].y>1.0-radius_s:
        pos[index].y=1.0-radius_s
        vel[index].y=vel[index].y*boundaryDamping
    if pos[index].y<-1.0+radius_s:
        pos[index].y=-1.0+radius_s
        vel[index].y=vel[index].y*boundaryDamping

    if pos[index].z>1.0-radius_s:
        pos[index].z=1.0-radius_s
        vel[index].z=vel[index].z*boundaryDamping
    if pos[index].z<-1.0+radius_s:
        pos[index].z=-1.0+radius_s
        vel[index].z=vel[index].z*boundaryDamping
    

def updateGrid(pos, worldOrigin, cellSize, gridSize, gridCells, maxParticlesPerCell, indx):
    pos_v=pos[indx]
    gridPose=calcGridPos(pos_v, worldOrigin, cellSize)
    addParticleToCell(gridPose, gridSize, gridCells, maxParticlesPerCell, indx)

def collide(old_pos, old_vel, worldOrigin, cellSize, gridSize,
            maxParticlesPerCell, radius, spring, damping, shear,
            attraction, colliderPos, i, gridCounters, dt):

    pos=old_pos[i]
    vel=old_vel[i]
    o_pos=pos
    o_vel=vel
    
    gridPose = calcGridPos(pos, worldOrigin, cellSize)

    force    = vector(0.0, 0.0, 0.0)

    for z  in range(-1, 2):
        for y in range(-1, 2):
            for x in range(-1, 2):
                force = force + collideCell(gridPose+vector(x, y, z), i, pos, vel,
                                             old_pos, old_vel, radius,spring, damping, shear,
                                             attraction, gridCounters, gridCells)
                

   
    n_vel[i] = vel+force
    

def collideCell(gridPose, index, pos, vel, old_pos, old_vel, radius,
                 spring, damping, shear, attraction, gridCounters, gridCell_):

    force = vector(0.0, 0.0, 0.0)

    if ( (gridPose.x<0) or (gridPose.x > gridSize.x-1) or
         (gridPose.y<0) or (gridPose.y > gridSize.y-1) or
         (gridPose.z<0) or (gridPose.z > gridSize.z-1) ):
        return force
                        
    gridHash = calcGridHash(gridPose, gridSize)
    
    particlesInCell = gridCounters[gridHash]
    particlesInCell = min(particlesInCell, (maxParticlesPerCell-1))
    
    for j in range(0, particlesInCell) :
        index2 = gridCell_[gridHash*maxParticlesPerCell+j]
        if index2 != index :
            
            
            pos2 = old_pos[index2]
            vel2 = old_vel[index2]

            projVec = collideSpheres(pos, vel, pos2, vel2, radius, radius, spring, damping, shear, attraction)
            force   = force + projVec
    
    return force

def collideSpheres(pos, vel, pos2, vel2, radius1, radius2, spring, damping, shear, attraction):
    force_cs    = vector(0.0,0.0,0.0)
    relPos=vector(0.0,0.0,0.0)
    relPos      = pos2-pos
 
    dist        = relPos.mag
 
    collideDist = radius1+radius2
    
           
    if dist < collideDist :
        
        norm      = relPos/dist
        
        relVel    = vel2-vel
        tanVel    = relVel-(np.dot(relVel, norm)*norm)

        #spring force
        force_cs  = -spring *(collideDist-dist)*norm
        
        
        #damping
        force_d   = damping*relVel
        frorce_cs = force_cs+force_d
        
        #shear force
        force_cs  = force_cs + shear*tanVel
        
        #attraction
        force_cs  = force_cs + attraction*relPos
        
        #force_cs  = force_cs + (alpha*(g-pos)/dt)
                   
    return force_cs

def calcGridPos(pos, worldOrigin, cellSize):
    gridPose   = vector(0, 0, 0)
    gridPose.x = math.floor((pos.x-worldOrigin.x)/cellSize.x)
    gridPose.y = math.floor((pos.y-worldOrigin.y)/cellSize.y)
    gridPose.z = math.floor((pos.z-worldOrigin.z)/cellSize.z)
    
    return gridPose

def calcGridHash(gridPose, gridSize):
    gridPose.x = max(0, min(gridPose.x, (gridSize.x-1)))
    gridPose.y = max(0, min(gridPose.y, (gridSize.y-1)))
    gridPose.z = max(0, min(gridPose.z, (gridSize.z-1)))

    gridHash   = int(((gridPose.z*gridSize.y)*gridSize.x)+(gridPose.y*gridSize.x)+gridPose.x)
    
    return gridHash

def addParticleToCell(gridPose, gridSize, gridCells, masParticlesPerCell, indx):
    gridHash              = calcGridHash(gridPose, gridSize)
    counter               = gridCounter[gridHash]
    gridCounter[gridHash] = gridCounter[gridHash]+1
     
    counter          = min(counter, (maxParticlesPerCell-1))
    indx2            = int(gridHash*maxParticlesPerCell+counter)
    gridCells[indx2] = indx

def calcCenterOfMass(pos, n):
    cm=np.zeros(3)
    for i in range(n):
        cm=cm+pos[i]
    cm=cm/n
    return cm

def calcApq(Apq, Aqq, pos_o, pos_d, cm_o, cm_d, n):
    
    for i in range(n):
        q=pos_o[i]-cm_o
        p=pos_d[i]-cm_d
    
        for j in range(3):
            for k in range(3):
                Apq[j,k]=Apq[j,k]+(p[j]*q[k])
                Aqq[j,k]=Apq[j,k]+(p[j]*q[k])

def maxElem(a): # Find largest off-diag. element a[k,l]
    n = len(a)
    aMax = 0.0
    k=0
    l=0
    for i in range(n-1):
        for j in range(i+1,n):
            #print a[i,j]
            #print ''
            if abs(a[i,j]) >= aMax:
                aMax = abs(a[i,j])
                k = i; l = j
    return aMax,k,l

    

def polarDecomposition(A, R, S):
    R=np.identity(3)
    AT=np.transpose(A)
    ATA=np.dot(AT,A)
   
    U=np.identity(3)
    jacobiRotate(ATA, U)

    #ATA^-1
    l0=ATA[0,0]
    if(l0<=0.0): l0=0.0
    else: l0=1.0/math.sqrt(l0)
    l1=ATA[1,1]
    if(l1<=0.0): l1=0.0
    else: l1=1.0/math.sqrt(l1)
    l2=ATA[2,2]
    if(l2<=0.0): l2=0.0
    else: l2=1.0/math.sqrt(l2)

    S1=np.zeros((3,3))
    S1[0,0]=l0*U[0,0]*U[0,0] + l1*U[0,1]*U[0,1] + l2*U[0,2]*U[0,2]
    S1[0,1]=l0*U[0,0]*U[1,0] + l1*U[0,1]*U[1,1] + l2*U[0,2]*U[1,2]
    S1[0,2]=l0*U[0,0]*U[2,0] + l1*U[0,1]*U[2,1] + l2*U[0,2]*U[2,2]
    S1[1,0]=S1[0,1]
    S1[1,1]=l0*U[1,0]*U[1,0] + l1*U[1,1]*U[1,1] + l2*U[1,2]*U[1,2]
    S1[1,2]=S1[2,1]=l0*U[1,0]*U[2,0] + l1*U[1,1]*U[2,1] + l2*U[1,2]*U[2,2]
    S1[2,0]=S1[0,2]
    S1[2,2]=l0*U[2,0]*U[2,0] + l1*U[2,1]*U[2,1] + l2*U[2,2]*U[2,2]
 
    #print S1[2,1]

    R=np.dot(A,S1)
    RT=np.transpose(R)
    S=np.dot(RT, A)
    #print R

    return R
    
    

def shapeMatching(o_pos, pos, g, n):
    o_pos_np=np.zeros((n,3))
    pos_np=np.zeros((n,3))
    for i in range(n):
        o_pos_np[i]=np.array([o_pos[i].x, o_pos[i].y, o_pos[i].z])
        pos_np[i]=np.array([pos[i].x, pos[i].y, pos[i].z])
    
    cm=calcCenterOfMass(pos_np, n)
    cm_o=calcCenterOfMass(o_pos_np, n)
    
    Apq=np.zeros((3,3))
    Aqq=np.zeros((3,3))
    
    calcApq(Apq, Aqq, o_pos_np, pos_np, cm_o, cm, n)
    
    R=0
    S=np.zeros((3,3))

    R=polarDecomposition(Apq, R, S)
    print R
    
    for i in range(n):
        p=o_pos_np[i]-cm_o
        Rx=np.zeros(3)+cm

        Rx[0]=R[0,0]*p[0]+R[0,1]*p[1]+R[0,2]*p[2]
        Rx[1]=R[1,0]*p[0]+R[1,1]*p[1]+R[1,2]*p[2]
        Rx[2]=R[2,0]*p[0]+R[2,1]*p[1]+R[2,2]*p[2]
        g_t=Rx+cm
        
        #print ''
        g[i]=vector(g_t[0], g_t[1], g_t[2])
        #print g[i]

def integrate_sm(pos, n_pos, g, vel, dt, n):
    dt1=1.0/dt
    
    for i in range(n):
        vel[i]=vel[i]+(g[i]-pos[i])*dt1
        pos[i]=n_pos[i]
    
    
display(title='Examples of particles',
        x=0, y=0, width=900, height=700,
        center=(0,0,0), background=(0,0,0))

#Initialize particle system
gsize=64
radius_s=1.0 / gsize
numParticles=9
gridSize=vector(gsize, gsize, gsize)
spacing = radius_s*2.0
jitter=radius_s*0.01
dt=0.5
globalDamping=1.0
gravity=vector(0.0, -0.0003, 0.0)
boundaryDamping=-0.5
numGridCells=int(gridSize.x*gridSize.y*gridSize.z)
maxParticlesPerCell=4
worldOrigin=vector(-1, -1, -1)
cellSize=vector(radius_s*2, radius_s*2, radius_s*2)
attraction=0.0
spring=0.5
damping=0.02
shear=0.1
colliderPos = vector(-2.0, -2.0, -2.0)
alpha=1.0

m_hPos=[vector(0.0,0.0,0.0) for i in range(numParticles)]
org_pos=[vector(0.0,0.0,0.0) for i in range(numParticles)]
m_hVel=[vector(0.0,0.0,0.0) for i in range(numParticles)]
n_vel=[vector(0.0,0.0,0.0) for i in range(numParticles)]
n_pos=[vector(0.0,0.0,0.0) for i in range(numParticles)]
org_pos=[vector(0.0,0.0,0.0) for i in range(numParticles)]
g=[vector(0.0,0.0,0.0) for i in range(numParticles)]

#Create the wall
side = 1.0
thk = 0.01
s2 = 2*side - thk
s3 = 2*side + thk
wallR = box (pos=( side, 0, 0), size=(thk, s2, s3),  color = color.red)
wallL = box (pos=(-side, 0, 0), size=(thk, s2, s3),  color = color.red)
wallB = box (pos=(0, -side, 0), size=(s3, thk, s3),  color = color.blue)
wallT = box (pos=(0,  side, 0), size=(s3, thk, s3),  color = color.blue)
wallBK = box(pos=(0, 0, -side), size=(s2, s2, thk), color = (0.7,0.7,0.7))

#Create the sphere
ball = []
s=int(math.ceil(pow(numParticles, (1.0/3.0))))
gridObj=vector(s, s, s)
for z in range(0, int(gridObj.z)):
    for y in range(0, int(gridObj.y)):
        for x in range(0, int(gridObj.x)):
            i = int((z*gridObj.y*gridObj.x) + (y*gridObj.x) + x)
            if i < numParticles:
                tmp=vector(x,y,z)
                m_hPos[i].x   = (x*spacing) + radius_s - 1.0 + (random.random()*2.0-1.0)*jitter
                m_hPos[i].y   = (y*spacing) + radius_s - 1.0 + (random.random()*2.0-1.0)*jitter;
                m_hPos[i].z   = (z*spacing) + radius_s - 1.0 + (random.random()*2.0-1.0)*jitter;
                                           
                m_hVel[i]   = vector(0.0,0.0,0.0);
 
                ball             = ball + [sphere (color = color.green, radius = radius_s)]
                ball[i].pos      = m_hPos[i]
                ball[i].mass     = 1.0
                ball[i].velocity = vector(0.0, 0.0, 0.0)

org_pos=m_hPos
iter =0
while 1:
    
    sleep(0.2)
    gridCells=[0 for i in range(numGridCells*maxParticlesPerCell)]
    gridCounter=[0 for i in range(0,numGridCells)]

    
    for i in range (numParticles):
        integrate(dt, m_hPos, m_hVel, i, gravity, globalDamping, boundaryDamping, radius_s)
    
    
    for i in range (numParticles):
        ball[i].pos=m_hPos[i]
        ball[i].velocity=m_hVel[i]

    
    for i in range (numParticles):
        updateGrid(m_hPos, worldOrigin, cellSize, gridSize, gridCells, maxParticlesPerCell, i)

    for i in range (numParticles):
        collide(m_hPos, m_hVel, worldOrigin, cellSize, gridSize,
                maxParticlesPerCell, radius_s, spring, damping,
                shear, attraction, colliderPos, i, gridCounter,dt)
        
    for i in range (numParticles):
        m_hVel[i]   = n_vel[i]
    
    

    iter = iter + 1
 
 
