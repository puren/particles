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
                

    '''
    pos2     = colliderPos
    vel2     = vector(0.0, 0.0, 0.0)
    radius_2=0.2
    '''
   
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

def jacobiRotate(A,R):
    #print A
    #print '****'
    #print 'jacobi'
    tol=1.0e-09
    aMax=1.0
    n = len(A)
    maxRot = 5*(n**2)       # Set limit on number of rotations
    p = identity(3)*1.0     # Initialize transformation matrix
    for i in range(maxRot): # Jacobi rotation loop
        '''
        for i in range(maxRot): # Jacobi rotation loop 
            aMax,k,l = maxElem(a)
            if aMax < tol: return diag(a)
            rotate(a,p,k,l)
        print 'Jacobi method did not converge'
        '''
        aMax,k,l = maxElem(A)
        #print A
        if aMax > tol:
            aDiff=A[k,k]-A[l,l]
            
            if abs(A[k,l]) < abs(aDiff)*1.0e-36:
                t = A[k,l]/aDiff
                
            else:
                d=(aDiff)/(2.0*A[k,l])
                t=1.0/(abs(d)+math.sqrt(d*d+1.0))
                if(d<0.0): t=-t
                
            c=1.0/math.sqrt(t*t+1)
            s=t*c
            #print t
            A[k,k]=A[k,k]+(t*A[k,l])
            A[l,l]=A[l,l]+(t*A[k,l])
            A[k,l]=A[l,k]=0.0
            #for i in range(2):
            
            Rkp=c*R[k,k]+s*R[k,l]
            Rkq=-s*R[k,k]+c*R[k,l]
            R[k,l]=Rkp
            R[k,k]=Rkq
            Rkp=c*R[l,k]+s*R[l,l]
            Rkq=-s*R[l,k]+c*R[l,l]
            R[l,l]=Rkp
            R[l,k]=Rkq
            '''
            print '---'
            print Rkp
            print Rkq
            print R
            print '---'
            '''
        else:
            #print R
            return R
           
    return R           

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
    
    '''
    print g[i]
    print pos[i]
    print ''
    '''

'''
def projectConstraints(pos, new_pos, noParticles):
    d=0.05
    isIn=0
    p12=vector(0.0,0.0,0.0)
    for i in range(0, noParticles):
        for j in range(i, noParticles):
    
            p12.x=pos[i*3]-pos[j*3]
            p12.y=pos[i*3+1]-pos[j*3+1]
            p12.z=pos[i*3+2]-pos[j*3+2]
            abs_p=p12.mag
            if abs_p-d<0 and abs_p>0:
                isIn=1
                #delta1=-(w1/(w1+w2))*(abs(p12)-d)*(p12/abs(p12))
                #delta2=(w2/(w1+w2))*(abs(p12)-d)*(p12/abs(p12))
                delta1=-0.5*(abs_p-d)*(p12.x/abs_p)
                delta2=0.5*(abs_p-d)*(p12.x/abs_p)
                new_pos[i*3]=pos[i*3]+delta1
                new_pos[j*3]=pos[j*3]+delta2

                #delta1=-(w1/(w1+w2))*(abs(p12)-d)*(p12/abs(p12))
                #delta2=(w2/(w1+w2))*(abs(p12)-d)*(p12/abs(p12))
                delta1=-0.5*(abs_p-d)*(p12.y/abs_p)
                delta2=0.5*(abs_p-d)*(p12.y/abs_p)
                new_pos[i*3+1]=pos[i*3+1]+delta1
                new_pos[j*3+1]=pos[j*3+1]+delta2

                #delta1=-(w1/(w1+w2))*(abs(p12)-d)*(p12/abs(p12))
                #delta2=(w2/(w1+w2))*(abs(p12)-d)*(p12/abs(p12))
                delta1=-0.5*(abs_p-d)*(p12.z/abs_p)
                delta2=0.5*(abs_p-d)*(p12.z/abs_p)
                new_pos[i*3+2]=pos[i*3+2]+delta1
                new_pos[j*3+2]=pos[j*3+2]+delta2


def jacobi(a,tol = 1.0e-9): # Jacobi method
 
    def maxElem(a): # Find largest off-diag. element a[k,l]
        n = len(a)
        aMax = 0.0
        for i in range(n-1):
            for j in range(i+1,n):
                if abs(a[i,j]) >= aMax:
                    aMax = abs(a[i,j])
                    k = i; l = j
        return aMax,k,l
 
    def rotate(a,p,k,l): # Rotate to make a[k,l] = 0
        n = len(a)
        aDiff = a[l,l] - a[k,k]
        if abs(a[k,l]) < abs(aDiff)*1.0e-36: t = a[k,l]/aDiff
        else:
            phi = aDiff/(2.0*a[k,l])
            t = 1.0/(abs(phi) + sqrt(phi**2 + 1.0))
            if phi < 0.0: t = -t
        c = 1.0/sqrt(t**2 + 1.0); s = t*c
        tau = s/(1.0 + c)
        temp = a[k,l]
        a[k,l] = 0.0
        a[k,k] = a[k,k] - t*temp
        a[l,l] = a[l,l] + t*temp
        for i in range(k):      # Case of i < k
            temp = a[i,k]
            a[i,k] = temp - s*(a[i,l] + tau*temp)
            a[i,l] = a[i,l] + s*(temp - tau*a[i,l])
        for i in range(k+1,l):  # Case of k < i < l
            temp = a[k,i]
            a[k,i] = temp - s*(a[i,l] + tau*a[k,i])
            a[i,l] = a[i,l] + s*(temp - tau*a[i,l])
        for i in range(l+1,n):  # Case of i > l
            temp = a[k,i]
            a[k,i] = temp - s*(a[l,i] + tau*temp)
            a[l,i] = a[l,i] + s*(temp - tau*a[l,i])
 
    n = len(a)
    maxRot = 5*(n**2)       # Set limit on number of rotations
    p = identity(n)*1.0     # Initialize transformation matrix
    for i in range(maxRot): # Jacobi rotation loop 
        aMax,k,l = maxElem(a)
        if aMax < tol: return diag(a)
        rotate(a,p,k,l)
    print 'Jacobi method did not converge'

def calcCenterOfMass(pos, n):
    cm=vector(0.0,0.0,0.0)
    for i in range(n):
        cm=cm+pos[i]
    cm=cm/n
    return cm

def multVec(p,q):
    p_arr=np.matrix([p.x, p.y, p.z])
    p_arr=np.transpose(p_arr)
    q_arr=[q.x, q.y, q.z]
    
    return p_arr*q_arr

def calcApq(pos_o, pos_d, cm_o, cm_d, n):
    A_pq=[[0.0 for i in range(3) ] for j in range(3)]
    for i in range(n):
        q=pos_o[i]-cm_o
        p=pos_d[i]-cm_d
        A_pq=A_pq+multVec(p,q)

    return A_pq

def shapeMatching(pos_o, pos_d, g_all, n):
    cm_o=calcCenterOfMass(pos_o,n)
    cm_d=calcCenterOfMass(pos_d,n)
    A_pq=calcApq(pos_o, pos_d, cm_o, cm_d, n)
    A_pq_T=np.transpose(A_pq)
    mult=dot(A_pq_T, A_pq)
    mult_diag=jacobi(mult,3)
    mult_diag=np.diag(mult_diag)
    S=np.sqrt(mult_diag)
    S_inv=np.linalg.inv(S)
    R=np.dot(A_pq, S_inv)

    arr_cm_o=[cm_o.x, cm_o.y, cm_o.z]
    for i in range(n):
        arr_x_o=[pos_o[i].x,pos_o[i].y,pos_o[i].z]
        arr_x_d=[pos_d[i].x,pos_d[i].y,pos_d[i].z]
        v1 = np.subtract(arr_x_o,arr_cm_o)
        v2 = np.dot(R, v1)
        g  = np.sum([v2,arr_x_d])
        g_all[i]=vector(g[0,0], g[0,1], g[0,2])
'''        
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
'''
n_pos=[vector(0.0,0.0,0.0) for i in range(numParticles)]
g_all=[vector(0.0,0.0,0.0) for i in range(numParticles)]
'''

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
    
    
    '''
    shapeMatching(org_pos, n_pos, g, numParticles)
    
    integrate_sm(m_hPos, n_pos, g, m_hVel, dt, numParticles)

    '''
    
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
 
 
