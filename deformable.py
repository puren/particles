import math
import random
import time, sys
import numpy as np
import matplotlib.pyplot as plt


def invert(A):
    d = A[0,0]*A[1,1] - A[0,1]*A[1,0]
    if d == 0.0:
        return False
    d = 1.0/d

    d00 = A[1,1]*d
    d01 = -A[0,1]*d
    d10 = -A[1,0]*d
    d11 = A[0,0]*d

    A[0,0] = d00
    A[0,1] = d01
    A[1,0] = d10
    A[1,1] = d11

    return True

def multVec(A, v):
    res = np.zeros(2)
    
    res[0] = A[0,0] * v[0] + A[0,1] * v[1]
    res[1] = A[1,0] * v[0] + A[1,1] * v[1]
    
    return res

def multiply(A, B):
    Res = np.zeros((2, 2))
    
    Res[0,0] = A[0,0] * B[0,0] + A[0,1]*B[1,0]
    Res[0,1] = A[0,0] * B[0,1] + A[0,1]*B[1,1]
    Res[1,0] = A[1,0] * B[0,0] + A[1,1]*B[1,0]
    Res[1,1] = A[1,0] * B[0,1] + A[1,1]*B[1,1]

    return Res

def multiplyTransposedLeft(left, right):
    res = np.zeros((2,2))
    
    res[0,0] = left[0,0] * right[0,0] + left[1,0]*right[1,0]
    res[0,1] = left[0,0] * right[0,1] + left[1,0]*right[1,1]
    res[1,0] = left[0,1] * right[0,0] + left[1,1]*right[1,0]
    res[1,1] = left[0,1] * right[0,1] + left[1,1]*right[1,1]

    return res

def jacobiRotate(A, R):
    d = (A[0,0] - A[1,1])/(2.0*A[0,1])
    t = 1.0/(abs(d)+math.sqrt(d*d + 1.0))

    if (d < 0.0): t = -t

    c = 1.0 / math.sqrt(t*t + 1)
    s = t*c

    A[0,0] = A[0,0] + t*A[0,1]
    A[1,1] = A[1,1] - t*A[0,1]
    A[0,1] = A[1,0] = 0.0

    for k in range(0,2):
        Rkp = c*R[k,0] + s*R[k,1]
        Rkq = -s*R[k,0] + c*R[k,1]
        
        R[k,0] = Rkp
        R[k,1] = Rkq
    
    return R, A

def eigenDecomposition(A, R):
    
    R = np.identity(2)
    R, A = jacobiRotate(A, R)

    return R, A

def polarDecomposition(A, R, S):
    
    R = np.identity(2)

    ATA = multiplyTransposedLeft(A,A)
    U   = np.zeros((2,2))
    U, ATA = eigenDecomposition(ATA, U)

    l0 = ATA[0,0]
    if (l0 <= 0.0): l0 = 0.0
    else: l0 = 1.0/math.sqrt(l0)
    l1 = ATA[1,1]
    if (l1 <= 0.0): l1 = 0.0
    else: l1 = 1.0/math.sqrt(l1)

    S1 = np.zeros((2,2))
    S1[0,0] = l0*U[0,0]*U[0,0] + l1*U[0,1]*U[0,1]
    S1[0,1] = l0*U[0,0]*U[1,0] + l1*U[0,1]*U[1,1]
    S1[1,0] = S1[0,1]
    S1[1,1] = l0*U[1,0]*U[1,0] + l1*U[1,1]*U[1,1]

    R = multiply(A, S1)
    S = multiplyTransposedLeft(R, A)

    return R, S

#set deformability parameters
dt                 = 0.01
gravity            = np.array([0.0, -9.8])
bounds_min         = np.array([-1.0, -1.0])
bounds_max         = np.array([1.0, 1.0])
alpha              = 0.1
beta               = 0.0
quadraticMatch     = False
volumeConservation = False
allowFlip          = True
alpha              = 1.0
rigidMatch              = True

#set particle parameters
gSize         = 64
radius        = 1.0 / gSize
mNumParticles = 9
jitter        = radius * 0.01
spacing       = radius * 2.0

mOriginalPos = np.zeros((mNumParticles,2))
mPos         = np.zeros((mNumParticles,2))
mNewPos      = np.zeros((mNumParticles,2))
mGoalPos     = np.zeros((mNumParticles,2))
mMasses      = np.ones(mNumParticles)
mVel         = np.zeros((mNumParticles,2))
mFixed       = [False for i in range(0, mNumParticles)]

s       = int(math.ceil(pow(mNumParticles, (1.0/3.0))))
gridObj = np.array([s, s, s])
for y in range(0, int(gridObj[1])):
    for x in range(0, int(gridObj[0])):
        i = int((y*gridObj[1]) + x)
        if i < mNumParticles:
            mOriginalPos[i,0] = (x*spacing) + radius - 0.5 + (random.random()*1.0-0.5)*jitter
            mOriginalPos[i,1] = (y*spacing) + radius - 0.5 + (random.random()*1.0-0.5)*jitter;
            mPos[i]     = mOriginalPos[i]
            mNewPos[i]  = mOriginalPos[i]
            mGoalPos[i] = mOriginalPos[i]

#draw the points
plt.ion()
line, = plt.plot(mOriginalPos[:, 0], mOriginalPos[:, 1], 'ro')
plt.ylim([-1,1])
plt.xlim([-1,1])
plt.show()
plt.draw()
time.sleep(1)

iter = 0
while 1:
    
    print 'iter: {0}'.format(iter)
    
    #external forces
    #print mPos
    for i in range(0, mNumParticles):
        if mFixed[i] == True:continue
        mVel[i]    += gravity*dt
        mNewPos[i]  = mPos[i] + mVel[i]*dt
        mGoalPos[i] = mOriginalPos[i]

    restitution = 0.9
    for i in range(0, mNumParticles):
        if mFixed[i] == True: continue
        if mNewPos[i,0] < bounds_min[0] or mNewPos[i,0] > bounds_max[0]:
            mNewPos[i,0] = mPos[i,0] - mVel[i,0]*dt*restitution
            mNewPos[i,1] = mPos[i,1]
        if mNewPos[i,1] < bounds_min[1] or mNewPos[i,1] > bounds_max[1]:
            mNewPos[i,1] = mPos[i,1] - mVel[i,1]*dt*restitution
            mNewPos[i,0] = mPos[i,0]

        #clamp
        if mNewPos[i,0] > bounds_max[0]: mNewPos[i,0] = bounds_max[0]
        if mNewPos[i,1] > bounds_max[1]: mNewPos[i,1] = bounds_max[1]
        if mNewPos[i,0] < bounds_min[0]: mNewPos[i,0] = bounds_min[0]
        if mNewPos[i,1] < bounds_min[1]: mNewPos[i,1] = bounds_min[1]

    #project positions
    cm         = np.zeros(2)
    originalCm = np.zeros(2)
    mass       = 0.0

    for i in range(mNumParticles):
        m = mMasses[i]
        
        if mFixed[i]: m = m * 100.0

        mass       = mass + m
        cm         = cm + mNewPos[i]*m
        originalCm = originalCm + mOriginalPos[i]*m

    cm         = cm/mass
    originalCm = originalCm/mass

    Apq = np.zeros((2, 2))
    Aqq = np.zeros((2, 2))

    for i in range(mNumParticles):
        p = mNewPos[i]-cm
        q = mOriginalPos[i]-originalCm
        m = mMasses[i]
        Apq[0, 0] = Apq[0, 0] + m * p[0] * q[0]
        Apq[0, 1] = Apq[0, 1] + m * p[0] * q[1]
        Apq[1, 0] = Apq[1, 0] + m * p[1] * q[0]
        Apq[1, 1] = Apq[1, 1] + m * p[1] * q[1]

        Aqq[0, 0] = Aqq[0, 0] + m * q[0] * q[0]
        Aqq[0, 1] = Aqq[0, 1] + m * q[0] * q[1]
        Aqq[1, 0] = Aqq[1, 0] + m * q[1] * q[0]
        Aqq[1, 1] = Aqq[1, 1] + m * q[1] * q[1]

    detApq = np.linalg.det(Apq)

    if allowFlip != True and detApq<0.0:
        Apq[0,1] = -Apq[0,1]
        Apq[1,1] = -Apq[1,1]

    R    = np.zeros((2,2))
    S    = np.zeros((2,2))
    R, S = polarDecomposition(Apq, R, S)

    if rigidMatch == True:
        A = Aqq
        invert(A)
        A = multiply(Apq, A)
        if volumeConservation == True:
            det = np.linalg.det(A)
            if det != 0.0:
                det = 1.0/math.sqrt(abs(det))
                if det > 2.0: det = 2.0
                A = A*det
        

        for i in range(mNumParticles):
            if mFixed[i] == True: continue
            q           = mOriginalPos[i] - originalCm
            mGoalPos[i] = multVec(R,q) + cm
            mNewPos[i]  = mNewPos[i] + (mGoalPos[i] - mNewPos[i]) * alpha
    elif quadraticMatch != True:
        A = Aqq
        invert(A)
        A = multiply(Apq, A)
        if volumeConservation == True:
            det = np.linalg.det(A)
            if det != 0.0:
                det = 1.0/math.sqrt(abs(det))
                if det > 2.0: det = 2.0
                A = A*det

        T = R * (1.0 - beta) + A * beta

        for i in range(mNumParticles):
            if mFixed[i] == True: continue
            q           = mOriginalPos[i] - originalCm
            mGoalPos[i] = multVec(T,q) + cm
            mNewPos[i]  = mNewPos[i] + (mGoalPos[i] - mNewPos[i]) * alpha

    #integrate
    dt1 = 1/dt
    for i in range(mNumParticles):
        mVel[i] = (mNewPos[i] - mPos[i]) * dt1
        mPos[i] = mNewPos[i]

    line.set_ydata(mPos[:, 1])
    line.set_xdata(mPos[:, 0])
    plt.draw()
    time.sleep(0.05)
    iter = iter + 1

print mPos


