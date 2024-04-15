import numpy as np
import matplotlib.pyplot as plt
from decimal import *
getcontext().prec = 3

mag=np.load('nfmms1bx.15.npy')                #magnetic field (x-component, mms1) 
time=np.load('Time.24.02.2019.mms1.npy')      #time

#for mms1
l1 = time.size

##################################################################################
# lambda values
lmda = np.logspace(-9,3,num=1000)
lmd = lmda.size

##################################################################################
#defining constants
a= np.zeros((3,5))
for i in range(3):
    a[i,i]=1
    a[i,i+1]=-2
    a[i,i+2]=1
a_t = np.transpose(a)
b=np.zeros((3,3))
for i in range(3):
    b[i,i]=4
for i in range(0,2):
    b[i,i+1]=1
    b[i+1,i]=1
diff1 = time1[l1-1]-time1[0]
h = diff1/(l1-1)
A = a/h
A_t = a_t/h
B = (h/6)*b

T = np.zeros([5,5])
for i in range(0,4):
    T[i,i+1]= 1
for i in range(1,4):
    T[i,i] = 4
T[0,0]=2
T[4,4]=2
for i in range(0,4):
    T[i+1,i] = 1
T_inv = np.linalg.inv(T)
x = np.zeros([5,5])
for i in range(0,4):
    x[i,i+1]= 1
for i in range(1,4):
    x[i,i] = 0
x[0,0]=-1
x[4,4]=1
for i in range(0,4):
    x[i+1,i] = -1
D = (3/h)*x
################################################################
#function for calculating exact smoothing parameter (lambda)
def closestx(spx, c):
    return (min(range(len(spx)), key = lambda y: abs(spx[y]-c)))

def lmdsx(z):
    f_arr = np.zeros(lmd)
    d_arr = np.zeros(lmd)
    P_arr = np.zeros(lmd)
    E_arr = np.zeros(lmd)

    for j in range(lmd):
        const = (lmda[j]**2)/h
        k = B+(const*(np.dot(A,A_t)))
        k_inv = np.linalg.inv(k)
        d = mag[z-2:z+3]
        c = np.dot(A,d)
        M = np.dot(k_inv,c)
        M_t = np.transpose(M)
        f = d - (const)*(np.dot(A_t,M))
        v = np.dot(B,M)
        intg = np.dot(M_t,v)
        P = np.sqrt(intg/(4*h))
        E = np.sqrt(np.sum((f-d)**2)/4)
        P_arr[j]=P
        E_arr[j]=E
    c = 0.95
    f = []
    spx = np.zeros(lmd)
    for j in range(lmd):
        ssp=np.divide(P_arr[j],P_arr[0])
        spx[j]=ssp
    p=closestx(spx, c)
    d_lmd= lmda[p]

    return d_lmd

#############################################
#calculating the derivative

fx1_arr = np.zeros(l1)
dx1_arr = np.zeros(l1)
mx1_arr = np.zeros(l1)


for i in range(2,l1-2):
    print(i-2)
    Lambdax = lmdsx(i)
    print("after calling smoothing",Lambdax)
    constx = (Lambdax**2)/h
    kx = B+(constx*(np.dot(A,A_t)))
    kx_inv = np.linalg.inv(kx)
    dx = mag[i-2:i+3]
    cx = np.dot(A,dx)
    Mx = np.dot(kx_inv,cx)
    fx = dx - (constx)*(np.dot(A_t,Mx))
    gx = np.dot(D,fx)
    mx = np.dot(T_inv,gx)
    fx1_arr[i]=fx[2]
    dx1_arr[i]=dx[2]
    mx1_arr[i]=mx[2]

np.save('nfderivativebx.15.mms1',mx1_arr)
np.save('nfsmooth.magneticbx.15.mms1',fx1_arr)

