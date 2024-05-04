#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from itertools import chain
import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors

J=np.load('nfJtz.npy')
l,m=J.shape

Jm=J

for s in range(10):
    Js= np.zeros_like(J)
    for i in range(1,l-1):
        for j in range(1,m-1):
            #print(J[i-1:i+1,j-1:j+1])
            Js[i,j]=np.average((0.25*np.sum(J[i-1:i+1,j-1:j+1]))+0.75*J[i,j])
    J=Js

Ly=25
Jsum = np.sum(J)/(l*m)
Javg = np.average(J)
Jrms = np.sqrt(np.average(J*J))
print(l,m)
Jth=2*Jrms

k=0
I_th=[]
J_th=[]
Jk= np.zeros_like(J)
J_abs=abs(J)

indx,indy = np.where(J_abs >= Jth)

a=indx.size
for i in range(a):
    ix=indx[i]
    iy=indy[i]
    Jk[ix,iy]=1
    
Jk_abs=abs(Jk)
            
#print(Jk)
    
blobs = np.column_stack((indx, indy))

min_samples = 8
K = max(2, min_samples)  
nn = NearestNeighbors(n_neighbors=K)
nn.fit(blobs)
distances, _ = nn.kneighbors(blobs)

kth_distances = distances[:, -1] 
sorted_distances = np.sort(kth_distances)[::-1] 

plt.figure(figsize=(10, 6))
plt.plot(range(1, len(sorted_distances) + 1), sorted_distances, 'bo-', linewidth=2)
plt.xlabel("Points")
plt.ylabel("Distance to K-th Nearest Neighbor")
plt.title("Elbow Method to Determine Optimal Epsilon for DBSCAN")
plt.xlim(0,6000)
plt.ylim(0,10)
plt.show()


# In[ ]:




