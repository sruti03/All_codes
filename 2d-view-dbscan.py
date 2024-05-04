#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from itertools import chain

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


# In[12]:


import numpy as np
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt


dbscan = DBSCAN(eps=2.5, min_samples=8)  # Adjust eps and min_samples as needed
dbscan.fit(blobs)
cluster_labels = dbscan.labels_

num_clusters = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
print("Number of clusters:", num_clusters)

# Plot the data points and centroids
fig, ax = plt.subplots()
#plt.figure(figsize=(4, 8))
ci=ax.imshow(J, cmap ='rainbow',interpolation ='nearest', origin ='lower')
ci=plt.colorbar(ci)

for cluster_label in set(cluster_labels):
    if cluster_label == -1:
        plt.scatter(blobs[cluster_labels == cluster_label][:, 1], blobs[cluster_labels == cluster_label][:, 0],
                    color='black', label='Noise',s=2)
    else:
        plt.scatter(blobs[cluster_labels == cluster_label][:, 1], blobs[cluster_labels == cluster_label][:, 0],
                    label=f'Cluster {cluster_label}',s=2)
        
core_sample_indices = dbscan.core_sample_indices_
core_points = blobs[core_sample_indices]

is_core = np.isin(np.arange(len(blobs)), core_sample_indices)
unique_clusters = np.unique(cluster_labels)        
boundary_points_by_cluster = {}
for cluster in unique_clusters:
    boundary_points = blobs[(cluster_labels == cluster) & (~is_core)]
    boundary_points_by_cluster[cluster] = boundary_points
        
for cluster in unique_clusters:
    boundary_points = boundary_points_by_cluster[cluster]
    plt.scatter(boundary_points[:, 1], boundary_points[:, 0], edgecolor='brown', color='brown', marker='x', label=f'Boundary Cluster {cluster}',s=2)

plt.title('K-means Clustering')
plt.xlabel(r"X($d_i$)",
            fontweight ='bold',
            size=10)

plt.ylabel(r"Y($d_i$)",
           fontweight ='bold',
           size=10)
#plt.grid(True)
#plt.xlim(0,50)
#plt.ylim(0,50)
plt.show()











