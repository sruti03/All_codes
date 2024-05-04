#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from sklearn.cluster import KMeans

#current density and its shape
J=np.load('nfJtz.npy')
l,m=J.shape

Jm=J

#smoothing the data
for s in range(10):
    Js= np.zeros_like(J)
    for i in range(1,l-1):
        for j in range(1,m-1):
            #print(J[i-1:i+1,j-1:j+1])
            Js[i,j]=np.average((0.25*np.sum(J[i-1:i+1,j-1:j+1]))+0.75*J[i,j])
    J=Js

#claculating the threshold
Jsum = np.sum(J)/(l*m)
Javg = np.average(J)
Jrms = np.sqrt(np.average(J*J))
print(l,m)
Jth=2*Jrms

J_abs=abs(J)

#points above threshold
indx,indy = np.where(J_abs >= Jth)

#creating array to visulaize regions above threshold
a=indx.size
Jk= np.zeros_like(J)
for i in range(a):
    ix=indx[i]
    iy=indy[i]
    Jk[ix,iy]=1
    
Jk_abs=abs(Jk)
            
#creating arrays for the clustering, with J values in z-axis    
#blobs = np.column_stack((indx, indy))
z=J[indx,indy]
Jblobs = np.column_stack((indx, indy,z))


# In[7]:


import numpy as np
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt


dbscan = DBSCAN(eps=2.5, min_samples=8)  # Adjust eps and min_samples as needed
dbscan.fit(Jblobs)
cluster_labels = dbscan.labels_

num_clusters = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
print("Number of clusters:", num_clusters)

# Plot the data points and centroids
#fig, ax = plt.subplots()
#plt.figure(figsize=(4, 8))
#ci=ax.imshow(J, cmap ='rainbow',interpolation ='nearest', origin ='lower')
#ci=plt.colorbar(ci)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

for cluster_label in set(cluster_labels):
    if cluster_label == -1:
        ax.scatter(Jblobs[cluster_labels == cluster_label][:, 1], Jblobs[cluster_labels == cluster_label][:, 0],Jblobs[cluster_labels == cluster_label][:, 2],
                    color='black', label='Noise',s=2)
    else:
        ax.scatter(Jblobs[cluster_labels == cluster_label][:, 1], Jblobs[cluster_labels == cluster_label][:, 0],Jblobs[cluster_labels == cluster_label][:, 2],
                    label=f'Cluster {cluster_label}',s=2)
        
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


# In[ ]:




