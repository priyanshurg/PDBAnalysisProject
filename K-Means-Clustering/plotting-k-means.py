#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 02:58:44 2023

@author: pgupta
"""
numcluster = 4 # Edit the number here between 3-6 to see plots of the k-means ...
                # ... clustering based on coordination number data 
                
                
                
import pickle
with open('xyzstore.pkl', 'rb') as f:  
    xyz2store, xyz3store, xyz4store, xyz5store, xyz6store = pickle.load(f)
    
    
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans

# Load the 3D coordinate data
Str = "xyz"+str(numcluster)+'store'
Vars = vars()
load_data = Vars[Str]
data = np.array(load_data)
kmeans = KMeans(n_clusters=numcluster, random_state=0)


kmeans.fit(data)
cluster_centers = kmeans.cluster_centers_
labels = kmeans.labels_


# Plot the 3D data points and cluster centers
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')

ox = 0
oy = 0
oz = 0

# Plot data points with color based on cluster labels
ax.scatter(data[:, 0], data[:, 1], data[:, 2], c=labels, cmap='viridis')

# Plot cluster centers with red crosses
ax.scatter(cluster_centers[:, 0], cluster_centers[:, 1], cluster_centers[:, 2], marker='x', s=200, c='red', label='Centroids')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.scatter(ox,oy,oz,c = 'black', label = 'Origin')

plt.legend()
plt.title("K-Means Clustering in 3D")


ax2.scatter(cluster_centers[:, 0], cluster_centers[:, 1], cluster_centers[:, 2], marker='x', s=200, c='red', label='Centroids')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')
ax2.scatter(ox,oy,oz,c = 'black', label = 'Origin')


plt.legend()
plt.title("K-Means Clustering in 3D")
plt.show()
