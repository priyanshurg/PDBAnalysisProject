#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 05:38:33 2023

@author: pgupta
"""
# %matplotlib widget 
from myfunctions import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
from random import *
import scipy as sp
from scipy.linalg import orthogonal_procrustes
import numpy as np
test= 0


'''
This is an implementation of a Partial Procrustes analysis where the shape of all the other 
sites is compared to that of the first file in the array. The sites have been partioned based
on the number of residues in them.

Unlike regular procrustes analysis, this implementation skips the scaling step during the optimisation.

'''

def partial_procrustes_method(A, B, scale):
    assert len(A) == len(B)

    N = A.shape[0];  # total points

    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)

    # center the points
    AA = A - np.tile(centroid_A, (N, 1))
    BB = B - np.tile(centroid_B, (N, 1))

    # dot is matrix multiplication for array
    if scale:
        H = np.transpose(BB) * AA / N
    else:
        H = np.transpose(BB) * AA

    U, S, Vt = np.linalg.svd(H)

    R = Vt.T * U.T


    if scale:
        varA = np.var(A, axis=0).sum()
        c = 1 / (1 / varA * np.sum(S))  # scale factor
        t = -R * (centroid_B.T * c) + centroid_A.T
    else:
        c = 1
        t = -R * centroid_B.T + centroid_A.T

    return c, R, t

scaling = False
current_directory = os.getcwd()

if test==0:
    pdb_files = ["AllData/"+filename for filename in os.listdir(current_directory+"/AllData/") if filename.endswith(".pdb")]
else:
    pdb_files = ["TestData/"+filename for filename in os.listdir(current_directory+"/TestData/") if filename.endswith(".pdb")]


coordination_number = []
pdb_names = []
i=1

# f0 = plt.figure("Procrustes_readings")
f1 = plt.figure("Coordination Number = 2")
f2 = plt.figure("Coordination Number = 3")
f3 = plt.figure("Coordination Number = 4")
f4 = plt.figure("Coordination Number = 5")
f5 = plt.figure("Coordination Number = 6")
# f6 = plt.figure("empty")
# ax0 = Axes3D(f0)
ax1 = Axes3D(f1)
ax2 = Axes3D(f2)
ax3 = Axes3D(f3)
ax4 = Axes3D(f4)
ax5 = Axes3D(f5)
# ax6 = Axes3D(f6)


i2 = 0
i3=0
i4=0
i5=0
i6=0
proscrutes2 = []
proscrutes3 = []
proscrutes4 = []
proscrutes5 = []
proscrutes6 = []

p2ref = []
p3ref = []
p4ref = []
p5ref = []
p6ref = []

xstore =[]
pr2label = []
pr3label = []
pr4label = []
pr5label = []
pr6label = []


for filename in pdb_files:
    aliasname = "structure"
    p53_1cpd = get_file(filename,aliasname)
    nh4s = []
    nh4atoms = []

    for residue in p53_1cpd.get_residues():
        if residue.resname == 'NH4':
            nh4s.append(residue)
            for atom in residue.get_atoms():
                nh4atoms.append(atom)
    for eachatom in nh4atoms:
        distance = 3.75
# =============================================================================
#         closestatoms = get_closest_atoms(p53_1cpd, eachatom, distance)
# =============================================================================
        closestatoms = get_closest_atoms_phi_sort(p53_1cpd, eachatom, distance)
        closestatoms = dict(sorted(closestatoms.items(), key=lambda item: item[1]))
        #print(filename,"atom index = ",eachatom.get_full_id()[3][1],eachatom.coord,len(closestatoms))
        atomsph_com, r_com,theta_com,phi_com = get_spherical_coords(closestatoms,eachatom)        
        atomcart_com, x_com,y_com,z_com = get_cartesian_coords(closestatoms,eachatom)
        coordination_number.append(len(atomcart_com))
        pdb_names.append(filename)
        #Fig(1) for coordination Number 2, Fig(2) for coordination Number 3, Fig(4) for Coordination Number 5, 
        # Fig(5) for coordination number 6
        if len(atomcart_com) == 2:
            if i2 ==0:
                for i in range(0,len(atomcart_com)):
                    p2ref.append([x_com[i],y_com[i],z_com[i]])
            pt2 = []
            for i in range(len(atomcart_com)):
                pt2.append([x_com[i],y_com[i],z_com[i]])

            A = np.matrix(p2ref)
            B = np.matrix(pt2)
            n = B.shape[0]
            s, ret_R, ret_t = partial_procrustes_method(A, B, scaling)
            #s, ret_R, ret_t = umeyama(A, B)

            # Find the error
            B2 = (ret_R * B.T) + np.tile(ret_t, (1, n))
            b = B2.T

            xplot = [b[i,0] for i in range(0,len(b))]
            yplot = [b[i,1] for i in range(0,len(b))]
            zplot = [b[i,2] for i in range(0,len(b))]
            # x0_com,y0_com,z0_com = [0,0,0]
            ax1.plot(xplot,yplot,zplot,"-o", color=[random(),random(),random()],label = filename+str(eachatom.get_full_id()[3][1]))
            # ax1.scatter(x0_com,y0_com,z0_com,color="black")
            i2+=1
        elif len(atomcart_com) == 3:
            if i3 ==0:
                for i in range(0,len(atomcart_com)):
                    p3ref.append([x_com[i],y_com[i],z_com[i]])
            pt3 = []
            for i in range(len(atomcart_com)):
                pt3.append([x_com[i],y_com[i],z_com[i]])
           #b = sp.spatial.procrustes(p3ref,pt3)

            A = np.matrix(p3ref)
            B = np.matrix(pt3)
            n = B.shape[0]
            s, ret_R, ret_t = partial_procrustes_method(A, B, scaling)
            #s, ret_R, ret_t = umeyama(A, B)

            # Find the error
            B2 = (ret_R * B.T) + np.tile(ret_t, (1, n))
            b = B2.T
            
            xplot = [b[i,0] for i in range(0,len(b))]
            yplot = [b[i,1] for i in range(0,len(b))]
            zplot = [b[i,2] for i in range(0,len(b))]
            # x0_com,y0_com,z0_com = [0,0,0]
            # xstore.append([x_com,y_com,z_com,filename])
            ax2.plot(xplot,yplot,zplot,"-o", color=[random(),random(),random()],label = filename+str(eachatom.get_full_id()[3][1]))
            # ax2.scatter(x0_com,y0_com,z0_com,color="black")
            i3+=1
        elif len(atomcart_com) == 4:
            if i4 ==0:
                for i in range(0,len(atomcart_com)):
                    p4ref.append([x_com[i],y_com[i],z_com[i]])
            pt4 = []
            for i in range(len(atomcart_com)):
                pt4.append([x_com[i],y_com[i],z_com[i]])

            A = np.matrix(p4ref)
            B = np.matrix(pt4)
            n = B.shape[0]
            s, ret_R, ret_t = partial_procrustes_method(A, B, scaling)
            #s, ret_R, ret_t = umeyama(A, B)

            # Find the error
            B2 = (ret_R * B.T) + np.tile(ret_t, (1, n))
            b = B2.T
            
            xplot = [b[i,0] for i in range(0,len(b))]
            yplot = [b[i,1] for i in range(0,len(b))]
            zplot = [b[i,2] for i in range(0,len(b))]
            
            # x0_com,y0_com,z0_com = [0,0,0]
            ax3.plot(xplot,yplot,zplot,"-o", color=[random(),random(),random()],label = filename+str(eachatom.get_full_id()[3][1]))
            # ax3.scatter(x0_com,y0_com,z0_com,color="black")
            i4+=1
        elif len(atomcart_com) == 5:
            if i5 ==0:
                for i in range(0,len(atomcart_com)):
                    p5ref.append([x_com[i],y_com[i],z_com[i]])
            pt5 = []
            for i in range(len(atomcart_com)):
                pt5.append([x_com[i],y_com[i],z_com[i]])

            A = np.matrix(p5ref)
            B = np.matrix(pt5)
            n = B.shape[0]
            s, ret_R, ret_t = partial_procrustes_method(A, B, scaling)
            #s, ret_R, ret_t = umeyama(A, B)

            # Find the error
            B2 = (ret_R * B.T) + np.tile(ret_t, (1, n))
            b = B2.T
            
            xplot = [b[i,0] for i in range(0,len(b))]
            yplot = [b[i,1] for i in range(0,len(b))]
            zplot = [b[i,2] for i in range(0,len(b))]
            # x0_com,y0_com,z0_com = [0,0,0]
            ax4.plot(xplot,yplot,zplot,"-o", color=[random(),random(),random()],label = filename+str(eachatom.get_full_id()[3][1]))
            # ax4.scatter(x0_com,y0_com,z0_com,color="black")         
            i5+=1
        elif len(atomcart_com) == 6:
            if i6 ==0:
                for i in range(0,len(atomcart_com)):
                    p6ref.append([x_com[i],y_com[i],z_com[i]])
            pt6 = []
            for i in range(len(atomcart_com)):
                pt6.append([x_com[i],y_com[i],z_com[i]])

            A = np.matrix(p6ref)
            B = np.matrix(pt6)
            n = B.shape[0]
            s, ret_R, ret_t = partial_procrustes_method(A, B, scaling)
            #s, ret_R, ret_t = umeyama(A, B)

            # Find the error
            B2 = (ret_R * B.T) + np.tile(ret_t, (1, n))
            b = B2.T
            
            xplot = [b[i,0] for i in range(0,len(b))]
            yplot = [b[i,1] for i in range(0,len(b))]
            zplot = [b[i,2] for i in range(0,len(b))]
            # x0_com,y0_com,z0_com = [0,0,0]
            ax5.plot(xplot,yplot,zplot,"-o", color=[random(),random(),random()],label = filename+str(eachatom.get_full_id()[3][1]))
            # ax5.scatter(x0_com,y0_com,z0_com,color="black")         
            i6+=1
        i+=1
ax1.legend(prop={'size': 2})
ax2.legend(prop={'size': 2})
ax3.legend(prop={'size': 2})
ax4.legend(prop={'size': 2})
ax5.legend(prop={'size': 2})

plt.show()

import pickle

pickle.dump(f1, open('Plots/CoordinationNum2.fig.pickle', 'wb'))
pickle.dump(f2, open('Plots/CoordinationNum3.fig.pickle', 'wb'))
pickle.dump(f3, open('Plots/CoordinationNum4.fig.pickle', 'wb'))
pickle.dump(f4, open('Plots/CoordinationNum5.fig.pickle', 'wb'))
pickle.dump(f5, open('Plots/CoordinationNum6.fig.pickle', 'wb'))