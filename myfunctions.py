#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 26 16:50:38 2023

@author: pgupta
"""
from __future__ import print_function
import math
from Bio import PDB
import numpy as np
import matplotlib.pyplot as plt
from Bio.Data.IUPACData import protein_letters_3to1

def prot3to1(resname):
    string = protein_letters_3to1[resname.capitalize()]
    return string


def get_file(filename,aliasname):
    repository = PDB.PDBList()
    parser = PDB.PDBParser()
    
    structure = parser.get_structure(aliasname, filename)
    return structure


def get_closest_atoms(pdb_struct, ref_atom, distance):
    atoms = {}
    rx, ry, rz = ref_atom.coord
    for atom in pdb_struct.get_atoms():
        if atom == ref_atom:
            continue
        
        x, y, z = atom.coord        
        my_dist = math.sqrt((x - rx)**2 + (y - ry)**2 + (z - rz)**2)
        if my_dist < distance:
            atoms[atom] = my_dist
    return atoms

def get_closest_atoms_phi_sort(pdb_struct, ref_atom, distance):
    atoms = {}
    rx, ry, rz = ref_atom.coord
    for atom in pdb_struct.get_atoms():
        if atom == ref_atom:
            continue
        
        x, y, z = atom.coord        
        my_dist = math.sqrt((x - rx)**2 + (y - ry)**2 + (z - rz)**2)
        if my_dist < distance and atom.get_full_id()[3][0]==' ':
            r, theta, phi = get_r_theta_phi(atom,ref_atom)
            atoms[atom] = phi
    return atoms

def vecdist(atom1, atom2):
    r1x, r1y, r1z = atom1.coord
    r2x, r2y, r2z = atom2.coord
    my_dist = math.sqrt((r1x-r2x)**2+(r1y-r2y)**2+(r1z-r2z)**2)
    return my_dist

def veclength(x,y,z):
    my_dist = math.sqrt((x)**2+(y)**2+(z)**2)
    return my_dist

def get_r_theta_phi(atom1,ref_atom):
    xnew, ynew, znew = atom1.coord - ref_atom.coord
    r = veclength(xnew,ynew,znew)
    theta = np.arccos(znew/r)
    phi = np.sign(ynew)*np.arccos(xnew/math.sqrt(xnew**2+ynew**2))
    return r, theta, phi    

def get_spherical_coords(closestatoms,ref_atom):
    r_store=[]
    theta_store=[]
    phi_store= []
    atomlist =[]
    for each in closestatoms:
        if each.element == "O" or each.element =="N":
            r,theta,phi = get_r_theta_phi(each,ref_atom)
            r_store.append(r)
            theta_store.append(theta)
            phi_store.append(phi)
            atomlist.append(each)
    return atomlist,r_store,theta_store,phi_store
        
def get_cartesian_coords(closestatoms,ref_atom):
    xlist,ylist,zlist = [],[],[]
    atomlist =[]
    for each in closestatoms:
        if each.element == "O" or each.element == "N":
            x,y,z = each.coord - ref_atom.coord
            xlist.append(x)
            ylist.append(y)
            zlist.append(z)
            atomlist.append(each)
    return atomlist,xlist,ylist,zlist
            



  
def get_cartesian_coords_com(closestatoms,ref_atom):
    xlist,ylist,zlist = [],[],[]
    atomlist =[]
    for each in closestatoms:
        if each.element == "O" or each.element == "N":
            x,y,z = each.coord - ref_atom.coord
            xlist.append(x)
            ylist.append(y)
            zlist.append(z)
            atomlist.append(each)
    x_com = np.mean(xlist)
    y_com = np.mean(ylist)
    z_com = np.mean(zlist)
    xl = [x-x_com for x in xlist]     
    yl = [y-y_com for y in ylist]     
    zl = [z-z_com for z in zlist]     
            
    return atomlist,xl,yl,zl


def histshow(listname):
    plt.figure()    
    nc = len(set(listname))
    plt.hist(listname,rwidth = 0.9, bins = nc)
    plt.xticks(np.linspace(0, nc-1, 2*nc+1)[1::2])