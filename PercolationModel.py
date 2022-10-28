# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 23:31:59 2021

@author: David Hayes (Python 3.8)

Description : This program models a percolating system on a square lattice and 
determines if a spanning cluster exists.
A Percolating system models the motion and filtering of fluids through porous 
materials, e.g. the filtration of water through soil, the strength of a porous 
network and the flow of oil through porous rock.
In defining the occurrence of a cluster, the following rules are applied. Two 
neighbouring sites are merged together to form one larger cluster, if they are 
both occupied. The cluster extends to all sites that neighbour another member 
of the cluster. A cluster which touches all four sides of the lattice is termed 
a spanning cluster. If a spanning cluster exists, a path will form, allowing 
fluid to percolate through the porous material.

Here, the sites in the lattice are occupied randomly, each time checking its 
neighbours and assigning a cluster number, until a spaning cluster is formed. 
The fraction of occupied sites in the lattice gives an estimate of the critical
concentration, p_c, for this particular percolation sequence.


"""


import numpy as np
import matplotlib.pyplot as plt
from random import randrange
#================================= def first =================================
def first(L,lat, clust):
    
    # Occupying a random site (First cluster)
    x = randrange(L) # Choosing random
    y = randrange(L) # x, y co-ords
   
    lat[x,y] = 1 # Site is occupied by cluster number 1
    clust.append(int(lat[x,y])) # Add cluster number to list
    
    return clust,lat
#================================= def second =================================
def second(L,lat, clust):
    # Picking a random site and checking its neighbours
    x = randrange(L)
    y = randrange(L)
    
    cn = neighbourhood(x,y,lat) # Check neighbourhood for clusters
    count = len(cn) # Count is the number of clusters in neighbourhood
    if count > 0:
        # label site as 1
        lat[x,y] = 1
    else:
        # label site as 2
        lat[x,y] = 2
        clust.append(int(lat[x,y]))
                    
    return clust, lat

#================================ def next_site ===============================
def next_site(L,lat, clust):
    # Choosing a random site and checking its neighbours
    x1 = randrange(L)
    y1 = randrange(L)
    
    cn = neighbourhood(x1,y1,lat) # Check neighbourhood for clusters
    clust, lat = label(x1,y1,lat,clust,cn) # Labelling site and relabelling lat
    
    return clust, lat

#================================ def next_site ===============================
def label(x1,y1,lat,clust,cn):
    # Labels the site lat[x1,y1]. Relabel it's neighbours and merges clusters
    # together if neccessary.
    
    count = len(cn) # Number of occupied neighbours
    if count == 0: # No neighbours are occupied     
           # label site with next cluster number
           lat[x1,y1] = np.amax(lat)+1
           clust.append(int(lat[x1,y1]))
           
    elif count == 1: # Only one cluster among all neighbouring occupied sites
        lat[x1,y1] = cn[0] # Set new site to this number
        
    elif count > 1: # More than one cluster number, bridges multiple cluster
        for i in range(-1,2):
            for j in range(-1,2):
                x, y = x1 + i, y1 + j
                
                if x != L and y != L and x!=-1 and y != -1:#Stop out of bounds
                    # If site is occupied, it needs to be relabelled
                    if lat[x,y] != 0 and lat[x,y] != min(cn):
                        if lat[x,y] in clust:
                            # Cluster number no longer exists
                            clust.remove(int(lat[x,y])) 
                        # Setting new cluster number to the smallest value
                        lat[x,y] = min(cn)        
        # Relabelling cluster numbers in the whole lattice
        for m in range(L):
            for n in range(L):
                if lat[m,n] in cn:
                    lat[m,n] = min(cn)
              
        lat[x1,y1] = min(cn)
        
    return clust, lat
#============================== def neighbourhood =============================
def neighbourhood(x1,y1,lat):
    # Find the occupied neighbours
    
    cn = [] #cluster numbers in neighbourhood of site
    for i in range(-1,2): # Checking neighbourhood
        for j in range(-1,2):
            x, y = x1 + i, y1 + j
            
            if x != L and y != L and x!=-1 and y != -1: # stop out of bounds
                if lat[x,y] != 0: # If site is occupied
                
                    if lat[x,y] not in cn: # If it is a new cluster number
                        cn.append(lat[x,y]) # add new cluster number to list
    return cn
   
#============================== def print_clusters ============================
def print_lattice(L,lat,clust):
  
   plt.figure(figsize = (6,6))
   plt.grid(color = 'grey',linestyle ='--')
   plt.xlim(-0.9,L), plt.ylim(-0.9,L)
   plt.xticks(range(0,L,1)),plt.yticks(range(0,L,1))
  
   c = ['b','k','r','c','g','purple','y', 'orange','grey','pink','brown',
        'beige','violet']
   
   for y in range(L):
       for x in range(L):
           i = int(lat[y,x])
           if i in clust :
                plt.plot(x,y,color = c[clust.index(i)], marker = 'o')
   
   
           
#================================= def edges ==================================
def spanning(L,lat):
    # Check if each edge has a common cluster number
    side1,side2,side3,side4 = [],[],[],[] # Store numbers on each side
    for i in range(0,L):
            side1.append(lat[0,i])
            side2.append(lat[i,L-1])
            side3.append(lat[L-1,i])
            side4.append(lat[i,0])
    # Storing common values between side 1 and side 2
    c1 = [value for value in side1 if value in side2]
    # Storing common values between side 3 and side 4
    c2 = [value for value in side3 if value in side4]
    # Checking for common elements
    c3 = list(set(c1).intersection(c2)) 
    
    for j in range(len(c3)):
        # if the common cluster number is > 0, all four sides are hit
        if c3[j]>0: 
           
            return True
    return False

#=================================== def pc ===================================
def pc(L,lat):
    # Calculates the fraction of occupied sites in the lattice.
    occupied = 0 # Used to count occupied sites
    size = L*L # size of the lattice
   
    for y in range(L):    
       for x in range(L):
           # If position in lattice is occupied, increase the count
           if lat[x,y] > 0: 
               occupied += 1
    
    return occupied/size
       
 #================================= def main ==================================   


                              
L = 20 # lattice dimensions (LxL) 
lat = np.zeros([L,L])# Creating an empty square lattice

clust = []  # Storing all cluster numbers

clust, lat = first(L,lat, clust) # Creating first cluster

clust,lat = second(L,lat, clust) # Occupying second site

while(spanning(L,lat)==False): # Continue with more sites until common
                               # cluster appears.
    clust, lat = next_site(L,lat, clust) # Occupying next site
  
print_lattice(L,lat,clust)
print('Fraction of occupied sites, Pc :',pc(L,lat))

