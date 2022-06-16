# -*- coding: utf-8 -*-
"""
Created on Fri Aug 03 18:28:43 2018
Edited on Tue Oct 15 12:37:14 2019

@author: agharag
"""

# from tkinter.filedialog import askopenfilename
from netCDF4 import Dataset
import numpy as np
import time
from scipy import spatial
from operator import itemgetter
import sys

''' read mesh and element conn '''

# mesh = input("Enter name of ADCIRC file (fort.63.nc): ")
mesh = "C:/Users/Thomas Thelen/OneDrive - North Carolina State " \
       "University/CarolinaBeach/SMS/20220610_NC9_Nov2021/ADCIRC/fort.63.nc"
# var = input('Enter name of variable (zeta): ')
var = 'zeta'

t0 = time.time()

ncf = Dataset(mesh, 'r')
elem = ncf.variables['element'][:]
xnode = ncf.variables['x'][:]
ynode = ncf.variables['y'][:]
var = ncf.variables[var][:]
N_ele = elem.shape[0]
nt = var.shape[0]
N_nod = len(xnode)
nodes = np.column_stack((xnode, ynode))

''' read 63 file (maxele, water level, Hs, Tp, Dir, WindUV) '''

## read station list
## text file format: StationName Latitude Longgitude
# station = input("Enter name of station file (Stations.txt): ")
station = "C:/Users/Thomas Thelen/OneDrive - North Carolina State University/CarolinaBeach/Model_Inputs/Stations.txt"
stat_list = []
stat_names = []
with open(station, 'r') as f_st:
    #f_st.readline()
    for line in f_st:
        m = line.split()
        stat_list.append([float(m[1]), float(m[2])])  # add lat and long to list
        stat_names.append(m[0])

''' find nearest triangles'''
t2 = time.time()

### nearest k nodes to each station in stat_list 
CKD = list(spatial.cKDTree(nodes).query(stat_list, k=4)[1])

triangles = []
for i in range(len(stat_list)):  # for each station
    tri = []
    ps = CKD[i]  # nearest k nodes
    for node in ps:
        tri.extend(list(np.where(elem == node + 1)[0]))  # find nearest elements/triangles
    triangles.append(set(tri))

t3 = time.time()
print(' ***** found nearby elements in ', round(t3 - t2, 3), 's')

## find exact triangle- and then interpolate and get data

GAMMA = []
ADJ = []
for i in range(len(stat_list)):
    GAMMA.append(0)
    ADJ.append(0) # we will assign these values later
    xp, yp = stat_list[i]
    print(stat_names[i])

    Gtot_min = 100 # initialize value for Gtot error check
    for j in triangles[i]:
        adjc = elem[j] - 1
        XX = [xp, xnode[adjc[0]], xnode[adjc[1]], xnode[adjc[2]]]
        YY = [yp, ynode[adjc[0]], ynode[adjc[1]], ynode[adjc[2]]]


        def Area(X, Y):
            return round(abs((X[1] * Y[2] - X[2] * Y[1]) - (X[0] * Y[2] - X[2] * Y[0]) + (X[0] * Y[1] - X[1] * Y[0])),
                         8)


        def Sub_Area(a, b, c):
            return Area(list(itemgetter(a, b, c)(XX)), list(itemgetter(a, b, c)(YY)))


        Tot_A = Area(XX[1:], YY[1:])
        G0 = Sub_Area(0, 2, 3) / Tot_A
        G1 = Sub_Area(1, 0, 3) / Tot_A
        G2 = Sub_Area(1, 2, 0) / Tot_A
        Gtot = G0 + G1 + G2
        Error_Gtot = abs(Gtot - 1)
        Error_Gtot_min = abs(Gtot_min - 1)

        if Error_Gtot < Error_Gtot_min:
            Gtot_min = Gtot
            print('Sum of G0 + G1 + G2 = {}'.format(Gtot_min))
            GAMMA[i] = [G0, G1, G2]
            ADJ[i] = adjc


t4 = time.time()
print(' ***** found exact element and weight function in ', round(t4 - t3, 3), 's')
print("Gamma weighting array is {}".format(GAMMA))

Stat_val = []
## output data
for s in range(len(stat_list)):
    G0, G1, G2 = GAMMA[s]
    n0, n1, n2 = ADJ[s]
    S = []
    print(stat_names[s])
    for t in range(nt):
        mval = G0 * var[t][n0].data + G1 * var[t][n1].data + G2 * var[t][n2].data #.data unmasks var vals
        mval = mval/(G0 + G1 + G2) # adjust for weighting coefficients where G0+G1+G2 not exactly equal to zero
        print(t, s, mval)
        S.append(mval)
    Stat_val.append(S)

t5 = time.time()
print(' ***** output ', round(t5 - t4, 3), 's')
