# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 16:36:12 2018

@author: nhermans
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt

def read_dat(Filename):
    f = open(Filename, 'r')
    #get headers
    headers = f.readlines()[0]
    headers.rstrip()
    headers = headers.split('\t')
    f.close()  
    #get data
    data = genfromtxt(Filename, skip_header = 1)
    return data, headers

def get_rotation_data(data, headers):
    """Open .dat file from magnetic tweezers, averages the least moving beads and substracts them from the signal. 
    Output is a 2D array with all the data
    ***kwargs:
        Beads = number of beads to use for averaging, default = 3
        MedianFilter = LowPass filter for applied to averaged signal, default = 5. Needs to be an odd number
    """
    T = data[:,headers.index('Time (s)')]
    Rot = data[:,headers.index('Stepper rot (turns)')]
    Z = data[:,headers.index('Z0'+' (um)')::4]
    X = data[:,headers.index('X0'+' (um)')::4]
    Y = data[:,headers.index('Y0'+' (um)')::4]
    T=T[Rot != 0]    
    X=X[Rot != 0]
    Y=Y[Rot != 0]
    Z=Z[Rot != 0]
    Rot=Rot[Rot != 0]
    return X,Y,Z, Rot, T

def get_force_data(data, headers):
    """Open .dat file from magnetic tweezers, averages the least moving beads and substracts them from the signal. 
    Output is a 2D array with all the data
    ***kwargs:
        Beads = number of beads to use for averaging, default = 3
        MedianFilter = LowPass filter for applied to averaged signal, default = 5. Needs to be an odd number
    """
    T = data[:,headers.index('Time (s)')]
    shift = data[:,headers.index('Stepper shift (mm)')]
    Z = data[:,headers.index('Z0'+' (um)')::4]
    X = data[:,headers.index('X0'+' (um)')::4]
    Y = data[:,headers.index('Y0'+' (um)')::4]
    F = calc_force(shift)
    return X,Y,Z, F, T

# calculate force
def calc_force(i):
    A = 85  # for 2.8 um beads (pN)
    l1 = 1.4  # decay length 1 (mm)
    l2 = 0.8  # decay length 2 (mm)
    f0 = 0.01  # force-offset (pN)
    return A * (0.7 * np.exp(-i / l1) + 0.3 * np.exp(-i / l2)) + f0
    
def default_pars():
    """Default fitting parameters, returns a {dict} with 'key'= paramvalue"""
    par = {}
    par['L_bp']= 4139
    par['Pitch_nm'] = 10.4
    par['MeltingTorque_pN'] = 0.16
    par['Torque'] = 0
    par['kBT_pN_nm'] = 4.2 #pn/nm 
    par['MeasurementERR (nm)'] = 5     #tracking inaccuracy in nm
    return par

def plot_sigma(Sigma,Z,fig1):
    ax1 = fig1.add_subplot(1, 2, 1)
    ax1.set_title('Twist')
    ax1.set_ylabel(r'Extension')
    ax1.set_xlabel(r'Sigma')
    ax1.scatter(Sigma, Z, color='black', lw=0.1, s=5)

def plot_force(F,Z,fig1):
    ax1 = fig1.add_subplot(1, 2, 1)
    ax1.set_title('Force Extension')
    ax1.set_ylabel(r'Extension')
    ax1.set_xlabel(r'Force')
    ax1.scatter(F, Z, color='black', lw=0.1, s=5)