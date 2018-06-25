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
#    T=T[Rot != 0]    
#    X=X[Rot != 0]
#    Y=Y[Rot != 0]
#    Z=Z[Rot != 0]
#    Rot=Rot[Rot != 0]
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
    par['P_nm'] = 50
    par['S_pN'] = 1200
    par['dsDNA_nm_bp']=0.34
    par['z0_nm'] = 0
    par['Pitch_nm'] = 10.4
    par['MeltingTorque_pN_nm'] = 0.16
    par['Torque'] = 0
    par['kBT_pN_nm'] = 4.2 #pn/nm 
    par['MeasurementERR (nm)'] = 5     #tracking inaccuracy in nm
    return par

def wlc(force,Pars): #in nm/pN, as fraction of L
    """Calculates WLC in nm/pN.
    Returns Z_WLC in nm"""
    f = np.array(force)
    return Pars['L_bp']*Pars['dsDNA_nm_bp']*(1 - 0.5*(np.sqrt(Pars['kBT_pN_nm']/(f*Pars['P_nm'])))+(f/Pars['S_pN'])) + Pars['z0_nm']

def wlc_fit(f, P, S, Pars):
    return Pars['L_bp']*Pars['dsDNA_nm_bp']*(1 - 0.5*(np.sqrt(Pars['kBT_pN_nm']/(f*P))+(f/S))) + Pars['z0_nm']

def offset_fit(f, z0, Pars):
    return Pars['L_bp']*Pars['dsDNA_nm_bp']*(1 - 0.5*(np.sqrt(Pars['kBT_pN_nm']/(f*Pars['P_nm'])))+(f/Pars['S_pN'])) + z0

def motor_pos(LogFile, MotorName):
    """Open the corresponding .log files from magnetic tweezers. Returns False if the file is not found"""
    try: 
        f = open(LogFile, 'r')
    except FileNotFoundError: 
        print(LogFile, '========> No valid logfile found')
        return 0   
    content = f.readlines()
    f.close()
    for lines in content:
        P =lines.split(' = ')
        if P[0]==MotorName:
            return round(float(P[1].strip('\n')),1)
    print("<<<<<<<<<<", MotorName, "not found >>>>>>>>>>")
    return 0

def plot_sigma(Sigma,Z,fig1):
    ax1 = fig1.add_subplot(1, 2, 1)
    ax1.set_title('Twist')
    ax1.set_ylabel(r'Extension')
    ax1.set_xlabel(r'Sigma')
    ax1.scatter(Sigma, Z, color='black', lw=0.1, s=5)

def plot_force(F,Z,fig1):
    ax2 = fig1.add_subplot(1, 2, 2)
    ax2.set_title('Force Extension')
    ax2.set_ylabel(r'Extension')
    ax2.set_xlabel(r'Force')
    ax2.scatter(F, Z, color='black', lw=0.1, s=5)
