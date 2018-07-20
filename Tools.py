# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 16:36:12 2018
@author: nhermans
"""
import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt
import json
import pandas as pd

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
    Z = data[:,headers.index('Z0'+' (um)')::4] * 1000
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
    T = np.array(data[:,headers.index('Time (s)')])
    shift = np.array(data[:,headers.index('Stepper shift (mm)')])
    Z = np.array(data[:,headers.index('Z0'+' (um)')::4] *1000)
    X = np.array(data[:,headers.index('X0'+' (um)')::4])
    Y = np.array(data[:,headers.index('Y0'+' (um)')::4])
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
    Pars = {}
    Pars['Filename'] = ''
    Pars['L_bp']= 4140#/3*4
    Pars['P_nm'] = 50
    Pars['S_pN'] = 1200
    Pars['dsDNA_nm_bp']=0.34
    Pars['z0_nm'] = 0
    Pars['Pitch_nm'] = 10.4
    Pars['MeltingTorque_pN_nm'] = 0.16
    Pars['Torque'] = 0
    Pars['kBT_pN_nm'] = 4.2 #pn/nm 
    Pars['MeasurementERR (nm)'] = 5     #tracking inaccuracy in nm
    Pars['date'] = ""  # date of measurement
    Pars['data'] = ""  # measurement
    Pars['bead'] = ""  # bead
    Pars['beads'] = ""  # total number of beads
    Pars['points'] = ""  # number of data-points
    Pars['points_exp'] = ""  # number of expected data-points
    Pars['points_frac'] = ""  # fraction of data-points
    Pars['X0_um'] = ""  # global X-position (um)
    Pars['Y0_um'] = ""  # global Y-position (um)
    Pars['Z0_um'] = ""  # global Z-position (um)
    Pars['dZ_um'] = ""  # absolute Z-extension after drift correction(um)
    Pars['Radius_nm'] = "" # Radius of xy rotation fit with a circle
    return Pars

def wlc(force,Pars): #in nm/pN, as fraction of L
    """Calculates WLC in nm/pN.
    Returns Z_WLC in nm"""
    f = np.array(force)
    return Pars['L_bp']*Pars['dsDNA_nm_bp']*(1 - 0.5*(np.sqrt(Pars['kBT_pN_nm']/(f*Pars['P_nm'])))+(f/Pars['S_pN']))

def ewlc_fit(f, P, S, z0, Pars):
    return Pars['L_bp']*Pars['dsDNA_nm_bp']*(1 - 0.5*(np.sqrt(Pars['kBT_pN_nm']/(f*P)))+(f/S)) + z0

def offset_fit(f, z0, Pars):
    return Pars['L_bp']*Pars['dsDNA_nm_bp']*(1 - 0.5*(np.sqrt(Pars['kBT_pN_nm']/(f*Pars['P_nm'])))+(f/Pars['S_pN'])) + z0

def fit_wlc(F,Z,Pars, MinFitForce, MaxFitForce):
    
    from scipy.optimize import curve_fit
    
    Z_fit = (Z)[np.all([F>MinFitForce,F<MaxFitForce], axis=0)]
    F_fit = F[np.all([F>MinFitForce,F<MaxFitForce], axis=0)]

    z0fit, z0fitcov = curve_fit(lambda f, z0: offset_fit(f,z0,Pars), F_fit[F_fit>10], Z_fit[F_fit>10], p0=Pars['z0_nm'], bounds=(-3000,10000))
    Pars['z0_nm'] = z0fit[0]
    popt,pcov = curve_fit(lambda f, P: ewlc_fit(f,P,Pars['z0_nm'],Pars['S_pN'] ,Pars), F_fit, Z_fit, p0=[Pars['P_nm']], bounds=([1],[100]))
    Pars['P_nm']=popt[0]
    popt_S,pcov_S = curve_fit(lambda f, S: ewlc_fit(f,Pars['P_nm'],S,Pars['z0_nm'],Pars), F_fit[F_fit>10], Z_fit[F_fit>10], p0=[Pars['S_pN']], bounds=([10],[2000]))
    Pars['S_pN']=popt_S[0]
    popt,pcov = curve_fit(lambda f,P,S,z0: ewlc_fit(f,P,S,z0,Pars), F_fit, Z_fit, p0=[Pars['P_nm'],Pars['S_pN'],Pars['z0_nm']], bounds=([1,10,-3000],[500,2000,10000]))
    Pars['P_nm']=popt[0]
    Pars['S_pN']=popt[1]
    Pars['z0_nm'] = popt[2]

    return popt, pcov
    

def read_log(LogFile, MotorName):
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
    
def fit_circle(X,Y):
    # coordinates of the barycenter
    x_m = np.mean(X)
    y_m = np.mean(Y)
    
    # calculation of the reduced coordinates
    u = X - x_m
    v = Y - y_m
    
    # linear system defining the center (uc, vc) in reduced coordinates:
    #    Suu * uc +  Suv * vc = (Suuu + Suvv)/2
    #    Suv * uc +  Svv * vc = (Suuv + Svvv)/2
    Suv  = sum(u*v)
    Suu  = sum(u**2)
    Svv  = sum(v**2)
    Suuv = sum(u**2 * v)
    Suvv = sum(u * v**2)
    Suuu = sum(u**3)
    Svvv = sum(v**3)
    
    # Solving the linear system
    A = np.array([ [ Suu, Suv ], [Suv, Svv]])
    B = np.array([ Suuu + Suvv, Svvv + Suuv ])/2.0
    uc, vc = np.linalg.solve(A, B)
    
    xc_1 = x_m + uc
    yc_1 = y_m + vc

    # Calcul des distances au centre (xc_1, yc_1)
    Ri_1     = np.sqrt((X-xc_1)**2 + (Y-yc_1)**2)
    R_1      = np.mean(Ri_1)
    #residu_1 = sum((Ri_1-R_1)**2)
    return R_1

def save_df(FileName, Pars):
    values = []
    for n, dict in enumerate(Pars):
        if n == 0:
            keys = list(dict.keys())
        values.append(list(dict.values()))

    df = pd.DataFrame(data=values)
    df.columns = keys

    df.to_csv(FileName, sep='\t')
