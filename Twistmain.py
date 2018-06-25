# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 16:07:05 2018

@author: nhermans

Script to plot Force Extension curves next to twist curves from the same data set

How to use:
1) Optional: correct drift with "subtract stuck bead script"
2) Select folder with .dat files, make sure the folder includes the corresponding .log files
3) Run script
4) Select the curves that make sense for futher analysis
    
"""

import os 
import Tools
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

folder = r'C:\Users\Klaas\Documents\Doxo Project\2018_06_14_Doxo_4kb_constrained\CorrectedDat'
filenames = os.listdir(folder)
os.chdir(folder)
Pars = Tools.default_pars()    
Rotationfile = False
FEfile = False
        
for Filenum, Filename in enumerate(filenames):
    if Filename[-4:] != '.dat':
        continue
    data, headers = Tools.read_dat(Filename)
    
    if 'Stepper rot (turns)' in headers: 
        X,Y,Z_rot, Rot, T = Tools.get_rotation_data(data, headers)
        Rotationfile = True
        MagnetPosition = Tools.motor_pos(Filename[:-4]+'.log', 'Shift Start (mm)')
        RotationForce = Tools.calc_force(MagnetPosition)
        Sigma = Rot/(Pars['L_bp']/Pars['Pitch_nm'])
    else: 
        X,Y,Z_ext, F, T = Tools.get_force_data(data, headers)
        FEfile = True
 
    if Rotationfile and FEfile and len(Z_rot.T)==len(Z_ext.T):   #only analyse when both FE and Twist are from the same FOV
        for i,row in enumerate(Z_rot.T[:]):
            Pars = Tools.default_pars() 
            Z_fit = np.array(Z_ext.T[i,:][F>0.5]*1000)
            F_fit = np.array(F[F>0.5])
            try: 
                z0fit, z0fitcov = curve_fit(lambda f, z0: Tools.offset_fit(f,z0,Pars), F_fit, Z_fit, p0=[Pars['z0_nm']])
                Pars['z0_nm'] = z0fit[0]
                popt,pcov = curve_fit(lambda f, P, S: Tools.wlc_fit(f,P,S,Pars), F_fit, Z_fit, p0=[Pars['P_nm'],Pars['S_pN']], bounds=[[1,100],[800,1600]])
            except RuntimeError:
                print('>>>>>>>> fit error! <<<<<<<')
                continue
            Pars['P_nm']=popt[0]
            Pars['S_pN']=popt[1]
            print('Offset = ',Pars['z0_nm'], 'P = ', Pars['P_nm'], 'S = ',Pars['S_pN'])
            fig1 = plt.figure()
            fig1.suptitle(Filename + '  bead ' + str(i), y=.99)  
            ax1 = fig1.add_subplot(1, 2, 1)
            ax1.set_title('Twist')
            ax1.set_ylabel(r'Extension')
            ax1.set_xlabel(r'Sigma')
            ax1.scatter(Sigma, row, color='blue',  s=5)
            ax1.text(np.min(Sigma), np.max(row), 'F = '+str(round(RotationForce,2))+' pN', fontsize=12, verticalalignment='center', horizontalalignment='left')
            ax2 = fig1.add_subplot(1, 2, 2)
            ax2.set_title('Force Extension')
            ax2.set_ylabel('Force')
            ax2.set_xlabel('Extension')
            ax2.scatter(Z_ext.T[i,:], F, color='red',  s=5)
            ax2.plot(Tools.wlc(F_fit, Pars)/1000,F_fit, color='black')  #plots the WLC
            ax2.text(np.min(Z_fit)/1000, np.max(F), 'P = ' + str(round(Pars['P_nm'],1)) + ' nm', fontsize=12, verticalalignment='center', horizontalalignment='right')
            fig1.show()
        Rotationfile, FEfile = False, False
