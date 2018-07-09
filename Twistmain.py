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

folder = r'G:\Klaas\Tweezers\Doxo Project\2018_06_14_Doxo_4kb_constrained\CorrectedDat\No Doxo'
MinFitForce = 0.3 #Force cut-off for the WLC fit in pN
MaxFitForce = 6 #pN

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
        MagnetPosition = Tools.motor_pos(Filename[:-4]+'.log', 'Shift Start (mm)')
        RotationForce = Tools.calc_force(MagnetPosition)
        X,Y,Z_rot, Rot, T = Tools.get_rotation_data(data, headers)
        Sigma = Rot/(Pars['L_bp']/Pars['Pitch_nm'])
        Rotationfile = Filename
    else: 
        X,Y,Z_ext, F, T = Tools.get_force_data(data, headers)
        FEfile = Filename
 
    if Rotationfile and FEfile and len(Z_rot.T)==len(Z_ext.T):   #only analyse when both FE and Twist are from the same FOV
        for i,row in enumerate(Z_rot.T[:]):
            Pars = Tools.default_pars() 
            if np.max(F)<MinFitForce:
                continue
            Z_fit = (Z_ext[:,i])[np.all([F>MinFitForce,F<MaxFitForce], axis=0)]
            F_fit = F[np.all([F>MinFitForce,F<MaxFitForce], axis=0)]
            try: 
                z0fit, z0fitcov = curve_fit(lambda f, z0: Tools.offset_fit(f,z0,Pars), F_fit, Z_fit, p0=Pars['z0_nm'], bounds=(-3000,10000))
                Pars['z0_nm'] = z0fit[0]
                popt,pcov = curve_fit(lambda f, P: Tools.ewlc_fit(f,P,Pars['z0_nm'],Pars['S_pN'] ,Pars), F_fit, Z_fit, p0=[Pars['P_nm']], bounds=([1],[100]))
                Pars['P_nm']=popt[0]
                popt_S,pcov_S = curve_fit(lambda f, S: Tools.ewlc_fit(f,Pars['P_nm'],S,Pars['z0_nm'],Pars), F_fit, Z_fit, p0=[Pars['S_pN']], bounds=([10],[10000]))
                Pars['S_pN']=popt_S[0]
                popt_S,pcov_S = curve_fit(lambda f,P,S,z0: Tools.ewlc_fit(f,P,S,z0,Pars), F_fit, Z_fit, p0=[Pars['P_nm'],Pars['S_pN'],Pars['z0_nm']], bounds=([1,10,-3000],[100,10000,10000]))
                Pars['P_nm']=popt_S[0]
                Pars['S_pN']=popt_S[1]
                Pars['z0_nm'] = popt_S[2]
            except RuntimeError:
                print('>>>>>>>> fit error! <<<<<<<')
                continue
            if Pars['P_nm'] < 10*np.sqrt(pcov_S[0,0]) or Pars['P_nm'] > 95 or Pars['P_nm'] < 5:
                continue
            print('Offset = ',Pars['z0_nm'], 'P = ', Pars['P_nm'])
            fig1 = plt.figure()
            fig1.suptitle('Bead ' + str(i), y=.99)  
            ax1 = fig1.add_subplot(1, 2, 1)
            ax1.set_title('Twist ' + Rotationfile[:-4])
            ax1.set_ylabel(r'Extension')
            ax1.set_xlabel(r'Sigma')
            ax1.set_ylim(-200, Pars['L_bp']*Pars['dsDNA_nm_bp']*1.2)
            ax1.scatter(Sigma, row - np.min(row), color='blue',  s=5)
            ax1.text(np.min(Sigma), Pars['L_bp']*Pars['dsDNA_nm_bp']*1.1, 'F = '+str(round(RotationForce,2))+' pN', fontsize=12, verticalalignment='center', horizontalalignment='left')
            ax2 = fig1.add_subplot(1, 2, 2)
            ax2.set_title('Force Extension ' + FEfile[:-4])
            ax2.set_ylabel('Force')
            ax2.set_xlabel('Extension')
            ax2.scatter(Z_ext.T[i,:]-Pars['z0_nm'], F, color='red',  s=5)
            ax2.plot(Tools.wlc(F, Pars),F, color='black')  #plots the WLC
            ax2.text(np.min(Z_ext.T[i,:]-Pars['z0_nm']), np.max(F), 'P = '+str(round(Pars['P_nm'],1))+'±'+ str(round(np.sqrt(pcov_S[0,0]),1))+' nm', fontsize=12, verticalalignment='center', horizontalalignment='left')
            #ax2.text(np.min(Z_ext.T[i,:]), np.max(F)/10*9, 'stdv = ' + str(round(np.sqrt(pcov[0,0]),1)), fontsize=12, verticalalignment='center', horizontalalignment='left')
            fig1.savefig(Filename[:-4]+'_'+str(i)+'.png')
        #Rotationfile, FEfile = False, False
