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
import numpy as np

folder = r'G:\Klaas\Tweezers\Doxo Project\2018_08_30_dimethylDoxo_167x25\2018_08_30_dimethylDoxo_167x2_CorrectedDat'
MinFitForce = 0.3 #Force cut-off for the WLC fit in pN
MaxFitForce = 65 #pN

filenames = os.listdir(folder)
os.chdir(folder)
Pars = Tools.default_pars()    
Rotationfile = False
FEfile = False
DataOut = [] #used to safe all results
DataOut.append(Tools.default_pars())
        
for Filenum, Filename in enumerate(filenames):
    if Filename[-4:] != '.dat':
        continue
    data, headers = Tools.read_dat(Filename)
    
    if 'Stepper rot (turns)' in headers: #Find out if file is FE or twist file
        RotationMotor = np.array(data[:,headers.index('Stepper rot (turns)')])
        if np.mean(RotationMotor) == 0:
            X_ext,Y_ext,Z_ext, F, T = Tools.get_force_data(data, headers)
            FEfile = Filename
            Framerate = Tools.read_log(Filename[:-4]+'.log', 'Framerate (fps)')
        else:
            MagnetPosition = Tools.read_log(Filename[:-4]+'.log', 'Shift Start (mm)')
            RotationForce = Tools.calc_force(MagnetPosition)
            X_rot,Y_rot,Z_rot, Rot, T = Tools.get_rotation_data(data, headers)
            Sigma = Rot/(Pars['L_bp']/Pars['Pitch_nm'])
            Rotationfile = Filename
            Framerate = Tools.read_log(Filename[:-4]+'.log', 'Framerate (fps)')
    else:
        X_ext,Y_ext,Z_ext, F, T = Tools.get_force_data(data, headers)
        FEfile = Filename
        Framerate = Tools.read_log(Filename[:-4]+'.log', 'Framerate (fps)')
        
    if Rotationfile and FEfile and len(Z_rot.T)==len(Z_ext.T):   #only analyse when both FE and Twist are from the same FOV
        for i,row in enumerate(Z_ext.T[:]):
            Pars = Tools.default_pars() 
            Pars['Filename']=FEfile
            Pars['date'] = FEfile.split("data",1)[0][:-1]  # date of measurement
            Pars['data'] = FEfile.split("data",1)[1][1:-4]  # measurement
            Pars['bead'] = i  # bead
            Pars['beads'] = len(Z_ext.T)  # total number of beads
            Pars['points'] = len(Z_ext)  # number of data-points
            Pars['points_exp'] = np.max(T) * Framerate  # number of expected data-points
            Pars['points_frac'] = Pars['points'] / Pars['points_exp'] # fraction of data-points
            Pars['X0_um'] = np.nanmean(X_ext[:,i])  # global X-position (um)
            Pars['Y0_um'] = np.nanmean(Y_ext[:,i])  # global Y-position (um)
            Pars['Z0_um'] = np.nanmean(Z_ext[:,i])  # global Z-position (um)
            Pars['dZ_um'] = abs(np.percentile(Z_ext[:,i],1)-np.percentile(Z_ext[:,i],99))  # absolute Z-extension after drift correction(um)

            if np.max(F)<MinFitForce:
                continue
            try: popt_S,pcov_S = Tools.fit_wlc(F,row,Pars, MinFitForce, MaxFitForce)
            except:
                print(FEfile, 'bead', i, '>>>>>>>> fit error! <<<<<<<')
                continue
            Pars['Radius_nm'] = Tools.fit_circle(X_rot[:,i],Y_rot[:,i]) # Radius of xy rotation fit with a circle
            
            fig1 = plt.figure()
            fig1.suptitle('Bead ' + str(i), y=.99)  
            ax1 = fig1.add_subplot(1, 2, 1)
            ax1.set_title('Twist ' + Rotationfile[:-4])
            ax1.set_ylabel(r'Extension')
            ax1.set_xlabel(r'Sigma')
            ax1.set_ylim(-200, Pars['L_bp']*Pars['dsDNA_nm_bp']*1.2)
            ax1.scatter(Sigma, Z_rot.T[i,:] - np.min(Z_rot.T[i,:]), color='blue',  s=5)
            ax1.text(np.min(Sigma), Pars['L_bp']*Pars['dsDNA_nm_bp']*1.1, 'F = '+str(round(RotationForce,2))+' pN', fontsize=12, verticalalignment='center', horizontalalignment='left')
            ax2 = fig1.add_subplot(1, 2, 2)
            ax2.set_title('Force Extension ' + FEfile[:-4])
            ax2.set_ylabel('Force')
            ax2.set_xlabel('Extension')
            ax2.scatter(Z_ext.T[i,:]-Pars['z0_nm'], F, color='red',  s=5)
            ax2.plot(Tools.wlc(F, Pars),F, color='black')  #plots the WLC
            ax2.text(np.min(Z_ext.T[i,:])-Pars['z0_nm'], np.max(F), 'P = '+str(round(Pars['P_nm'],1))+'±'+ str(round(np.sqrt(pcov_S[0,0]),1))+' nm', fontsize=12, verticalalignment='center', horizontalalignment='left')
            ax2.text(np.min(Z_ext.T[i,:])-Pars['z0_nm'], np.max(F)/10*9, 'S = ' + str(round(Pars['S_pN'],1)), fontsize=12, verticalalignment='center', horizontalalignment='left')
            fig1.savefig(Filename[:-4]+'_'+str(i)+'.png')
            plt.close()
            ###Remove bad fits:              
            if Pars['P_nm']>50*np.sqrt(pcov_S[0,0]) or Pars['P_nm'] < 15 or Pars['P_nm'] > 150:
                continue        
            ###Save data without duplicates            
            if DataOut!=[] and DataOut[-1]['Filename']!=Filename:
                DataOut.append(Pars)
    
    else:
        for i,row in enumerate(Z_ext.T[:]):
            Pars = Tools.default_pars() 
            Pars['Filename']=FEfile
            Pars['date'] = FEfile.split("data",1)[0][:-1]  # date of measurement
            Pars['data'] = FEfile.split("data",1)[1][1:-4]  # measurement
            Pars['bead'] = i  # bead
            Pars['beads'] = len(Z_ext.T)  # total number of beads
            Pars['points'] = len(Z_ext)  # number of data-points
            Pars['points_exp'] = np.max(T) * Framerate  # number of expected data-points
            Pars['points_frac'] = Pars['points'] / Pars['points_exp'] # fraction of data-points
            Pars['X0_um'] = np.nanmean(X_ext[:,i])  # global X-position (um)
            Pars['Y0_um'] = np.nanmean(Y_ext[:,i])  # global Y-position (um)
            Pars['Z0_um'] = np.nanmean(Z_ext[:,i])  # global Z-position (um)
            Pars['dZ_um'] = abs(np.percentile(Z_ext[:,i],1)-np.percentile(Z_ext[:,i],99))  # absolute Z-extension after drift correction(um)

            if np.max(F)<MinFitForce:
                continue
            try: popt_S,pcov_S = Tools.fit_wlc(F,row,Pars, MinFitForce, MaxFitForce)
            except:
                print(FEfile, 'bead', i, '>>>>>>>> fit error! <<<<<<<')
                continue
            
            ###Remove bad fits:  
            #if Pars['P_nm']>50*np.sqrt(pcov_S[0,0]) or Pars['P_nm'] < 15 or Pars['P_nm'] > 150:
                #continue        
            ###Save data without duplicates
            if DataOut[-1]['Filename']!=Filename:
                DataOut.append(Pars)

#            fig1 = plt.figure()
#            fig1.suptitle('Bead ' + str(i), y=.99)  
#            ax2 = fig1.add_subplot(1,1,1)
#            ax2.set_title('Force Extension ' + FEfile[:-4])
#            ax2.set_ylabel('Force')
#            ax2.set_xlabel('Extension')
#            ax2.scatter(Z_ext.T[i,:] - Pars['z0_nm'], F, color='red',  s=5)
#            ax2.plot(Tools.wlc(F, Pars),F, color='black')  #plots the WLC
#            ax2.text(np.min(Z_ext.T[i,:]-Pars['z0_nm']), np.max(F), 'P = '+str(round(Pars['P_nm'],1))+'±'+ str(round(np.sqrt(pcov[0,0]),1))+' nm', fontsize=12, verticalalignment='center', horizontalalignment='left')
#            ax2.text(np.min(Z_ext.T[i,:]-Pars['z0_nm']), np.max(F)/10*9.3, 'S = ' + str(round(Pars['S_pN'])), fontsize=12, verticalalignment='center', horizontalalignment='left')
#            fig1.savefig(FEfile[:-4]+'_'+str(i)+'.png')
#            plt.close()
            #print('Offset = ',Pars['z0_nm'], 'P = ', Pars['P_nm'], 'S = ', Pars['S_pN'])
       
Tools.save_df("DataOut.txt", DataOut) #Saves all parameters from Pars.
