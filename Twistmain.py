# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 20:58:16 2018

@author: Klaas
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 16:07:05 2018

@author: nhermans

"""

import os 
import Tools
import matplotlib.pyplot as plt

folder = r'G:\Klaas\Tweezers\Doxo Project\2018_06_15_197x25_Doxo_Constrained\CorrectedDat\FOV3'
filenames = os.listdir(folder)
os.chdir(folder)
Pars = Tools.default_pars()    
Rotationfile = False
FEfile = False
        
for Filenum, Filename in enumerate(filenames):
    if Filename[-4:] != '.dat':
        continue
    data, headers = Tools.read_dat(Filename)
    
    try: 
        X,Y,Z_rot, Rot, T = Tools.get_rotation_data(data, headers)
        Rotationfile = True
        Sigma = Rot/(Pars['L_bp']/Pars['Pitch_nm'])
    except: 
        X,Y,Z_ext, F, T = Tools.get_force_data(data, headers)
        FEfile = True
 
    if Rotationfile and FEfile:
        for i,row in enumerate(Z_rot.T[:]):
            fig1 = plt.figure()
            fig1.suptitle('bead '+ str(i), y=.99)  
            ax1 = fig1.add_subplot(1, 2, 1)
            ax1.set_title('Twist')
            ax1.set_ylabel(r'Extension')
            ax1.set_xlabel(r'Sigma')
            ax1.scatter(Sigma, row, color='blue', lw=0.1, s=5)
            ax2 = fig1.add_subplot(1, 2, 2)
            ax2.set_title('Force Extension')
            ax2.set_ylabel('Force')
            ax2.set_xlabel('Extension')
            ax2.scatter(Z_ext.T[i,:], F, color='blue', lw=0.1, s=5)
            fig1.show()
        Rotationfile, FEfile = False, False
