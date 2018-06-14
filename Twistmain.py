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

folder = r'G:\Klaas\Tweezers\Doxo Project\2018_06_14_Doxo_4kb_constrained\DoxFOV\Dox'
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
            fig1.suptitle(i, y=.99)  
            Tools.plot_sigma(Sigma,row,fig1)
            Tools.plot_force(F,Z_ext.T[i,:],fig1)
