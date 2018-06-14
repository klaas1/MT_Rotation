# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 16:07:05 2018

@author: nhermans
"""

import os 
import Tools
import numpy as np
import matplotlib.pyplot as plt

folder = r'G:\Klaas\Tweezers\Doxo Project\2018_06_14_Doxo_4kb_constrained\DoxFOV\Dox'
filenames = os.listdir(folder)
os.chdir(folder)
Pars = Tools.default_pars()    
        
for Filenum, Filename in enumerate(filenames):
    if Filename[-4:] != '.dat':
        continue
    data, headers = Tools.read_dat(Filename)
    try: 
        X,Y,Z, Rot, T = Tools.get_rotation_data(data, headers)
        Rotationfile = True
    except: 
        X,Y,Z, F, T = Tools.get_force_data(data, headers)
        Rotationfile = False
    Sigma = Rot/(Pars['L_bp']/Pars['Pitch_nm'])
    
    fig1 = plt.figure()
    #fig1.suptitle(i, y=.99)    
    for i,row in enumerate(Z.T[:]):
        fig1 = plt.figure()
        if Rotationfile: Tools.plot_sigma(Sigma,row,fig1)
        else: Tools.plot_force(F,row,fig1)
        
