# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 20:29:14 2019

@author: Trevor Ho
"""

import os
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle

# Provide the Directories containing FCS files and the output Directories for saving data
dataRootDir = r'W:\Data storage & Projects\PhD Project_Yiyu\Transposase Project'
dataFolderDir = 'FC011'
fileName = 'IBM_FC011R2PI24_Data'
pickleFilename = fileName + '.pickle'
pickleDir = os.path.join(dataRootDir,dataFolderDir,pickleFilename)


# Load Gated data
pickle_handle = open(pickleDir,"rb")
plate = pickle.load(pickle_handle)

# Create a list of wellNames
wellNames = []
for i in range(8):
    char = chr(i+65)
    for j in range(12):        
        wellNames.append(char+str(j+1))    

#%%

fig, ax = plt.subplots(8, 12, sharex=True)
fig.set_size_inches(12, 8)

for i in range(8):
    for j in range(12):
        if plate.get(wellNames[(12*i+j)]) is not None:
            rfpdata = plate[wellNames[(12*i+j)]].data['RFP2-H']
            rfplist = rfpdata.values
            rfpLogMin, rfpLogMax = np.log10(rfplist.min()),np.log10(rfplist.max())
            rfpnewBins = np.logspace(rfpLogMin, rfpLogMax,100)
            sns.distplot(rfplist, bins=rfpnewBins, kde=False, color='r', ax= ax[i,j])
            plt.xscale('log')
            plt.xlim(1, 1e6)
        
        if i == 7 and j == 0: 
            ax[i,j].axes.get_yaxis().set_visible(False)
        else:
            ax[i,j].axes.get_xaxis().set_visible(False)
            ax[i,j].axes.get_yaxis().set_visible(False)

#%%
figFilename = fileName + '_rfpHisto.png'
figDir = os.path.join(dataRootDir,dataFolderDir,figFilename)
fig.savefig(figDir)