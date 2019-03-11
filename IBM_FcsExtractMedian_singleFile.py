11# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 21:24:20 2019

@author: Trevor Ho
"""

try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle

import os
from FlowCytometryTools.core.gates import CompositeGate
from FlowCytometryTools import FCMeasurement, ThresholdGate
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt


# Provide the Directories containing FCS files and the output Directories for saving data
dataRootDir = r'W:\Data storage & Projects\PhD Project_Yiyu\Transposase Project'
dataFolderDir = 'BM003'
fcsFolderDir = 'IBM_BM003_TestRun1_FCS'
filename = '20190214_IBM_BM003_TestRun1_Experiment_Group_IBMc101+IBMc120'

# Obtain filename handle
def fileNameParser(string=str):
        splitFN = string.split('Experiment_Group_')
        return splitFN[1]
    
fscid = fileNameParser(filename)

# Processing
fcsfilename = filename + '.fcs'
datafile = os.path.join(dataRootDir,dataFolderDir,fcsFolderDir,fcsfilename)
dataSample = FCMeasurement(ID=fscid, datafile=datafile)

# Gating
fsc_gate = ThresholdGate(1000.0, ['FSC-H'], region='above')
ssc_gate = ThresholdGate(1000.0, ['SSC-H'], region='above')
rfpA_gate = ThresholdGate(1.0, ['RFP2-A'], region='above')
fscssc_gate = CompositeGate(fsc_gate,'and',ssc_gate)
dataSample = dataSample.gate(fscssc_gate)
dataSample = dataSample.gate(rfpA_gate)


#fig1 = plt.figure('figure 1')
#rfpdata = dataSample.data['RFP2-H']
#rfplist = rfpdata.values
#rfpLogMin, rfpLogMax = np.log10(rfplist.min()),np.log10(rfplist.max())
#rfpnewBins = np.logspace(rfpLogMin, rfpLogMax,100)
#sns.distplot(rfplist, bins=rfpnewBins, kde=False, color='r')
#plt.xscale('log')
#plt.xlim(1, 1e6)

data = dataSample.data[['FSC-H','RFP2-H']]
df = pd.DataFrame(data, columns=["x", "y"])
#sns.kdeplot(df.x, df.y)
#sns.jointplot(x="RFP2-H", y="FSC-H", data=data, kind="kde");