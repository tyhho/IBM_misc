# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 17:31:09 2019

@author: Trevor Ho
"""
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle

import os
from FlowCytometryTools.core.gates import CompositeGate
from FlowCytometryTools import FCPlate, ThresholdGate
import pandas as pd

# Parser Function for FlowCytometryTools
def fcsNameParser(string=str):
        splitFN = string.split('Experiment_Group_')
        splitFN2 = splitFN[1].split('.fcs')
        return splitFN2[0]

# Specify folder location
    # Each folder must contain only fcs files that end with well location
dataRootDir = r'W:\Data storage & Projects\PhD Project_Yiyu\Transposase Project'
dataFolderDir = 'FC011'

# Create dict with information of FCS folder name and growth method (for deciding what gates to use)
plateList = ['IBM_FC011R2PI5',
             'IBM_FC011R2PI24'
             ]

#%% Core Processing Codes

for plateNameCore in plateList:
    #plateNameCore = 'TCI_FC004_Run1_GM1_pTCI028' #for debugging
    fcsFolderDir = plateNameCore + '_FCS'
    outputFilename = plateNameCore + '_Data'
    
    datadir = os.path.join(dataRootDir,dataFolderDir,fcsFolderDir)
    plate = FCPlate.from_dir(ID='Plate', path=datadir, parser=fcsNameParser, position_mapper='name')
    
    # Gating
    fsc_gate = ThresholdGate(1000.0, ['FSC-H'], region='above')
    ssc_gate = ThresholdGate(1000.0, ['SSC-H'], region='above')
    
    rfpA_gate = ThresholdGate(1.0, ['RFP2-A'], region='above')
    fscssc_gate = CompositeGate(fsc_gate,'and',ssc_gate)
    plate = plate.gate(fscssc_gate)
    plate = plate.gate(rfpA_gate)
    
    # Calculate Median from data
    def calculate_median_Y2(well):
        return well.data['RFP2-H'].median()
    
    dfAll = {} 
    dfAll[0] = plate.apply(calculate_median_Y2)
    dfAll[1] = plate.counts()
    
    # Save Median values to an excel file
    sheetName = ['RFP','Count']
    outputFilenameXlsx = outputFilename + '.xlsx'
    outputDir = os.path.join(dataRootDir,dataFolderDir,outputFilenameXlsx)
    writer = pd.ExcelWriter(outputDir, engine='xlsxwriter')
    
    for i in range(2):
        dfAll[i].to_excel(writer, sheet_name=sheetName[i])
    writer.save()
    
    # Export the plate data as a Pickle Object
    outputFilenamePickle = outputFilename + '.pickle'
    pickleDir = os.path.join(dataRootDir,dataFolderDir,outputFilenamePickle)
    
    with open(pickleDir, 'wb') as handle:
        pickle.dump(plate, handle, protocol=pickle.HIGHEST_PROTOCOL)
