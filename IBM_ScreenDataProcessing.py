# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 14:54:38 2019

@author: Trevor Ho
"""

import os
import pandas as pd
import numpy as np
import IBM_CustomFunctions as cf

# Function for renaming Index
def renameIndex(indexStr=str):
    indexStr = indexStr[0] + str(int(indexStr[1:3]))
    return indexStr

#TODO: Specify directories for input and output data
    # all csv files that contain df data of fluorescence with inducer conc must be in this folder
    # this is also the folder where the output will be deposited
dataRootDir = r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir = 'BM005'
dataFolder2Dir = 'Screen1'
outputCSV = 'BM005_Screen1_pooledData.csv'

#%%
# Create dict for handling multiple files at once
filenameCoreList = {'PR_BM005_Screen1_PI18_': ['P1','P2']
                    }

# Create list for induction OFF & ON states
inductionList = ['-','+']


for fNVar, plateList in filenameCoreList.items():
    mergedData = pd.DataFrame(columns=[])
    for plateNo in plateList:
        plateData = pd.DataFrame(columns=[])
        
        for indCon in inductionList:
            
            excelFN = fNVar + plateNo + indCon + '.xlsx'
            dataDir = os.path.join(dataRootDir, dataFolderDir,dataFolder2Dir,excelFN)
            data = pd.read_excel(dataDir,sheet_name = 'End point',skiprows=12)
            data = data.drop(['Content'], axis=1)
            
            #Background Correction
            #TODO: Provide excel file to look for control data
            ctrlExcelFN = fNVar+ 'Blank.xlsx'
            ctrlDir = os.path.join(dataRootDir, dataFolderDir,dataFolder2Dir,ctrlExcelFN)
            ctrlData = pd.read_excel(ctrlDir,sheet_name = 'End point',index_col=0,skiprows=12)

            #TODO: Provide the well location of the blank
            blankOD = ctrlData.iloc[0][1]
            blankFluo = ctrlData.iloc[0][2]
            
            data['Raw Data (600 1)'] -= blankOD
            data['Raw Data (584 2)'] -= blankFluo
            
            plateData['Well'] = data['Well']
            
            if indCon == '-':
                plateData['-Ind OD'] = data['Raw Data (600 1)']
                plateData['-Ind Fluo'] = data['Raw Data (584 2)']
            else:
                plateData['+Ind OD'] = data['Raw Data (600 1)']
                plateData['+Ind Fluo'] = data['Raw Data (584 2)']
            
            # Update the index so it matches the conventional nomenclature
            indexList = plateData['Well'].tolist()
            for n, ind in enumerate(indexList):
                indexList[n]= renameIndex(ind)
            plateData['Well'] = pd.Series(indexList)
            plateData['Plate'] = plateNo
            
            # Read Metadata Excelfile
            metadataNameCore = 'BM005_Screen1_' + plateNo + '_PRPlateMetadata.xlsx'
            metafilename = metadataNameCore
            metadataDir = os.path.join(dataRootDir,dataFolderDir,dataFolder2Dir,metafilename)
            metadf = cf.metadata_to_metadf(metadataDir)

        plateData = plateData.merge(metadf,left_on='Well',right_index=True)
        mergedData = mergedData.append(plateData,ignore_index=False,sort=False)
        mergedData['-Ind Fluo/OD'] =  mergedData['-Ind Fluo']/mergedData['-Ind OD']
        mergedData['+Ind Fluo/OD'] =  mergedData['+Ind Fluo']/mergedData['+Ind OD']
        mergedData['IndFoldChange'] = mergedData.apply(lambda mergedData :mergedData['-Ind Fluo/OD']/mergedData['+Ind Fluo/OD'] if mergedData['-Ind Fluo/OD']!=0 else np.nan,axis=1) 

    # Export Data File, one for each induction time point
    outputCSV = fNVar + 'pooledData.csv'
    csvDir = os.path.join(dataRootDir,dataFolderDir,dataFolder2Dir,outputCSV)
    mergedData.to_csv(csvDir,header=True)

#mergedData = mergedData.sort_values(by=['+C Fluo/OD'],ascending=False)