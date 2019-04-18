# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 14:54:38 2019

@author: Trevor Ho
"""

import os
import pandas as pd
import numpy as np

# Function for renaming Index
def renameIndex(indexStr=str):
    indexStr = indexStr[0] + str(int(indexStr[1:3]))
    return indexStr

#TODO: Specify directories for input and output data
    # all csv files that contain df data of fluorescence with inducer conc must be in this folder
    # this is also the folder where the output will be deposited
dataRootDir = r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir = 'BM004'
dataFolder2Dir = 'Screen4'
outputCSV = 'BM004_Screen4_pooledData.csv'

#TODO: Provide excel file to look for control data
ctrlExcelFN = 'PR_BM004_Screen4_PI24_Blank.xlsx'
ctrlDir = os.path.join(dataRootDir, dataFolderDir,dataFolder2Dir,ctrlExcelFN)
ctrlData = pd.read_excel(ctrlDir,sheet_name = 'End point',index_col=0,skiprows=12)
#TODO: Provide the well location of the blank
blankOD = ctrlData.iloc[0][1]
blankFluo = ctrlData.iloc[0][2]


#%%
# Create dict for handling multiple files at once
filenameCoreList = {'PR_BM004_Screen4_PI24_': ['1']
                    }

# Create list for induction OFF & ON states
inductionList = ['-','+']


for fNVar, plateList in filenameCoreList.items():
    mergedData = pd.DataFrame(columns=['ID','Plate','Well','+C OD','+C Fluo'])
    for plateNo in plateList:
        plateData = pd.DataFrame(columns=['ID','Plate','Well','+C OD','+C Fluo'])
        
        for indCon in inductionList:
            excelFN = fNVar + plateNo + indCon + '.xlsx'
            dataDir = os.path.join(dataRootDir, dataFolderDir,dataFolder2Dir,excelFN)
            data = pd.read_excel(dataDir,sheet_name = 'End point',skiprows=12)
            data = data.drop(['Content'], axis=1)
            data['Raw Data (600 1)'] -= blankOD
            data['Raw Data (584 2)'] -= blankFluo
            
            plateData['Well'] = data['Well']
            
            if indCon == '-':
                plateData['-C OD'] = data['Raw Data (600 1)']
                plateData['-C Fluo'] = data[ 'Raw Data (584 2)']
            else:
                plateData['+C OD'] = data['Raw Data (600 1)']
                plateData['+C Fluo'] = data[ 'Raw Data (584 2)']
            
            # Update the index so it matches the conventional nomenclature
            indexList = plateData['Well'].tolist()
            for n, ind in enumerate(indexList):
                indexList[n]= renameIndex(ind)
            plateData['Well'] = pd.Series(indexList)
            plateData['Plate'] = plateNo
            plateData['ID'] = ( int(plateNo) -1 ) * 96 + np.linspace(1,96,96)
   
        mergedData = mergedData.append(plateData,ignore_index=False,sort=False)
        mergedData['-C Fluo/OD'] =  mergedData['-C Fluo']/mergedData['-C OD']
        mergedData['+C Fluo/OD'] =  mergedData['+C Fluo']/mergedData['+C OD']
        mergedData['IndFoldChange'] = mergedData.apply(lambda mergedData :mergedData['+C Fluo/OD']/mergedData['-C Fluo/OD'] if mergedData['-C Fluo/OD']!=0 else np.nan,axis=1) 


mergedData = mergedData.sort_values(by=['+C Fluo/OD'],ascending=False)
        
# Export Data File, one for each induction time point
csvDir = os.path.join(dataRootDir,dataFolderDir,dataFolder2Dir,outputCSV)
mergedData.to_csv(csvDir,header=True)