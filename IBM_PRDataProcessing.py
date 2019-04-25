# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 17:56:09 2019

@author: Trevor Ho

This script process all plate reader data excel files and output a single csv file
that contains information on corrected OD600, corrected red fluorescence
and corrected fluorescence/OD600, it also adds metadata to the file regarding
the Run number, the SampleID and the condition

This script is tailored for FC014 Run 2

"""

import os
import pandas as pd
import IBM_CustomFunctions as cf

# TODO: Specify folder location
dataRootDir=r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir='FC014'

# TODO: Prepare file of metadata
# Check that the blank well has been labeled as "Blank"
metafilename = 'PRMD_IBM_FC014R2.xlsx'
outputCSV = 'IBM_FC014R2-4_PRData.csv'

# Read PR Metadata Excelfile
metadataDir = os.path.join(dataRootDir,dataFolderDir,metafilename)
metadata_df = cf.metadata_to_metadf(metadataDir)

#%%

alldata = pd.DataFrame(columns = [])

folderDir = os.path.join(dataRootDir, dataFolderDir)

'''
Strategy: 
    Loop through every single file starting with 'PR', check if filename contains the key
    If yes, process it and map the metadata onto it

'''

keylist = [
        'FC014R2',
        'FC014R3',
        'FC014R4'
           ]

for file in os.listdir(folderDir):
     filename = os.fsdecode(file)

     if filename.endswith(".xlsx") and filename.split('_')[0] == 'PR' \
     and (not keylist or any(key in filename for key in keylist)):
         
         # Read Excel file
         excelDir = os.path.join(folderDir,filename)
         data = pd.read_excel(excelDir,sheet_name = 'End point',index_col=0,skiprows=12)
         data = data.drop(['Content','Raw Data (600 1)','Raw Data (584 2)'], axis=1)
         
         # Update the index so it matches the conventional nomenclature
         data = cf.renameDfWellIndex(data)
         
         # Extract information from filename
         filename = filename.split('.xlsx')[0][3:]    #Removes 'PR' from filename for downstream analysis
         
         run_no = int(filename.split('R')[1][0])
                  
         # Check if induction time and plate no info are in the filename
         if filename.find('PI')>=0:
             ind_time_frag = filename.split('PI')[1]
             if ind_time_frag.find('P') >=0:
                 ind_time, plate_no = int(ind_time_frag.split('P')[0]), int(ind_time_frag.split('P')[1])
             else:
                 ind_time = int(ind_time_frag)
         elif filename.find('P') >=0:
             plate_no = int(filename.split('P')[1])
         
         
         data.rename(columns={'Blank corrected based on Raw Data (600 1)': 'PR_Corrected OD600',
                              'Blank corrected based on Raw Data (584 2)': 'PR_Corrected Fluo (a.u.)',
                              'RF/OD600': 'PR_Corrected Red Fluo/OD600 (a.u.)'}, inplace=True)
         data['Run'] = run_no
         data['Post-induction (hrs)'] = ind_time
         data = data.merge(metadata_df,left_index=True,right_index=True)
         
         # TODO: Check that blank info is in metadata excel file
         blank_list = data.index[data['SampleID']=='Blank'].tolist()
         data = data.drop(blank_list,axis=0)
         
         data['PR_Well'] = data.index
         alldata = alldata.append(data,ignore_index=True,sort=None)

#%%
# Create dict for handling multiple files at once
    # For each file, provide the key to find blank OD and background fluorescence
#fileList = {'': ['PR_TCI_FC006R1P1','A01'],
#            '': ['PR_TCI_FC006R1P1','A01']
#                    }

#for filenameCore, ctrlInfo in fileList.items():
    # Loop starts here
    # Process data (median fluorescence) from 96 well plate format into Seaborn-friendly format

#    # Read control data from xlsx file
#    ctrlExcelFN = ctrlInfo[0] + '.xlsx' # extract control filename
#    ctrlExcelDir = os.path.join(dataRootDir, dataFolderDir,ctrlExcelFN)
#    ctrlData = pd.read_excel(ctrlExcelDir,sheet_name = 'End point',index_col=0,skiprows=12)
#    ctrlData = ctrlData.drop(['Content'], axis=1)

#    blankIndex = ctrlData.index.get_loc(ctrlInfo[1])
#    blankOD = ctrlData.iloc[blankIndex][0]  #extract blank OD
#    blankGF = ctrlData.iloc[blankIndex][1]  #extract blank green fluorescence
#    blankRF = ctrlData.iloc[blankIndex][2]  #extract blank red fluorescence
    
    # Read data to be processed from xlsx file
#    filename = filenameCore + '.xlsx'
#    dataDir = os.path.join(dataRootDir, dataFolderDir,filename)
#    data = pd.read_excel(dataDir,sheet_name = 'End point',index_col=0,skiprows=12)
#    data = data.drop(['Content'], axis=1)
    
    # Background correction
#    data['Raw Data (600 1)'] -= blankOD
#    data['Raw Data (485 2)'] -= blankGF
#    data['Raw Data (584 3)'] -= blankRF
    
    # Update the index so it matches the conventional nomenclature
#    indexList = data.index.tolist()
#    for n, ind in enumerate(indexList):
#        indexList[n]= renameIndex(ind)
#    data.set_index(pd.Index(indexList),inplace=True)
#    
#    data.rename(columns={'Raw Data (600 1)': 'PR_Corrected OD600',
#                                    'Raw Data (485 2)': 'PR_Corrected Green Fluorescence (a.u.)',
#                                    'Raw Data (584 3)': 'PR_Corrected Red Fluorescence (a.u.)'}, inplace=True)
#
    
# Save file as CSV
csvDir = os.path.join(folderDir,outputCSV)
alldata.to_csv(csvDir,header=True)