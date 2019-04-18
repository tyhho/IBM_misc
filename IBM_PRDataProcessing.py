# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 17:56:09 2019

@author: Trevor Ho

This script process all plate reader data excel files and output a single csv file
that contains information on corrected OD600, corrected red fluorescence
and corrected fluorescence/OD600, it also adds metadata to the file regarding
the Run number, the SampleID and the condition

This script is tailored for FC014 Run 1

It should be noted that this way of processing is not sustainable because the
conditions of the plates will only get more and more complex. Needs to rework
the code to generate something that automatically match plate number and 
plate metadata file.
"""

import os
import pandas as pd

# Function for renaming Index
def renameIndex(indexStr=str):
    indexStr = indexStr[0] + str(int(indexStr[1:3]))
    return indexStr

# TODO: Specify folder location
# TODO: Prepare file of metadata
dataRootDir=r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir='FC014'
metafilename = 'IBM_FC014R1_PRPlateMetadata.xlsx'
outputCSV = 'IBM_FC014R1_PRData.csv'


#%% Since there is only a single metafile in this experiment

# Read PR Metadata Excelfile
metadataDir = os.path.join(dataRootDir,dataFolderDir,metafilename)
metaxls= pd.ExcelFile(metadataDir)
metadata = {sheet:metaxls.parse(sheet) for sheet in metaxls.sheet_names}

# Process the df of metadata & replace the old dictionary with new metadata dataframes
for meta_property, metadf_96format in metadata.items():
    metadf = pd.DataFrame(columns = [])
    for char in ['A','B','C','D','E','F','G','H']:
        metadf_96_row = metadf_96format.transpose()[char]
        metadf_96_row = metadf_96_row.to_frame(meta_property)
        indexList = metadf_96_row.index.tolist()
        
        # Retrieve information of well and set it as index
        for columnIndex in range(len(indexList)):
            indexList[columnIndex] = char + str(indexList[columnIndex])
        metadf_96_row['Well']=indexList
        metadf_96_row.set_index('Well', inplace=True)
        
        # Append to major dataframe of the single metadata
        metadf = metadf.append(metadf_96_row,sort=False)
    metadata[meta_property] = metadf

#%%

alldata = pd.DataFrame(columns = [])

folderDir = os.path.join(dataRootDir, dataFolderDir)
# This script loops through all PR data automatically instead of mapping each file with a "blank finder"
for file in os.listdir(folderDir):
     filename = os.fsdecode(file)
     if filename.endswith(".xlsx") and filename[0:2] == 'PR':
         
         filenameParts_info = filename.split('_')[3]
         run_no = filenameParts_info[1]
         ind_time = filenameParts_info.split('PI')[1].split('P')[0]
         ind_cond_code = filenameParts_info.split('PI')[1].split('P')[1][0]
         excelDir = os.path.join(folderDir,filename)
         data=pd.read_excel(excelDir,sheet_name = 'End point',index_col=0,skiprows=12)
         data = data.drop(['Content','Raw Data (600 1)','Raw Data (584 2)'], axis=1)
         
         # Update the index so it matches the conventional nomenclature
         indexList = data.index.tolist()
         for n, ind in enumerate(indexList):
             indexList[n]= renameIndex(ind)
         data.set_index(pd.Index(indexList),inplace=True)
         data.rename(columns={'Blank corrected based on Raw Data (600 1)': 'PR_Corrected OD600',
                              'Blank corrected based on Raw Data (584 2)': 'PR_Corrected Fluo (a.u.)',
                              'RF/OD600': 'PR_Corrected Red Fluo/OD600 (a.u.)'}, inplace=True)
         data['Run'] = run_no
         data['Post-induction (hrs)'] = ind_time
         if ind_cond_code == '1':
             data['Induction'] = '1 mM arabinose'
         elif ind_cond_code == '2':
             data['Induction'] = '1 mM arabinose + 0.1 mM caffeine'
         else:
             data['Induction'] = 'error'
         
            # Append all metadata to dataframe with plate reader data
         for meta_property, metadf_96format in metadata.items():
             data = data.merge(metadf,left_index=True,right_index=True)
        
        # TODO: Check the row of blank that should be removed
         # Remove data of the blank well 'H12'
         data = data.drop('H12',axis=0)
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