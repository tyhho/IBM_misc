# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 19:49:17 2019

@author: Trevor Ho

Provides Custom Functions to commonly uesd scripts in IBM

"""

import pandas as pd


'''Read a single metadata file and convert it into a dictionary of metadata for 96 well plate'''

def metadata_to_metaDict(metadataDir):

    metaxls= pd.ExcelFile(metadataDir)
    dict_of_metadata = {sheet:metaxls.parse(sheet) for sheet in metaxls.sheet_names}
    
    # Process the df of metadata & create a metadata dictionary
        # Key = metadata field
        # Value = metadata dataframes
    for meta_property, metadf_96format in dict_of_metadata.items():
        metadf = pd.DataFrame(columns = [])
        for char in ['A','B','C','D','E','F','G','H']:
            metadf_96_row = metadf_96format.transpose()[char]
            metadf_96_row = metadf_96_row.to_frame(meta_property)
            
            # Retrieve information of well and set it as index
            indexList = metadf_96_row.index.tolist()
            for columnIndex in range(len(indexList)):
                indexList[columnIndex] = char + str(indexList[columnIndex])
            metadf_96_row['Well']=indexList
            metadf_96_row.set_index('Well', inplace=True)
            
            # Append to major dataframe of the single metadata
            metadf = metadf.append(metadf_96_row,sort=False)
        dict_of_metadata[meta_property] = metadf

    return dict_of_metadata

'''Read a single metadata file and convert it into a dataframe for merging'''

def metadata_to_metadf(metadataDir):
    dict_of_metadata = metadata_to_metaDict(metadataDir)
    meta_df_list = list(dict_of_metadata.values())
    all_meta_df = meta_df_list[0]
    meta_df_list = meta_df_list[1:]
    for meta_df in meta_df_list:
        all_meta_df = all_meta_df.merge(meta_df,left_index=True,right_index=True,sort=False)
    return all_meta_df


'''Rename index in a dataframe from plate reader exported Excel file to normal well index (e.g 'A01' to 'A1')'''
# Note: only take df as input

def renameDfWellIndex(data):
     indexList = data.index.tolist()
     for n, ind in enumerate(indexList):
         indexList[n]= ind[0] + str(int(ind[1:3]))
     data.set_index(pd.Index(indexList),inplace=True)
     return data

#%% Rename files
#
#import os
#
## TODO: Specify folder location
#dataRootDir=r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
#dataFolderDir=
#
#folderDir = os.path.join(dataRootDir, dataFolderDir)
## This script loops through all PR data automatically instead of mapping each file with a "blank finder"
#for file in os.listdir(folderDir):
#     filename = os.fsdecode(file)
#     if filename.endswith(".RUC") and filename[0:2] == 'PR':
#         filenameParts = filename.split('_')
#         newFilename = filenameParts[0] + '_' + filenameParts[2] + '_' + filenameParts[3] + filenameParts[4]
#         oldFileDir = os.path.join(dataRootDir,dataFolderDir,filename)
#         newFileDir = os.path.join(dataRootDir,dataFolderDir,newFilename)

#TODO: only uncomment the following line until the newFilename was expected
##         os.rename(oldFileDir,newFileDir)