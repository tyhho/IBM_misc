# -*- coding: utf-8 -*-

"""
Created on Thu Mar 28 17:56:09 2019

@author: Trevor Ho

This script process all plate reader data excel files and output a single csv file
that contains information on corrected OD600, corrected red fluorescence
and corrected fluorescence/OD600, it also adds metadata to the file regarding
the Run number, the SampleID and the condition

"""

import os
import pandas as pd
import CustomFunctions as cf
import fnmatch

# TODO: Specify folder location
root_path = '..'
data_folder = 'FC007'
output_fn = 'IBM_FC007R4,5,7_PRData.csv'
# outputCSV = 'IBM_FC033R1,4,5_PRData.csv'


# TODO: Prepare file of metadata
# Check that the blank well has been labeled as "Blank"

# TODO: Use master_blank_info if all files use the same blank info
master_blank_info = []    #provide fn and well location of where the blank should be
alldata = pd.DataFrame(columns = [])
data_path = os.path.join(root_path, data_folder)
allFiles = os.listdir(data_path)

#%% Strategy 1

# Create dict for handling multiple files at once
    # The key specifies the filename search criteria
    # For each file, provide 2 info in a list
    # 1. metadata file
    # 2. info for blank, can be:
        # a. null, not giving a item 2 in list, the script will look for blank using the metadata provided and raise error if blank is not found
        # b. master_blank_info, when one blank info is used consistently across all files in one processing run 
        # c. a customized list of [filename, well] that specifies where the blank should be looked up
            # if 'filename' here is given as 'customFunction', then it will call the findBlankXlsx function to look for the blank
        # Having a blank info will force the script to use the blank given instead of the blank in metadata file
    
fileList = {
            'PR_IBM_FC007R[4,5,7]*P1.xlsx': ['PRMD_IBM_FC007R4P1.xlsx'],
            # 'PR_IBM_FC033R[1,4,5]*P1.xlsx': ['PRMD_IBM_FC033P1.xlsx'],
            # 'PR_IBM_FC033R[1,4,5]*P2.xlsx': ['PRMD_IBM_FC033P2.xlsx'],
            # 'PR_IBM_FC033R[1,4,5]*P3.xlsx': ['PRMD_IBM_FC033P3.xlsx'],
            # 'PR_IBM_FC033R[1,4,5]*P4.xlsx': ['PRMD_IBM_FC033P4.xlsx']
            }


# Process data (median fluorescence) from 96 well plate format into Seaborn-friendly format

for fnSearchSeq, metaInfo in fileList.items():

    # Get metadata
    metafn = metaInfo[0]
    metadataDir = os.path.join(root_path, data_folder, metafn)
    metadata_df = cf.metadata_to_metadf(metadataDir)
    
    # Try to get 'Blank' from metadatafile
    from_meta_blank_list = metadata_df.index[metadata_df['SampleID']=='Blank'].tolist()
    useBlankInSelf = True
    
    # Decides whether to use blank info or not
    try:
        blankInfo = metaInfo[1]
        # check if blankInfo was provided but given in wrong format
        if type(blankInfo) != list or len(blankInfo) != 2:
            raise ValueError('Blank Info should be a list with 1st item as filename and 2nd item as well location')
        # If not, assume everything is ok and continue to extract blank info
        blank_fn = blankInfo[0] # extract filename containing blank
        blank_file_path = os.path.join(root_path, data_folder, blank_fn)
        dataInBlank = pd.read_excel(blank_file_path, sheet_name = 'End point', index_col=0, skiprows=12)
        dataInBlank = dataInBlank.drop(['Content'], axis=1)

        blankIndex = dataInBlank.index.get_loc(blankInfo[1])
        blankOD = dataInBlank.iloc[blankIndex][0]  #extract blank OD
        blankRF = dataInBlank.iloc[blankIndex][1]  #extract blank red fluorescence
        
        del blank_fn, blank_file_path, dataInBlank, blankIndex
        useBlankInSelf = False
        
    except IndexError:
        # If there is no blank info given and no 'Blank' found in metadata
        if not from_meta_blank_list:
            raise ValueError('No \'Blank\' present in metadata. Blank Info needed')
            
    # Get all files that match the search criteria
    matched_files = fnmatch.filter(allFiles,fnSearchSeq)
    
    # Process each file
    for matched_fn in matched_files:
        
        # Read data
        data_path = os.path.join(root_path, data_folder, matched_fn)
        data = pd.read_excel(data_path, sheet_name = 'End point', index_col=0, skiprows=12)
        data = data.drop(['Content'], axis=1)
        
        # Update the index so it matches the conventional nomenclature
        data = cf.renameDfWellIndex(data)
        
        # Handle blanks 
        if useBlankInSelf == True:
            # Drop raw data that did not receive blank correction
            data = data.drop(['Raw Data (600 1)','Raw Data (584 2)'], axis=1)
            # Drop blank data
            data = data.drop(from_meta_blank_list,axis=0)
            # Rename cols
            data.rename(columns={'Blank corrected based on Raw Data (600 1)': 'PR_Corrected OD600',
                                 'Blank corrected based on Raw Data (584 2)': 'PR_Corrected Red Fluorescence (a.u.)',
                                 'RF/OD600': 'PR_Corrected Red Fluo/OD600 (a.u.)'}, inplace=True)

        else:
            # blank correction
            data['Raw Data (600 1)'] -= blankOD
            data['Raw Data (584 2)'] -= blankRF
            # Rename cols
            data.rename(columns={'Raw Data (600 1)': 'PR_Corrected OD600',
                                 'Raw Data (584 2)': 'PR_Corrected Red Fluorescence (a.u.)'}, inplace=True)
            # Drop any possible remaining blank data
            if from_meta_blank_list:
                data = data.drop(from_meta_blank_list,axis=0)
        
            # Calculate Fluo/OD
            data['PR_Corrected Red Fluo/OD600 (a.u.)'] = data['PR_Corrected Red Fluorescence (a.u.)'] / data['PR_Corrected OD600']
        
        # Add extra metadata based on info in filename
        # Extract information from filename
        fnInfo = matched_fn.split('.xlsx')[0][3:]    #Removes 'PR' from filename for downstream analysis
        run_no = int(fnInfo.split('R')[1][0])
                  
        # Check if induction time and plate no info are in the filename
        if fnInfo.find('PI')>=0:
            ind_time_frag = fnInfo.split('PI')[1]
            if ind_time_frag.find('P') >=0:
                try:
                    plate_no = int(ind_time_frag.split('P')[1])
                except ValueError:
                    plate_no = ind_time_frag.split('P')[1]
            try:
                ind_time = int(ind_time_frag.split('P')[0])
            except ValueError:
                ind_time = ind_time_frag.split('P')[0]
        elif fnInfo.find('P') >=0:
            # pass
            plate_no = int(fnInfo.split('P')[1])
        
        data['Run'] = run_no
        data['Post-induction (hrs)'] = ind_time
        
        try:
            data['PR_Plate'] = plate_no
        except NameError:
            pass

        # Append metadata
        data = data.merge(metadata_df,left_index=True,right_index=True)
        
        data['PR_Well'] = data.index
        alldata = alldata.append(data,ignore_index=True,sort=None)

#%% Strategy 2

# For everyfile pattern that matches the filename, look for corresponding metadatafile and blank excel file

# metadatafnCore = 'PRMD_IBM_FC029R4'
metadatafnCore = 'PRMD_IBM_FC007'
blank_well = 'F08'
blank_plate = '1'

fnSearchSeqList = [
            'PR_IBM_FC007R[1-3]*.xlsx',
            # 'PR_IBM_FC034R[1-3]*P[2-5].xlsx',
#            'PR_IBM_FC021R[2,3,4]*P2.xlsx',
#            'PR_IBM_FC021R[3,4,5]*P3.xlsx',
#            'PR_IBM_FC021R[4,6,7]*P4.xlsx'
#            'PR_IBM_FC021R[2-7]*.xlsx'
            ]

# fnSearchSeqList = [
#             'PR_IBM_FC029R[4,5,7]*.xlsx',
#            'PR_IBM_FC021R[2,3,4]*P2.xlsx',
#            'PR_IBM_FC021R[3,4,5]*P3.xlsx',
#            'PR_IBM_FC021R[4,6,7]*P4.xlsx'
#            'PR_IBM_FC021R[2-7]*.xlsx'
            # ]


#fnSearchSeq = 'PR_IBM_FC021R[2-7]*.xlsx'

# Process data (median fluorescence) from 96 well plate format into Seaborn-friendly format

for fnSearchSeq in fnSearchSeqList:
    
    # Get all files that match the search criteria
    matchedFiles = fnmatch.filter(allFiles,fnSearchSeq)
    
    # Process each file
    for matchedfn in matchedFiles:
        
        # Get metadata
        metafn = cf.findMetaXlsx(matchedfn, metadatafnCore)
        metadataDir = os.path.join(root_path, data_folder, metafn)
        metadata_df = cf.metadata_to_metadf(metadataDir)
        
        # Get blank info
        # If not, assume everything is ok and continue to extract blank info
        blankFN = cf.findBlankXlsx(matchedfn, blank_plate) #get filename using custom function
        blankDir = os.path.join(root_path, data_folder, blankFN)
        dataInBlank = pd.read_excel(blankDir,sheet_name = 'End point',index_col=0,skiprows=12)
        dataInBlank = dataInBlank.drop(['Content'], axis=1)
    
        blankIndex = dataInBlank.index.get_loc(blank_well)
        blankOD = dataInBlank.iloc[blankIndex][0]  #extract blank OD
        blankRF = dataInBlank.iloc[blankIndex][1]  #extract blank red fluorescence
        
        del blankFN, blankDir, dataInBlank,blankIndex
        
        # Read datas
        dataDir = os.path.join(root_path, data_folder, matchedfn)
        data = pd.read_excel(dataDir, sheet_name = 'End point', index_col=0, skiprows=12)
        data = data.drop(['Content'], axis=1)
        
        # Update the index so it matches the conventional nomenclature
        data = cf.renameDfWellIndex(data)
        
        # Blank correction
        data['Raw Data (600 1)'] -= blankOD
        data['Raw Data (584 2)'] -= blankRF
        # Rename cols
        data.rename(columns={'Raw Data (600 1)': 'PR_Corrected OD600',
                             'Raw Data (584 2)': 'PR_Corrected Red Fluorescence (a.u.)'}, inplace=True)
    
        # Calculate Fluo/OD
        data['PR_Corrected Red Fluo/OD600 (a.u.)'] = data['PR_Corrected Red Fluorescence (a.u.)'] / data['PR_Corrected OD600']
        
        # Add extra metadata based on info in filename
        # Extract information from filename
        fnInfo = matchedfn.split('.xlsx')[0][3:]    #Removes 'PR' from filename for downstream analysis
        run_no = int(fnInfo.split('R')[1][0])
                  
        # Check if induction time and plate no info are in the filename
        if fnInfo.find('PI')>=0:
            ind_time_frag = fnInfo.split('PI')[1]
            if ind_time_frag.find('P') >=0:
                try:
                    plate_no = int(ind_time_frag.split('P')[1])
                except ValueError:
                    plate_no = ind_time_frag.split('P')[1]
            try:
                ind_time = int(ind_time_frag.split('P')[0])
            except ValueError:
                ind_time = ind_time_frag.split('P')[0]
        elif fnInfo.find('P') >=0:
            plate_no = int(fnInfo.split('P')[1])
        
        data['Run'] = run_no
        data['Post-induction (hrs)'] = ind_time
        data['PR_Plate'] = plate_no
    
        # Append metadata
        data = data.merge(metadata_df,left_index=True,right_index=True)
        
        data['PR_Well'] = data.index
        alldata = alldata.append(data,ignore_index=True,sort=None)
    
alldata = alldata[alldata['SampleID'] != 'Blank']

# %%
# Save file as CSV
output_csv_path = os.path.join(root_path, data_folder, output_fn)
alldata.to_csv(output_csv_path)

json_path = output_csv_path.split(".csv")[0] + ".json"
alldata.to_json(json_path)