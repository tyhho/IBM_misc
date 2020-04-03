# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 16:35:21 2020

@author: User
"""
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
    
import CustomFunctions as cf
import fnmatch
import os
import pandas as pd


# mode 1: exported data already contained corrected information
root_path = '..'
expt_folder = 'FC015'
output_fn = 'IBM_FC015R1_PRData.csv'



file_patterns_to_search = ['PR_IBM_FC015R1PI23P1.xlsx']
metadata_fn_core = 'PRMD_IBM_FC015R1'


json_output_fn = 'PR_IBM_FC015R1_PRData.json'

def extract_pr_data_mode1(root_path: str, expt_folder: str,
                      metadata_fn_core: str, file_patterns_to_search: list,
                      json_output_fn: str):
    
    ''' Process exported xlsx files from BMG and save as a json file
    
    All .xlsx files of the same experiment should be placed under the same
    experiment folder. This function processs all those xlsx filesand put them
    into a table.
    
    It then adds metadata to the table by finding the metadata files within
    the same experiment folder. The metadata file is matched according to the
    xlsx filename.
    
    This function automatically saves the final dataframe as a json object
    under the experiment folder.
    
    Parameters:
        root_path (str): the root directory where the experiment folder is,
            abs path (raw string) or relative path
            
        expt_folder (str): name of experiment folder
        
        metadata_fn_core (str): core components of metadata filename, 
            see CustomFunctions.find_meta_xlsx()
        
        file_patterns_to_search (list): a list of file patterns to be 
            searched, each should be a str starting with '.PR_' and ending 
            with '.xlsx'
        
        json_output_fn (str): filename of json, must end in '.json'
                
    Returns: dataframe containing all plate reading data of one experiment
    
    '''
    
    expt_folder_path = os.path.join(root_path, expt_folder)
    files_in_expt_folder = os.listdir(expt_folder_path)
    
    for file_patterns_to_search in file_patterns_to_search:
        matched_fns = fnmatch.filter(files_in_expt_folder, file_patterns_to_search)
        pr_data_all_plates = pd.DataFrame()
        
        for matched_fn in matched_fns:
            
            metadata_fn = cf.find_meta_xlsx(matched_fn, metadata_fn_core)
            metadata_file_path = os.path.join(expt_folder_path, metadata_fn)
            metadata = cf.metadata_to_metadf(metadata_file_path)
            
            pr_xlsx_path = os.path.join(expt_folder_path, matched_fn)
            pr_xlsx_data = pd.read_excel(pr_xlsx_path, sheet_name = 'End point', index_col=0, skiprows=12)
            pr_xlsx_data = pr_xlsx_data.drop(['Content'], axis=1)
            
            # Update the index so it matches the conventional nomenclature
            pr_xlsx_data = cf.renameDfWellIndex(pr_xlsx_data)
            
            # Merge PR data with metadata
            pr_data = pr_xlsx_data.merge(metadata, left_index=True, right_index=True)
            
            # Drop raw data that did not receive blank correction
            pr_data = pr_data.drop(['Raw Data (600 1)','Raw Data (584 2)'], axis=1)
            # Drop blank data
            pr_data = pr_data.query('SampleID != "Blank"')
            # Rename columns
            pr_data.rename(
                columns={
                'Blank corrected based on Raw Data (600 1)': 'PR_Corrected OD600',
                'Blank corrected based on Raw Data (584 2)': 'PR_Corrected Red Fluorescence (a.u.)',
                'RF/OD600': 'PR_Corrected Red Fluo/OD600 (a.u.)'},
                inplace=True)
            
            # Add extra metadata based on info in filename
            
            # Extract information from filename
            
            fn_info = matched_fn.split('.xlsx')[0][3:] # Removes 'PR' from fn
            run_no = int(fn_info.split('R')[1][0])
                      
            # Check if induction time and plate no info are in the filename
            if fn_info.find('PI')>=0:
                ind_time_frag = fn_info.split('PI')[1]
                if ind_time_frag.find('P') >=0:
                    try:
                        plate_no = int(ind_time_frag.split('P')[1])
                    except ValueError:
                        plate_no = ind_time_frag.split('P')[1]
                try:
                    ind_time = int(ind_time_frag.split('P')[0])
                except ValueError:
                    ind_time = ind_time_frag.split('P')[0]
            elif fn_info.find('P') >=0:
                # pass
                plate_no = int(fn_info.split('P')[1])
            
            pr_data.insert(len(pr_data.columns), 'Run', run_no)
            pr_data.insert(len(pr_data.columns), 'Post-induction (hrs)', ind_time)
            
            try:
                pr_data.insert(len(pr_data.columns), 'PR_Plate', plate_no)
            except NameError:
                pass
           
            pr_data.insert(len(pr_data.columns), 'PR_Well', pr_data.index)
            pr_data_all_plates = pr_data_all_plates.append(pr_data, 
                                                           ignore_index=True, 
                                                           sort=None)
        
    json_output_path = os.path.join(root_path, expt_folder, json_output_fn)
    pr_data_all_plates.to_json(json_output_path)
    return pr_data_all_plates
            
pr_data = extract_pr_data_mode1(root_path, expt_folder,
                      metadata_fn_core, file_patterns_to_search,
                      json_output_fn)
#%% 

# Get all files that match the search criteria
 
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

    # Drop raw data that did not receive blank correction
    data = data.drop(['Raw Data (600 1)','Raw Data (584 2)'], axis=1)
    # Drop blank data
    data = data.drop(from_meta_blank_list,axis=0)
    # Rename cols
    data.rename(columns={'Blank corrected based on Raw Data (600 1)': 'PR_Corrected OD600',
                         'Blank corrected based on Raw Data (584 2)': 'PR_Corrected Red Fluorescence (a.u.)',
                         'RF/OD600': 'PR_Corrected Red Fluo/OD600 (a.u.)'}, inplace=True)
    
            # Add extra metadata based on info in filename
        # Extract information from filename
        fn_info = matched_fn.split('.xlsx')[0][3:]    #Removes 'PR' from filename for downstream analysis
        run_no = int(fn_info.split('R')[1][0])
                  
        # Check if induction time and plate no info are in the filename
        if fn_info.find('PI')>=0:
            ind_time_frag = fn_info.split('PI')[1]
            if ind_time_frag.find('P') >=0:
                try:
                    plate_no = int(ind_time_frag.split('P')[1])
                except ValueError:
                    plate_no = ind_time_frag.split('P')[1]
            try:
                ind_time = int(ind_time_frag.split('P')[0])
            except ValueError:
                ind_time = ind_time_frag.split('P')[0]
        elif fn_info.find('P') >=0:
        # pass
        plate_no = int(fn_info.split('P')[1])
    
    data['Run'] = run_no
    data['Post-induction (hrs)'] = ind_time
    
    try:
        data['PR_Plate'] = plate_no
    except NameError:
        pass

#%%
    
fileList = {
            'PR_IBM_FC007R[4,5,7]*P1.xlsx': ['PRMD_IBM_FC007R4P1.xlsx'],
            # 'PR_IBM_FC033R[1,4,5]*P1.xlsx': ['PRMD_IBM_FC033P1.xlsx'],
            # 'PR_IBM_FC033R[1,4,5]*P2.xlsx': ['PRMD_IBM_FC033P2.xlsx'],
            # 'PR_IBM_FC033R[1,4,5]*P3.xlsx': ['PRMD_IBM_FC033P3.xlsx'],
            # 'PR_IBM_FC033R[1,4,5]*P4.xlsx': ['PRMD_IBM_FC033P4.xlsx']
            }



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
        


        data['PR_Well'] = data.index
        alldata = alldata.append(data,ignore_index=True,sort=None)

