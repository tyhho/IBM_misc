# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 23:59:32 2020

@author: User
"""
import os
import pandas as pd
import fnmatch
import CustomFunctions as cf

def extract_pr_data_mode2(root_path: str, expt_folder: str,
                      metadata_fn_core: str, file_patterns_to_search: list,
                      blank_plate: int, blank_well: str,
                      json_output_fn: str):
    
    ''' Process exported xlsx files from BMG and save as a json file. Mode 2
    
    Mode 2: 1 or more exported xlsx files only have raw fluo and OD data and
    no blank data on its own. Blank data is supplied elsewhere from another
    plate, which is still within the list of files processed.
    
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
    
    pr_data_all_plates = pd.DataFrame()
        
    for file_patterns_to_search in file_patterns_to_search:
        matched_fns = fnmatch.filter(files_in_expt_folder, file_patterns_to_search)
            
        for matched_fn in matched_fns:
            
            # Get metadata
            metadata_fn = cf.find_meta_xlsx(matched_fn, metadata_fn_core)
            metadata_file_path = os.path.join(expt_folder_path, metadata_fn)
            metadata = cf.metadata_to_metadf(metadata_file_path)
            
            # Get corresponding blank of the same expt and the same induction time
            blank_data_fn = cf.findBlankXlsx(matched_fn, blank_plate)
            blank_data_file_path = os.path.join(root_path, expt_folder, blank_data_fn)
            blank_data = pd.read_excel(blank_data_file_path,
                                       sheet_name = 'End point',
                                       index_col=0, skiprows=12)
            
            blank_data = blank_data.drop(['Content'], axis=1)
            blank_ID = blank_data.index.get_loc(blank_well)
            blank_OD = blank_data.iloc[blank_ID][0]
            blank_fluo = blank_data.iloc[blank_ID][1]
            
            # Get plate reade
            pr_xlsx_path = os.path.join(expt_folder_path, matched_fn)
            pr_xlsx_data = pd.read_excel(pr_xlsx_path,
                                         sheet_name = 'End point',
                                         index_col=0, skiprows=12)
            pr_xlsx_data = pr_xlsx_data.drop(['Content'], axis=1)
            
            # Update the index so it matches the conventional nomenclature
            pr_xlsx_data = cf.renameDfWellIndex(pr_xlsx_data)
            
            # Blank correction and calculate fluo/OD600
            pr_xlsx_data.insert( len(pr_xlsx_data.columns), 
                                'PR_corrected OD600',
                                pr_xlsx_data['Raw Data (600 1)'] - blank_OD)
            pr_xlsx_data.insert( len(pr_xlsx_data.columns), 
                                'PR_corrected fluorescence (a.u.)',
                                pr_xlsx_data['Raw Data (584 2)'] - blank_fluo)
            pr_xlsx_data.insert( len(pr_xlsx_data.columns), 
                                'PR_corrected fluorescence/OD600 (a.u.)',
                                pr_xlsx_data['PR_corrected fluorescence (a.u.)'] / \
                                pr_xlsx_data['PR_corrected OD600']
                                )
                
            # Merge PR data with metadata
            pr_data = pr_xlsx_data.merge(metadata, left_index=True, right_index=True)
            
            # Drop raw data that did not receive blank correction
            pr_data = pr_data.drop(['Raw Data (600 1)','Raw Data (584 2)'], axis=1)
            
            # Drop blank data
            pr_data = pr_data.query('`sample ID` != "Blank"')
            
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
            
            pr_data.insert(len(pr_data.columns), 'run', run_no)
            pr_data.insert(len(pr_data.columns), 'post-induction (hrs)', ind_time)
            
            try:
                pr_data.insert(len(pr_data.columns), 'PR_plate', plate_no)
            except NameError:
                pass
           
            pr_data.insert(len(pr_data.columns), 'PR_well', pr_data.index)
            pr_data_all_plates = pr_data_all_plates.append(pr_data, 
                                                           ignore_index=True, 
                                                           sort=None)
        
    json_output_path = os.path.join(root_path, expt_folder, json_output_fn)
    pr_data_all_plates.to_json(json_output_path)
    return pr_data_all_plates