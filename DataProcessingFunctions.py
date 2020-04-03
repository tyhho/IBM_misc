# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 15:11:47 2020

@author: Trevor Ho
"""

import fnmatch
import glob
import os

import numpy as np
import pandas as pd
import CustomFunctions as cf
from scipy import stats

from FlowCytometryTools.core.gates import CompositeGate
from FlowCytometryTools import IntervalGate, FCPlate


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
            pr_data = pr_data.query('`sample ID` != "Blank"')
            # Rename columns
            pr_data.rename(
                columns={
                'Blank corrected based on Raw Data (600 1)': 'PR_corrected OD600',
                'Blank corrected based on Raw Data (584 2)': 'PR_corrected Red Fluorescence (a.u.)',
                'RF/OD600': 'PR_corrected Red Fluo/OD600 (a.u.)'},
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

def extract_fcs_data(root_path: str, expt_folder: str,
                     metadata_fn_core: str, folder_patterns_to_search: list,
                     json_output_fn: str):
    
    ''' Process .fcs files, gate, extract RFP2-H values and save as a dataframe
    
    All .fcs files of the same plate should be placed under a folder ending
    with '_FCS' in an experiment folder. An experiment folder can contain
    multiple fcs plate folders. This function processs all those folders and
    their fcs files inside by gating them, extracting their single-cell fluo
    values and then putting them into a table.
    
    It then adds metadata to the table by finding the metadata files within
    the same experiment folder. The metadata file is matched according to the
    plate folder name.
    
    This function automatically saves the final dataframe as a json object
    under the experiment folder.
    
    Parameters:
        root_path (str): the root directory where the experiment folder is,
            abs path (raw string) or relative path
            
        expt_folder (str): name of experiment folder
        
        json_output_fn (str): filename of json, must end in '.json'
        
        metadata_fn_core (str): core components of metadata filename, 
            see CustomFunctions.find_meta_xlsx()
            
        folder_patterns_to_search (list): a list of folder patterns to be 
            searched, each should be a str ending with '_FCS'
        
    Returns: dataframe containing all cytometry data of one experiment

    '''
    #Core Processing Codes
    
    cyto_data = pd.DataFrame()
    
    for folder_pattern_to_search in folder_patterns_to_search:
        
        directories_to_search = os.path.join(root_path, expt_folder,folder_pattern_to_search)
        matched_folders = glob.glob(directories_to_search)
        
        for matched_folder in matched_folders:
            
            plate_name_core = matched_folder.rsplit('\\')[-1].split('_FCS')[0]
                    
            # Read Metadata Excelfile
            metadata_fn = cf.find_meta_xlsx(plate_name_core + '.', metadata_fn_core)
            # metadata_fn = "PRMD_IBM_EXP012.xlsx" # for debugging
            metadata_path = os.path.join(root_path, expt_folder, metadata_fn)
            metadf = cf.metadata_to_metadf(metadata_path)
    
            
            # plate_name_core = 'IBM_EXP012R1' #for debugging
            fcs_folder = plate_name_core + '_FCS'
            datadir = os.path.join(root_path,expt_folder,fcs_folder)
            plate = FCPlate.from_dir(ID='Plate', path=datadir, parser=cf.fcsNameParser, position_mapper='name')
            
            # Gating
            fsc_gate = IntervalGate((1e3,1e5), ['FSC-H'], region='in')
            ssc_gate = IntervalGate((1e3,1e5), ['SSC-H'], region='in')
            
    
            rfpA_gate = IntervalGate((1,1e6), ['RFP2-A'], region='in')
            rfpH_gate = IntervalGate((1,1e6), ['RFP2-H'], region='in')
            
            fsc_ssc_gate = CompositeGate(fsc_gate,'and',ssc_gate)
            rfp_AH_gate = CompositeGate(rfpA_gate,'and',rfpH_gate)
            plate = plate.gate(fsc_ssc_gate)
            plate = plate.gate(rfp_AH_gate)
            
            cyto_data_entries = pd.DataFrame()
            for well,fc_data in plate.items():
                entry = pd.DataFrame(
                    {
                        'well': well,
                        'FC_fluorescence (a.u.)': [np.array(fc_data.data['RFP2-H'].to_list())],
                        'FC_median fluorescence (a.u.)': fc_data.data['RFP2-H'].median(),
                        'FC_Count':fc_data.counts
                        }
                    )
                cyto_data_entries = cyto_data_entries.append(entry, ignore_index=True)
            
            del plate, well, fc_data, entry # delete var "plate" because it is taking too much space
            
            # Infer metadata from code
            cyto_data_entries['Run'] = int(plate_name_core.split('R')[1].split('PI')[0])  # add run number to df
            
            if "PI" in plate_name_core:
                post_induction_time = plate_name_core.split('PI')[1].split('P')[0]
            
                # In the event that post_induction_time is not an integer
                try:
                    cyto_data_entries['post-induction (hrs)'] = int(post_induction_time)
                except ValueError:
                    cyto_data_entries['post-induction (hrs)'] = post_induction_time
            
            # Add FC_Plate number to dataset, if it is present
            try:
                fc_plate_no = plate_name_core.split('PI')[1] .split('P')[1]
                cyto_data_entries['FC_plate'] = int(fc_plate_no)
            except IndexError or ValueError:
                pass
            
            cyto_data_entries = cyto_data_entries.merge(metadf, left_on= 'well' ,right_index=True)
            cyto_data = cyto_data.append(cyto_data_entries, ignore_index=True)
    
    json_output_path = os.path.join(root_path, expt_folder, json_output_fn)
    cyto_data.to_json(json_output_path)
    return cyto_data
    # csv_data = cyto_data.drop(['FC_fluorescence (a.u.)'], axis=1)
    # csv_output_path = json_output_path.split('.json')[0] + '.csv'
    # csv_data.to_csv(csv_output_path)

def merge_cyto_and_pr_data(root_path: str, expt_folder: str,
                           pr_data_fn: str, cyto_data_fn: str,
                           json_output_fn: str):
    
    ''' Merge plate reader and cytometry data and save as new json file
    
    This function automatically saves the final dataframe as a json object
    and as a csv file but wihtout single-cell fluo under the experiment folder.
    
    Parameters:
        root_path (str): the root directory where the experiment folder is,
            abs path (raw string) or relative path
            
        expt_folder (str): name of experiment folder
        
        pr_input_fn (str): filename of plate reader json file, must end in '.json'
        
        cyto_input_fn (str): filename of plate reader json file, must end in '.json'

        json_output_fn (str): filename of output json, must end in '.json'
        
        
    Returns: dataframe containing all cytometry and plate reader data of one experiment

    '''

    pr_data_path = os.path.join(root_path, expt_folder, pr_data_fn)
    pr_data = pd.read_json(pr_data_path)
    
    cyto_data_path = os.path.join(root_path, expt_folder, cyto_data_fn)
    cyto_data = pd.read_json(cyto_data_path)
    
    # Export data into both json and 
    merged_data = cyto_data.merge(pr_data)

    json_output_path = os.path.join(root_path, expt_folder, json_output_fn)
    merged_data.to_json(json_output_path)
    csv_data = merged_data.drop(['FC_fluorescence (a.u.)'], axis=1)
    csv_output_path = json_output_path.split('.json')[0] + '.csv'
    csv_data.to_csv(csv_output_path)
    
    return merged_data


#%%

def json_cal_stats(root_path: str, expt_folder: str,
                   json_data_fn: str):
    
    ''' Caculate statics for cytometer data from a json file and save as new json file
    
    This function automatically saves the final dataframe as a json object
    and as a csv file but wihtout single-cell fluo under the experiment folder.
    The filenames will end as '_Stats.json' or '_Stats.csv'
    
    This function:
        1. Combine all single-cell fluo data and obtain a single median and rSD
        2. Calcualtes mean and std of median fluo from multiple runs
    
    Parameters:
        root_path (str): the root directory where the experiment folder is,
            abs path (raw string) or relative path
            
        expt_folder (str): name of experiment folder
        
        json_data_fn (str): filename of in json file, must end in '.json'
        This file should be generated from merged cytometer and plate reader files
        
    Returns: dataframe containing calculated statiscs and concatanated single-cell
        fluorescence from one experiment

    '''

    # Read data
    json_data_path = os.path.join(root_path, expt_folder, json_data_fn)
    json_data = pd.read_json(json_data_path)
    
    # Generate deduplicated sample_ID, pi_times and inductions
    sample_IDs = json_data['sample ID'].drop_duplicates().to_list()
    inductions = json_data['induction'].drop_duplicates().to_list()
    pi_times = json_data['post-induction (hrs)'].drop_duplicates().to_list()
    query_line = '`sample ID` == @sample_ID and induction == @induction and `post-induction (hrs)` == @pi_time' 
    
    ## Generate data for mean and sd, according to the induction
    
    stat_data = pd.DataFrame(columns=[])
    for pi_time in pi_times:
        for induction in inductions:
            for sample_ID in sample_IDs:
                sub_data = json_data.query(query_line)
                
                fluo = []
                for single_sample_fluo in sub_data['FC_fluorescence (a.u.)'].to_list():
                    fluo += single_sample_fluo
                fluo = np.array(fluo)
                
                stat_entry = pd.DataFrame(
                    {
                    'sample ID': sample_ID,
                    'induction': induction,
                    'post-induction (hrs)': pi_time,
                    'fluorescence (a.u.)': [fluo],
                    'mean of median fluorescence (a.u.)': sub_data['FC_median fluorescence (a.u.)'].mean(),
                    'std of median fluorescence (a.u.)': sub_data['FC_median fluorescence (a.u.)'].std(),
                    'median fluorescence (a.u.)': np.median(fluo),
                    'rSD (a.u.)':stats.median_absolute_deviation(fluo)
                    })
                
                stat_data = stat_data.append(stat_entry, ignore_index=True)
                
    # Export statiscs to json and csv file
    stat_data_path = json_data_path.split('.json')[0] + '_Stats.json'
    stat_data.to_json(stat_data_path)
    
    csv_stat_data = stat_data.drop(['fluorescence (a.u.)'], axis=1)
    csv_output_path = stat_data_path.split('.json')[0] + '.csv'
    csv_stat_data.to_csv(csv_output_path)

    return stat_data