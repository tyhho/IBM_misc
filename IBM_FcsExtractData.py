# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 15:11:47 2020

@author: Trevor Ho
"""

import numpy as np
import glob
import os
import pandas as pd
import CustomFunctions as cf
from FlowCytometryTools.core.gates import CompositeGate
from FlowCytometryTools import IntervalGate, FCPlate

# TODO: Specify folder location
    # Each folder must contain only fcs files that end with well location
root_path = '..'
expt_folder = 'FC007'

# TODO: Specify the source of plate reader data to merge with flow cytometry data
# pr_data_filename = 'IBM_FC007R4,5,7_PRData.json'

pr_data_filename = 'IBM_FC007R1-3_PRData.json'

# TODO: Specify the output filename for the combined data
json_output_fn = 'IBM_FC007R1-3_cyto_&_PR_data.json'

# TODO: Specify metadata file core
metadata_fn_core = 'PRMD_IBM_FC007'

#%%
# TODO: Specify folder sequence for processing

folder_patterns_to_search = [
            # 'IBM_FC007R[4,5,7]*_FCS',
            'IBM_FC007R[1-3]*_FCS',
            # 'IBM_FC021R[2,3,4]*P2_FCS',
            # 'IBM_FC021R[2-7]*_FCS'
            ]

# Get all folder that match the search criteria

#%%

# for folder_pattern_to_search in folder_patterns_to_search:
    
#     directories_to_search = os.path.join(root_path, expt_folder,folder_pattern_to_search)
#     matched_folders = glob.glob(directories_to_search)

#     for matched_folder in matched_folders:
        
#         plate_name_core = matched_folder.rsplit('\\')[-1].split('_FCS')[0]
#         metadata_fn = cf.findMetaXlsx(plate_name_core + '.', metadata_fn_core)
#         print(metadata_fn)

#%% Core Processing Codes

cyto_df = pd.DataFrame()

for folder_pattern_to_search in folder_patterns_to_search:
    
    directories_to_search = os.path.join(root_path, expt_folder,folder_pattern_to_search)
    matched_folders = glob.glob(directories_to_search)
    
    for matched_folder in matched_folders:
        
        plate_name_core = matched_folder.rsplit('\\')[-1].split('_FCS')[0]
                
        # Read Metadata Excelfile
        metadata_fn = cf.findMetaXlsx(plate_name_core + '.', metadata_fn_core)
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
        
        cyto_df_entries = pd.DataFrame()
        for well,fc_data in plate.items():
            entry = pd.DataFrame(
                {
                    'Well': well,
                    'FC_fluorescence (a.u.)': [np.array(fc_data.data['RFP2-H'].to_list())],
                    'FC_median fluorescence (a.u.)': fc_data.data['RFP2-H'].median(),
                    'FC_Count':fc_data.counts
                    }
                )
            cyto_df_entries = cyto_df_entries.append(entry, ignore_index=True)
        
        del plate, well, fc_data, entry # delete var "plate" because it is taking too much space
        
        cyto_df_entries['Run'] = int(plate_name_core.split('R')[1].split('PI')[0])  # add run number to df
        
        if "PI" in plate_name_core:
            post_induction_time = plate_name_core.split('PI')[1].split('P')[0]
        
            # In the event that post_induction_time is not an integer
            try:
                cyto_df_entries['Post-induction (hrs)'] = int(post_induction_time)
            except ValueError:
                cyto_df_entries['Post-induction (hrs)'] = post_induction_time
        
        # Add FC_Plate number to dataset, if it is present
        try:
            fc_plate_no = plate_name_core.split('PI')[1] .split('P')[1]
            cyto_df_entries['FC_Plate'] = int(fc_plate_no)
        except IndexError or ValueError:
            pass
        
        cyto_df_entries = cyto_df_entries.merge(metadf, left_on= 'Well' ,right_index=True)
        
        cyto_df = cyto_df.append(cyto_df_entries, ignore_index=True)

#%%
# Merge all data from plate reader into cytometer

pr_data_path = os.path.join(root_path, expt_folder, pr_data_filename)
pr_data = pd.read_json(pr_data_path)

# Export data into both json and 
cyto_df = cyto_df.merge(pr_data)
json_output_path = os.path.join(root_path, expt_folder, json_output_fn)
cyto_df.to_json(json_output_path)
csv_data = cyto_df.drop(['FC_fluorescence (a.u.)'], axis=1)
csv_output_path = json_output_path.split('.json')[0] + '.csv'
csv_data.to_csv(csv_output_path)