# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 21:32:10 2020

@author: Trevor Ho
"""

import numpy as np
import glob
import os
import pandas as pd
import IBM_CustomFunctions as cf
from FlowCytometryTools.core.gates import CompositeGate
from FlowCytometryTools import ThresholdGate, FCPlate

# TODO: Specify folder location
    # Each folder must contain only fcs files that end with well location
dataRootDir = '..'
dataFolderDir = 'FC007'

# TODO: Specify the source of plate reader data to merge with flow cytometry data
pr_data_filename = 'IBM_FC007R1-3_PRData.json'

# TODO: Specify the output filename for the combined data
json_output_fn = 'IBM_FC007R1-3_cyto_&_PR_data.json'

# TODO: Specify metadata file core
metadatafnCore = 'PRMD_IBM_FC007'

#%%
# TODO: Specify folder sequence for processing

coreSearchSeqList = [
            # 'IBM_FC007R[4,5,7]*_FCS',
            'IBM_FC007R[1-3]*_FCS',
            # 'IBM_FC021R[2,3,4]*P2_FCS',
            # 'IBM_FC021R[2-7]*_FCS'
            ]

# Get all folder that match the search criteria


#%%

for coreSearchSeq in coreSearchSeqList:
    
    dirSearchSeq = os.path.join(dataRootDir, dataFolderDir,coreSearchSeq)
    matchedFolderList = glob.glob(dirSearchSeq)
    
    for matchedFolder in matchedFolderList:
        
        
        plateNameCore = matchedFolder.rsplit('\\')[-1].split('_FCS')[0]
        metafilename = cf.findMetaXlsx(plateNameCore + '.', metadatafnCore)
        print(metafilename)



#%% Core Processing Codes

final_df = pd.DataFrame()

for coreSearchSeq in coreSearchSeqList:
    
    dirSearchSeq = os.path.join(dataRootDir, dataFolderDir,coreSearchSeq)
    matchedFolderList = glob.glob(dirSearchSeq)
    
    for matchedFolder in matchedFolderList:
        
        plateNameCore = matchedFolder.rsplit('\\')[-1].split('_FCS')[0]
                
        # Read Metadata Excelfile
        metafilename = cf.findMetaXlsx(plateNameCore + '.', metadatafnCore)
        # metafilename = "PRMD_IBM_EXP012.xlsx" # for debugging
        metadataDir = os.path.join(dataRootDir,dataFolderDir,metafilename)
        metadf = cf.metadata_to_metadf(metadataDir)

        
        # plateNameCore = 'IBM_EXP012R1' #for debugging
        fcsFolderDir = plateNameCore + '_FCS'
        datadir = os.path.join(dataRootDir,dataFolderDir,fcsFolderDir)
        plate = FCPlate.from_dir(ID='Plate', path=datadir, parser=cf.fcsNameParser, position_mapper='name')
        
        # Gating
        fsc_gate = ThresholdGate(1000.0, ['FSC-H'], region='above')
        ssc_gate = ThresholdGate(1000.0, ['SSC-H'], region='above')
        
        rfpA_gate = ThresholdGate(1.0, ['RFP2-A'], region='above')
        fscssc_gate = CompositeGate(fsc_gate,'and',ssc_gate)
        plate = plate.gate(fscssc_gate)
        plate = plate.gate(rfpA_gate)
        
        cyto_df = pd.DataFrame()
        for well,fc_data in plate.items():
            entry = pd.DataFrame(
                {
                    'Well': well,
                    'FC_fluorescence (a.u.)': [np.array(fc_data.data['RFP2-H'].to_list())],
                    'FC_median fluorescence (a.u.)': fc_data.data['RFP2-H'].median(),
                    'FC_Count':fc_data.counts
                    }
                )
            cyto_df = cyto_df.append(entry, ignore_index=True)
        
        
        del plate, well, fc_data, entry # delete var "plate" because it is taking too much space
        
        cyto_df['Run'] = int(plateNameCore.split('R')[1].split('PI')[0])  # add run number to df
        
        if "PI" in plateNameCore:
            post_induction_time = plateNameCore.split('PI')[1].split('P')[0]
        
            # In the event that post_induction_time is not an integer
            try:
                cyto_df['Post-induction (hrs)'] = int(post_induction_time)
            except ValueError:
                cyto_df['Post-induction (hrs)'] = post_induction_time
        
        # Add FC_Plate number to dataset, if it is present
        try:
            fc_plate_no = plateNameCore.split('PI')[1] .split('P')[1]
            cyto_df['FC_Plate'] = int(fc_plate_no)
        except IndexError or ValueError:
            pass
        
        cyto_df = cyto_df.merge(metadf, left_on= 'Well' ,right_index=True)
        
        final_df = final_df.append(cyto_df, ignore_index=True)

#%%
# Merge all data from plate reader into cytometer

pr_data_path = os.path.join(dataRootDir, dataFolderDir, pr_data_filename)
pr_data = pd.read_json(pr_data_path)

#%%

# #%%
final_df = final_df.merge(pr_data)
json_output_path = os.path.join(dataRootDir, dataFolderDir, json_output_fn)
final_df.to_json(json_output_path)
csv_data = final_df.drop(['FC_fluorescence (a.u.)'], axis=1)
csv_output_path = json_output_path.split('.json')[0] + '.csv'
csv_data.to_csv(csv_output_path)