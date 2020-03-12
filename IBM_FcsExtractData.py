# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 21:32:10 2020

@author: Trevor Ho
"""

import pickle
import numpy as np
import glob
import os
import pandas as pd
import IBM_CustomFunctions as cf
from FlowCytometryTools.core.gates import CompositeGate
from FlowCytometryTools import ThresholdGate, FCPlate

# TODO: Specify folder location
    # Each folder must contain only fcs files that end with well location
dataRootDir = r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir = 'EXP012'

# TODO: Specify the source of plate reader data to merge with flow cytometry data
#pr_data_filename = 'IBM_FC021R2-7_Finalized_PRData.csv'

# TODO: Specify the output filename for the combined data
all_doi_filename = 'IBM_EXP012_annotated_cyto_&_PR_data.csv'

# TODO: Specify metadata file core
metadatafnCore = 'PRMD_IBM_EXP012'

#%%
# TODO: Specify folder sequence for processing

# coreSearchSeqList = [
            # 'IBM_FC032R1PI*_FCS',
#            'IBM_FC021R[2,3,4]*P2_FCS',
#            'IBM_FC021R[2-7]*_FCS'
            # ]

# Get all folder that match the search criteria

#%% Core Processing Codes

all_doi_df = pd.DataFrame(columns=[])

# for coreSearchSeq in coreSearchSeqList:
    
#     dirSearchSeq = os.path.join(dataRootDir, dataFolderDir,coreSearchSeq)
#     matchedFolderList = glob.glob(dirSearchSeq)
    
    
#     for matchedFolder in matchedFolderList:
        
#         plateNameCore = matchedFolder.rsplit('\\')[-1].split('_FCS')[0]
        
# Read Metadata Excelfile
# metafilename = cf.findMetaXlsx(plateNameCore + '.',metadatafnCore)
metafilename = "PRMD_IBM_EXP012.xlsx" # for debugging
metadataDir = os.path.join(dataRootDir,dataFolderDir,metafilename)
metadf = cf.metadata_to_metadf(metadataDir)


plateNameCore = 'IBM_EXP012R1' #for debugging
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

#%%
# dfAll = {} 
# dfAll['Count'] = plate.counts()
# dfAll = {**dfAll, **metadata} # merge the data and metadata together

# Save Median values to an excel file
# outputFilenameXlsx = plateNameCore + '_Data.xlsx'
# outputDir = os.path.join(dataRootDir,dataFolderDir,outputFilenameXlsx)
# writer = pd.ExcelWriter(outputDir, engine='xlsxwriter')

# for sheetName, data in dfAll.items():
#     data.to_excel(writer, sheet_name=sheetName)
# writer.save()

# Produces a separate csv file that merges the metadata and the median
    # It also adds Experiment Run data and also Condition as specified by the filename
    # A single csv file will be generated in the end for all experiment runs

# Extract median red fluorescence and counts for each well and save as a dict, data of interest(doi)

cyto_df = pd.DataFrame()
for well,fc_data in plate.items():
    entry = pd.DataFrame(
        {
            'Well': well,
            'fluorescence (a.u.)': [np.array(fc_data.data['RFP2-H'].to_list())],
            'median fluorescence (a.u.)': fc_data.data['RFP2-H'].median(),
            'Count':fc_data.counts
            }
        )
    cyto_df = cyto_df.append(entry, ignore_index=True)


del plate, well, fc_data, entry # delete var "plate" because it is taking too much space

#%%
# Create df from data of interest

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

#%%
# Merge all data from plate reader into cytometer

# #%%
final_df = cyto_df.merge(alldata)
final_df_minus_np = final_df.drop(['fluorescence (a.u.)'], axis=1)
all_doi_dir = os.path.join(dataRootDir,dataFolderDir,all_doi_filename)
final_df_minus_np.to_csv(all_doi_dir)
final_df.to_json((all_doi_dir[:-4]+'.json'))

get_final = pd.read_json((all_doi_dir[:-4]+'.json'))