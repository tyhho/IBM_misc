# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 17:31:09 2019

@author: Trevor Ho
"""
#try:
#    import cPickle as pickle
#except ModuleNotFoundError:
#    import pickle

import os
from FlowCytometryTools.core.gates import CompositeGate
from FlowCytometryTools import FCPlate, ThresholdGate
import pandas as pd
import IBM_CustomFunctions as cf

# TODO: Specify folder location
    # Each folder must contain only fcs files that end with well location
dataRootDir = r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir = 'FC022'

# TODO: Specify the source of plate reader data to merge with flow cytometry data
pr_data_filename = 'IBM_FC022R1_PRData.csv'

# TODO: Specify the output filename for the combined data
all_doi_filename = 'IBM_FC022R1_FCmedian&metadata&PRData.csv'

# TODO: Specify subfolders containing FCS files
# Create dict with information of FCS folder name (key) and also Metadata file (value)

# TODO: In the future, this needs to be done automatically
plateList = {'IBM_FC018R5PI5P1':'FCMD_IBM_FC018P1',
             'IBM_FC018R5PI5P2':'PRMD_IBM_FC018P2',
             'IBM_FC018R5PI5P3':'PRMD_IBM_FC018P3',
             'IBM_FC018R5PI5P4':'PRMD_IBM_FC018P4',
             'IBM_FC018R5PI5P1':'FCMD_IBM_FC018P1',
             'IBM_FC018R5PI24P2':'PRMD_IBM_FC018P2',
             'IBM_FC018R5PI24P3':'PRMD_IBM_FC018P3',
             'IBM_FC018R5PI24P4':'PRMD_IBM_FC018P4',
             'IBM_FC018R6PI5P1':'FCMD_IBM_FC018P1',
             'IBM_FC018R6PI5P2':'PRMD_IBM_FC018P2',
             'IBM_FC018R6PI5P3':'PRMD_IBM_FC018P3',
             'IBM_FC018R6PI5P4':'PRMD_IBM_FC018P4',
             'IBM_FC018R6PI5P1':'FCMD_IBM_FC018P1',
             'IBM_FC018R6PI24P2':'PRMD_IBM_FC018P2',
             'IBM_FC018R6PI24P3':'PRMD_IBM_FC018P3',
             'IBM_FC018R6PI24P4':'PRMD_IBM_FC018P4',
             'IBM_FC018R7PI5P1':'FCMD_IBM_FC018P1',
             'IBM_FC018R7PI5P2':'PRMD_IBM_FC018P2',
             'IBM_FC018R7PI5P3':'PRMD_IBM_FC018P3',
             'IBM_FC018R7PI5P4':'PRMD_IBM_FC018P4',
             'IBM_FC018R7PI5P1':'FCMD_IBM_FC018P1',
             'IBM_FC018R7PI24P2':'PRMD_IBM_FC018P2',
             'IBM_FC018R7PI24P3':'PRMD_IBM_FC018P3',
             'IBM_FC018R7PI24P4':'PRMD_IBM_FC018P4'
             }

#%% Core Processing Codes

all_doi_df = pd.DataFrame(columns=[])

for plateNameCore,metadataNameCore in plateList.items():
    # Read Metadata Excelfile
    metafilename = metadataNameCore + '.xlsx'
    metadataDir = os.path.join(dataRootDir,dataFolderDir,metafilename)
    metaxls= pd.ExcelFile(metadataDir)
    metadata = {sheet:metaxls.parse(sheet, index_col=0) for sheet in metaxls.sheet_names}    #import all sheets in metadata file into a dict, property name=keys, metadata df = values
    
    #plateNameCore = '' #for debugging
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

    # Produces 96 well layout for intuitive observation
    # Calculate Median from data
    def calculate_median_Y2(well):
        return well.data['RFP2-H'].median()
    
    dfAll = {} 
    dfAll['RFP'] = plate.apply(calculate_median_Y2)
    dfAll['Count'] = plate.counts()
    dfAll = {**dfAll, **metadata} # merge the data and metadata together
    
    # Save Median values to an excel file
    outputFilenameXlsx = plateNameCore + '_Data.xlsx'
    outputDir = os.path.join(dataRootDir,dataFolderDir,outputFilenameXlsx)
    writer = pd.ExcelWriter(outputDir, engine='xlsxwriter')
    
    for sheetName, data in dfAll.items():
        data.to_excel(writer, sheet_name=sheetName)
    writer.save()
    
    #%% Produces a separate csv file that merges the metadata and the median
        # It also adds Experiment Run data and also Condition as specified by the filename
        # A single csv file will be generated in the end for all experiment runs

    # Extract median red fluorescence and counts for each well and save as a dict, data of interest(doi)
    doi_dict = {well[0]: [well[1].data['RFP2-H'].median(),well[1].counts] \
                                  for well in plate.items()}    
    
    del plate # delete var "plate" because it is taking too much space
    
    # Create df from data of interest
    doi_df = pd.DataFrame.from_dict(doi_dict,orient='index',columns=['median fluorescence (a.u.)','Count'])
    doi_df['Run'] = int(plateNameCore.split('R')[1].split('PI')[0])  # add run number to df
    post_induction_time = plateNameCore.split('PI')[1] .split('P')[0]
    
    # In the event that post_induction_time is not an integer
    try:
        doi_df['Post-induction (hrs)'] = int(post_induction_time)
    except ValueError:
        doi_df['Post-induction (hrs)'] = post_induction_time
    
    # TODO: need to add method to add FC_Plate number to dataset
    
    # Process the df of metadata & merge into the main dataframe
    for meta_property, metadf_96format in metadata.items():
        metadf = pd.DataFrame(columns = [])
        for char in ['A','B','C','D','E','F','G','H']:
            metadf_96_row = metadf_96format.transpose()[char]
            metadf_96_row = metadf_96_row.to_frame(meta_property)
            indexList = metadf_96_row.index.tolist()
            
            # Retrieve information of the wells and set it as index
            for columnIndex in range(len(indexList)):
                indexList[columnIndex] = char + str(indexList[columnIndex])
            metadf_96_row['Well']=indexList
            metadf_96_row.set_index('Well', inplace=True)
            
            # Append to major dataframe
            metadf = metadf.append(metadf_96_row,sort=False)
            
        doi_df = doi_df.merge(metadf,left_index=True,right_index=True)
    doi_df['FC_Well'] = doi_df.index
    all_doi_df = all_doi_df.append(doi_df,ignore_index=True,sort=False)
    
    #%%
# Merge all data from plate reader into cytometer
pr_data_dir = os.path.join(dataRootDir,dataFolderDir,pr_data_filename)
all_pr_data = pd.read_csv(pr_data_dir,index_col=0)
final_all_doi_df = all_doi_df.merge(all_pr_data)

all_doi_dir = os.path.join(dataRootDir,dataFolderDir,all_doi_filename)
final_all_doi_df.to_csv(all_doi_dir)