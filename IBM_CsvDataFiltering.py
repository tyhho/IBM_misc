# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 17:51:20 2019

@author: Trevor Ho

Reorganize data of median and sd such that the fluorescence of induction are arranged side by side.
Then calculate fold change and return the list of samples to keep. Apply that to the original dataframe.
"""

import os
import pandas as pd

# TODO: Specify folder location

dataRootDir = r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir = 'FC023'
rawDataFN = 'IBM_FC023R3-5_FCmedian&metadata&PRData_InductionProcessed.csv'
rawStatDataFN = 'IBM_FC023R3-5_FCmedianStats&metadata.csv'
newDataFN = 'IBM_FC023R3-5_FCmedian&metadata&PRData_InductionProcessed_FoldChangeFiltered.csv'
newStatDataFN = 'IBM_FC023R3-5_FCmedianStats&metadata_FoldChangeFiltered.csv'

statCsvDir = os.path.join(dataRootDir,dataFolderDir,rawStatDataFN)
raw_stat_df = pd.read_csv(statCsvDir,index_col=0)


control_list = ['IBMc252 + IBMc249',
        'IBMc259 + IBMc249',
        'IBMc259 + IBMc226']

# Remove controls
stat_df = raw_stat_df[~(raw_stat_df['SampleID'].isin(control_list)) &
                      (raw_stat_df['Post-induction (hrs)']==24) ]

inductions = ['no induction','1 mM arabinose','25 μM DAPG','1 mM arabinose + 25 μM DAPG']

# Create base df for future merge
reorgData = stat_df['SampleID'].drop_duplicates().sort_values().to_frame()

# Reorganize df
for induction in inductions:
    stat_df_sub = stat_df[stat_df['Induction']==induction]
    stat_df_sub = stat_df_sub[['SampleID','mean of median fluorescence (a.u.)']]
    stat_df_sub = stat_df_sub.rename(columns={'mean of median fluorescence (a.u.)':induction})
    reorgData = reorgData.merge(stat_df_sub)

#%%
# Calculate fold change

denominator = '1 mM arabinose + 25 μM DAPG'
foldChangeColNames = []

for induction in inductions:
    if induction != denominator:
        colName = induction + ' / ' + denominator
        foldChangeColNames.append(colName)
        reorgData[colName] = reorgData[induction] / reorgData[denominator]

# 
foldChangeThreshold = 2

all_samples_to_drop = []

for col in foldChangeColNames:
    data_to_drop = reorgData[reorgData[col] < foldChangeThreshold]
    samples_to_drop = data_to_drop['SampleID'].to_list()
    all_samples_to_drop = all_samples_to_drop + samples_to_drop

new_stat_df = raw_stat_df[~(raw_stat_df['SampleID'].isin(all_samples_to_drop))]

newStatDataDir = os.path.join(dataRootDir,dataFolderDir,newStatDataFN)
new_stat_df.to_csv(newStatDataDir)

#%% Return to the original source data csv and perform filtering

rawCsvDir = os.path.join(dataRootDir,dataFolderDir,rawDataFN)
raw_df = pd.read_csv(rawCsvDir,index_col=0)
new_df = raw_df[~(raw_df['SampleID'].isin(all_samples_to_drop))]

newDataDir = os.path.join(dataRootDir,dataFolderDir,newDataFN)
new_df.to_csv(newDataDir)