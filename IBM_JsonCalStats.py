# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 00:07:32 2020

@author: Trevor Ho
"""


import os
import pandas as pd
import json
import numpy as np
from  scipy import stats

# TODO: Specify folder location
    # csv file containing data to be plotted should be in this folder

root_path = '..'
data_folder = 'FC007'
json_data_fn = 'IBM_FC007R1-5,7_cyto_&_PR_data_InductionRelabelled.json'
stat_output_fn = json_data_fn.split('.json')[0] + '_Stats.json'

# TODO: Set induction and induction time information
# inductions = ['+ DMSO','+ 10 μM 4-HT']
# inductions = ['no induction','1 mM arabinose','25 μM DAPG','1 mM arabinose + 25 μM DAPG']
# inductions = ['no induction','1 mM arabinose','1 μM AHL','1 mM arabinose + 1 μM AHL']
# pi_times = [5,24]

# Read data
json_data_path = os.path.join(root_path, data_folder, json_data_fn)
json_data = pd.read_json(json_data_path)

#%%


# Generate deduplicated SampleID, pi_times and inductions
sampleIDs = json_data['SampleID'].drop_duplicates().to_list()
inductions = json_data['Induction'].drop_duplicates().to_list()
pi_times = json_data['Post-induction (hrs)'].drop_duplicates().to_list()
query_line = 'SampleID == @sampleID and Induction == @induction and `Post-induction (hrs)` == @pi_time' 

## Generate data for mean and sd, according to the induction

# sampleIDs = ['IBMc052']
# induction = ['no induction']
# pi_time = [5, 24]

#%%
stat_data = pd.DataFrame(columns=[])
for pi_time in pi_times:
    for induction in inductions:
        for sampleID in sampleIDs:
            sub_data = json_data.query(query_line)
            
            fluo = []
            for single_sample_fluo in sub_data['FC_fluorescence (a.u.)'].to_list():
                fluo += single_sample_fluo
            fluo = np.array(fluo)
            
            stat_entry = pd.DataFrame(
                {
                'SampleID': sampleID,
                'Induction': induction,
                'Post-induction (hrs)': pi_time,
                'fluorescence (a.u.)': [fluo],
                'mean of median fluorescence (a.u.)': sub_data['FC_median fluorescence (a.u.)'].mean(),
                'std of median fluorescence (a.u.)': sub_data['FC_median fluorescence (a.u.)'].std(),
                'median fluorescence (a.u.)': np.median(fluo),
                'rSD (a.u.)':stats.median_absolute_deviation(fluo)
                })
            
            stat_data = stat_data.append(stat_entry, ignore_index=True)

#%% Export statiscs to file
stat_data_path = os.path.join(root_path, data_folder, stat_output_fn)
stat_data.to_json(stat_data_path)

csv_stat_data = stat_data.drop(['fluorescence (a.u.)'], axis=1)
csv_output_path = stat_data_path.split('.json')[0] + '.csv'
csv_stat_data.to_csv(csv_output_path)
#%% Filter out data points where SD were too large

# Single out criteria for filtering

pi24_data = stat_data[(stat_data['Post-induction (hrs)'] ==24) & (stat_data['Induction'] == '1 mM arabinose + 25 μM DAPG')]
pi24_data = pi24_data[pi24_data['SampleID'].str.contains("IBMs")]
pi24_data['cv'] = pi24_data.loc[:,('std of median fluorescence (a.u.)')] / pi24_data.loc[:,('mean of median fluorescence (a.u.)')]

fc_data_induced = stat_data[(stat_data['Post-induction (hrs)'] ==24) & (stat_data['Induction'] == '1 mM arabinose + 25 μM DAPG')]
fc_data_induced = fc_data_induced[['SampleID','mean of median fluorescence (a.u.)']]
fc_data_induced.rename(columns={'mean of median fluorescence (a.u.)':'induced'}, inplace=True)

fc_data_uninduced = stat_data[(stat_data['Post-induction (hrs)'] ==24) & (stat_data['Induction'] == 'no induction')]
fc_data_uninduced = fc_data_uninduced[['SampleID','mean of median fluorescence (a.u.)']]
fc_data_uninduced.rename(columns={'mean of median fluorescence (a.u.)':'uninduced'}, inplace=True)

fc_data_withCtrl = fc_data_induced.merge(fc_data_uninduced, left_on='SampleID', right_on='SampleID')
fc_data = fc_data_withCtrl[fc_data_withCtrl['SampleID'].str.contains("IBMs")]
fc_data['fold_change'] = fc_data.loc[:,('uninduced')] / fc_data.loc[:,('induced')]

fold_change_threshold = 2
cv_threshold = 0.3

#%%
data_to_drop_cv = pi24_data.query(f'cv > {cv_threshold}')

data_to_drop_fc = fc_data.query(f'fold_change < {fold_change_threshold}')

samples_to_drop = data_to_drop_cv['SampleID'].to_list() + data_to_drop_fc['SampleID'].to_list()

#%%
filtered_stat_data =  stat_data[~(stat_data['SampleID'].isin(samples_to_drop))]
filtered_stat_data.reset_index(drop=True, inplace=True)

#%% Export filtered statiscs to file
filteredStatFN = dataFN.split('.csv')[0] + '_FilteredStats.csv'
filteredOutputDataFP = os.path.join(dataRootDir,dataFolderDir,filteredStatFN)
filtered_stat_data.to_csv(filteredOutputDataFP)

#%% Remerge with fluo data
filtered_fluo_data = fluo_data[~(fluo_data['SampleID'].isin(samples_to_drop))]
filteredFluoFN = dataFN.split('.csv')[0] + '_Filtered.csv'
filteredFluoFP = os.path.join(dataRootDir,dataFolderDir,filteredFluoFN)
filtered_fluo_data.to_csv(filteredFluoFP)