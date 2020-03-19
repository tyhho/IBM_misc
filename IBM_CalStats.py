# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 17:14:54 2020

@author: Trevor Ho
"""

import os
import pandas as pd


# TODO: Specify folder location
    # csv file containing data to be plotted should be in this folder

dataRootDir=r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir='FC034'
dataFN = 'IBM_FC034R1-4_median&metadata&PRData_InductionRelabelled.csv'
statFN = dataFN.split('.csv')[0] + '_Stats.csv'

# TODO: Set induction and induction time information
# inductions = ['+ DMSO','+ 10 μM 4-HT']
inductions = ['no induction','1 mM arabinose','25 μM DAPG','1 mM arabinose + 25 μM DAPG']
pi_times = [5,24]

# Read data
dataFP = os.path.join(dataRootDir,dataFolderDir,dataFN)
fluo_data = pd.read_csv(dataFP,index_col=0)

## Generate data for mean and sd, according to the induction
stat_data = pd.DataFrame(columns=[])
for pi_time in pi_times:
    for induction in inductions:
        subdata = fluo_data[(fluo_data['Induction']==induction) & (fluo_data['Post-induction (hrs)']==pi_time)]
        subdata.set_index(['SampleID'],inplace=True)
        for index in set(subdata.index):
            subframe = subdata.loc[index]
            medianFluos = subframe['median fluorescence (a.u.)']
            mean_medianFluo = medianFluos.mean()
            sd_medianFluo = medianFluos.std()
            rowInfo = pd.DataFrame([[index,mean_medianFluo,sd_medianFluo,induction,pi_time]],
                                   columns=['SampleID','mean of median fluorescence (a.u.)',
                                            'std of median fluorescence (a.u.)',
                                            'Induction','Post-induction (hrs)'])
            stat_data = stat_data.append(rowInfo,sort=False)

#%% Export statiscs to file
outputDataFP = os.path.join(dataRootDir,dataFolderDir,statFN)
stat_data.to_csv(outputDataFP)

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