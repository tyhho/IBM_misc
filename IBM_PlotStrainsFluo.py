# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 19:36:40 2019

@author: Trevor Ho
"""
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rc('ps', fonttype=42)

import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

sns.reset_defaults

# Specify folder location
    # csv file containing data to be plotted should be in this folder

dataRootDir=r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir='FC023'
datafilename = 'IBM_FC023R3-5_FCmedian&metadata&PRData_InductionProcessed.csv'
figName = 'IBM_FC013R3-5_Result.eps'

csvDir = os.path.join(dataRootDir,dataFolderDir,datafilename)
data = pd.read_csv(csvDir,index_col=0)

## Generate data for mean and sd, according to the induction
inductions = ['no induction','1 mM arabinose','25 μM DAPG','1 mM arabinose + 25 μM DAPG']
pi_times = [5,24]

stat_data = pd.DataFrame(columns=[])
for pi_time in pi_times:
    for induction in inductions:
        subdata = data[(data['Induction']==induction) & (data['Post-induction (hrs)']==pi_time)]
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
outputFilename = 'IBM_FC023R3-5_FCmedianStats&metadata.csv'
outputDir = os.path.join(dataRootDir,dataFolderDir,outputFilename)
stat_data.to_csv(outputDir)

#%%
## Order the data from highest to lowest fluorescence based on induced fluorescence of PI24

# Preparations
control_list = ['IBMc252 + IBMc249',
        'IBMc259 + IBMc249',
        'IBMc259 + IBMc226']

# Take out control data for ordering and focus on PI24 when induced
order_df = stat_data[(~stat_data['SampleID'].isin(control_list)) & (stat_data['Post-induction (hrs)']==24) & \
               (stat_data['Induction']=='no induction')]

# Keep only the median fluorescence column for sorting
order_df = order_df[['SampleID','mean of median fluorescence (a.u.)']]
# Sort by fluorescence and assign a sorting index to each sampleID
order_df = order_df.sort_values(by=['mean of median fluorescence (a.u.)'],ascending=False)
order_df['SortIndex'] = list(range(len(control_list)+1,len(control_list)+len(order_df)+1))

# Save order of sample IDs as labels
x_tick_labels_samples = order_df['SampleID'].tolist()
x_tick_labels = ['background','reporter only','TetR + reporter']
#x_tick_labels = x_tick_labels + x_tick_labels_samples

# Return to the data to assign sorting index to sampleID
controls_order_df = pd.DataFrame(control_list,columns=['SampleID'])
controls_order_df['SortIndex'] = list(range(1,len(control_list)+1))

order_df = order_df.drop('mean of median fluorescence (a.u.)',axis=1)
order_df = order_df.append(controls_order_df,ignore_index=True,sort=False)
data = pd.merge(data,order_df)
stat_data = pd.merge(stat_data,order_df)

stat_data = stat_data.sort_values(by=['SortIndex'])


#%% Plot

# Remove unnecessary data
data2 = data[['median fluorescence (a.u.)','Induction','SortIndex','Post-induction (hrs)']]
stat_data2 =stat_data[['mean of median fluorescence (a.u.)',
                              'std of median fluorescence (a.u.)',
                              'Induction','SortIndex','Post-induction (hrs)']]

# Create Canvas
sns.set(style='ticks')
fig, ax = plt.subplots(2, 1, sharex=True, figsize=(8, 6), dpi = 200)

pi_dict = {0:5, 1:24}

for ax_ID, pi_time in pi_dict.items():
    
    # Create df to isolate data of interest
    subdata = data2[data2['Post-induction (hrs)']==pi_time]
    sub_stat_data = stat_data2[stat_data2['Post-induction (hrs)']==pi_time]
    
    # Plot data
    sns.scatterplot(x='SortIndex', y='median fluorescence (a.u.)', 
#                   hue_order=['no induction','1 mM arabinose'],
                   ci='sd',
                   hue ='Induction',
#                   palette = ['#808080','b','g','r'],
                   data=subdata,
                   ax=ax[ax_ID],zorder=10)
    
    # Calculate control data
    negdata = subdata[subdata['SortIndex']==1]
    negFluos = negdata['median fluorescence (a.u.)']
    neg_mean = negFluos.mean()
    neg_sd = negFluos.std()
    
    posdata = subdata[subdata['SortIndex']==2]
    posFluos = posdata['median fluorescence (a.u.)']
    pos_mean = posFluos.mean()
    pos_sd = posFluos.std()
    
    neg_sd_span = True
    if neg_sd < neg_mean*0.05:
        neg_sd_span = False
    
    pos_sd_span = True
    if pos_sd < pos_mean*0.05:
        
        pos_sd_span = False
    
    # Plot negative control
    if neg_sd_span == True:
        ax[ax_ID].axhspan(neg_mean - neg_sd, neg_mean + neg_sd, facecolor='#E8E8E8', zorder=0)
    else:
        ax[ax_ID].axhline(neg_mean, color='#E8E8E8', zorder=0)
    ax[ax_ID].axhline(neg_mean, color='k', dashes=[2,2], linestyle = '--', zorder=0)
    
    # Plot positive control
    if pos_sd_span == True:
        ax[ax_ID].axhspan(pos_mean - pos_sd, pos_mean + pos_sd, facecolor='#eda3ba', zorder=0)
    else:
        ax[ax_ID].axhline(pos_mean, color='#eda3ba', zorder=0)
    ax[ax_ID].axhline(pos_mean, color='r', dashes=[2,2], linestyle = '--', zorder=0)
    
    # Remove legend from plot
    ax[ax_ID].get_legend().remove()


#    for induction in inductions: # in case scatter plot needs to be used
#        sub_stat_data2 = sub_stat_data[sub_stat_data['Induction']=='1 mM arabinose']
#        sortIDVec = sub_stat_data2['SortIndex'].tolist()
#        meanVec = sub_stat_data2['mean of median fluorescence (a.u.)'].tolist()
#        errVec = sub_stat_data2['std of median fluorescence (a.u.)'].tolist()
    
    ax[ax_ID].set(yscale="log")
    ax[ax_ID].set_ylim(1e2, 1e5)
    ax[ax_ID].set_xticklabels(x_tick_labels)
    ax[ax_ID].set_xlabel('')
    plt.xticks(rotation='vertical')
    
    # Distinguish between controls and strains
    ax[ax_ID].axvline(4.5, color='k', dashes=[2,2], linestyle = '--', zorder=6)

# Customize y labels of the two different graphs
ax[0].set_ylabel('5 hours')
ax[1].set_ylabel('24 hours')

ax[0].text(-7, 1000, 'median fluorescence (a.u.)', fontsize=12, rotation=90)

# Customize legend for fluorescence
# frameon=False removes the frame where the labels were in
# Get the hue category label and removes it    
#handles, labels = ax[1].get_legend_handles_labels()
ax[0].legend(
#        handles=handles[1:], labels=labels[1:],
      loc='upper right',frameon=False)
plt.gcf().subplots_adjust(bottom=0.3)

'''Save Figure'''
outputDir = os.path.join(dataRootDir,dataFolderDir,figName)
#fig.savefig(outputDir)