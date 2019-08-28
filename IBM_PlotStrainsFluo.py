# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 19:36:40 2019

@author: Trevor Ho
"""
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rc('ps', fonttype=42)
from matplotlib import gridspec


import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

sns.reset_defaults

# TODO: Specify folder location
    # csv file containing data to be plotted should be in this folder

dataRootDir=r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir='FC023'
dataFN = 'IBM_FC023R3-5_FCmedian&metadata&PRData_InductionProcessed_FoldChangeFiltered.csv'
statFN = 'IBM_FC023R3-5_FCmedianStats&metadata_FoldChangeFiltered.csv'
figName = 'IBM_FC013R3-5_Result_FoldChangeFiltered.pdf'

# TODO: Set induction information
inductions = ['no induction','1 mM arabinose','25 μM DAPG','1 mM arabinose + 25 μM DAPG']
pi_times = [5,24]

# Read data
csvDir = os.path.join(dataRootDir,dataFolderDir,dataFN)
data = pd.read_csv(csvDir,index_col=0)

## Generate data for mean and sd, according to the induction
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
outputDir = os.path.join(dataRootDir,dataFolderDir,statFN)
stat_data.to_csv(outputDir)

#%%
## Order the data from highest to lowest fluorescence based on induced fluorescence of PI24

# TODO: Update control list
control_list = ['IBMc252 + IBMc249',
        'IBMc259 + IBMc249',
        'IBMc259 + IBMc226']
# TODO: Set x tick labels of controls
ctrl_tick_labels = ['background','reporter','TetR + reporter']

# Take out control data for ordering and focus on PI24 when induced
# TODO: Inspect data and set up criteria for ordering
ordered_df = stat_data[(~stat_data['SampleID'].isin(control_list)) & (stat_data['Post-induction (hrs)']==24) & \
               (stat_data['Induction']=='1 mM arabinose + 25 μM DAPG')]

# Keep only the median fluorescence column for sorting
ordered_df = ordered_df[['SampleID','mean of median fluorescence (a.u.)']]
# Sort by fluorescence and assign a sorting index to each sampleID
ordered_df = ordered_df.sort_values(by=['mean of median fluorescence (a.u.)'],ascending=False)
ordered_df['SortIndex'] = list(range(1,len(ordered_df)+1))

# Return to the data to assign sorting index to sampleID
controls_order_df = pd.DataFrame(control_list,columns=['SampleID'])
controls_order_df['SortIndex'] = list(range(1,len(control_list)+1))

ordered_df = ordered_df.drop('mean of median fluorescence (a.u.)',axis=1)
ordered_df = ordered_df.append(controls_order_df,ignore_index=True,sort=False)
data = pd.merge(data,ordered_df)
stat_data = pd.merge(stat_data,ordered_df)

#%%

# TODO: Input information for constructing hlines

# hline_input instruction:
    # key = str for df.query. Note that str inside str are double-quoted
    # value = list of colors to use, 1st color = line color, 2nd color = shade color
hline_input = {
        'SampleID == "IBMc252 + IBMc249"': ['k','#E8E8E8'],
        'SampleID == "IBMc259 + IBMc249"': ['#F96495','#eda3ba'],
        'SampleID == "IBMc259 + IBMc226" & (Induction == "no induction" | Induction == "25 μM DAPG")': ['#A952A5','#E9D4E9'],
        'SampleID == "IBMc259 + IBMc226" & (Induction == "1 mM arabinose" | Induction == "1 mM arabinose + 25 μM DAPG")': ['#A952A5','#E9D4E9'],
        }

#%% Calculate stats for hlines
pi_dict = {0:5, 1:24}
hline_data = pd.DataFrame(columns=[])

for key, pi_time in pi_dict.items():
    for queryLine, color in hline_input.items():
        hRefSubData = data.query(queryLine)
        hRefSubData = hRefSubData[hRefSubData['Post-induction (hrs)']==pi_time]
        hRefSubData = hRefSubData['median fluorescence (a.u.)']        
        hline_data_row_dict = {
                'Post-induction (hrs)': [pi_time],
                'mean of median fluorescence (a.u.)': [hRefSubData.mean()],
                'std of median fluorescence (a.u.)': [hRefSubData.std()],
                'span for std': [(hRefSubData.std() > hRefSubData.mean()*0.1)],
                'line color': [color[0]],
                'shade color': [color[1]]
                }
        hline_data = hline_data.append(pd.DataFrame.from_dict(hline_data_row_dict))

#%% Plot

# Remove unnecessary data
candidate_stat_data = stat_data[(~stat_data['SampleID'].isin(control_list))]
ctrl_stat_data = stat_data[stat_data['SampleID'].isin(control_list)]

# Create Canvas
sns.set(style='ticks')
fig, ax = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(8, 6), dpi = 200)

# TODO: Decide color for inductions
color_mapper = {
        'no induction':'#F96495',
        '1 mM arabinose':'#55A0FB',
        '25 μM DAPG':'#099963',
        '1 mM arabinose + 25 μM DAPG':'#1D1D1B'
        }

# Calculate width ratio by number of elements
ctrl_len = len(control_list)
candidate_len = len(candidate_stat_data)/len(inductions)/len(pi_dict)

# Grid spec
gs = gridspec.GridSpec(2, 2, width_ratios=[ctrl_len+10, candidate_len+4]) 

# Plot graph of candidate strains
for col_ax_ID in [0,1]:
    for row_ax_ID, pi_time in pi_dict.items():
        
        # Create df to isolate data of interest
        if col_ax_ID == 0:
            stat_data_to_plot = ctrl_stat_data[ctrl_stat_data['Post-induction (hrs)']==pi_time]
        elif col_ax_ID == 1:
            stat_data_to_plot = candidate_stat_data[candidate_stat_data['Post-induction (hrs)']==pi_time]
        
        ax[row_ax_ID,col_ax_ID] = plt.subplot(gs[row_ax_ID,col_ax_ID])
        # Plot data
        sns.scatterplot(x='SortIndex', y='mean of median fluorescence (a.u.)',
                        hue='Induction',
                       hue_order=inductions,
                       palette=color_mapper,
                       data=stat_data_to_plot,
                       ax=ax[row_ax_ID,col_ax_ID],zorder=10)
        
        ax[row_ax_ID,col_ax_ID].set(yscale="log")
        ax[row_ax_ID,col_ax_ID].set_ylim(1e2, 3e5)
       
        ax[row_ax_ID,col_ax_ID].set_xlabel('')
        ax[row_ax_ID,col_ax_ID].set_xticklabels('')
    
        # Plot error bars
        for induction in inductions:
            stat_data_to_plot_ind = stat_data_to_plot[stat_data_to_plot['Induction']==induction]
            sortIDVec = stat_data_to_plot_ind['SortIndex'].tolist()
            meanVec = stat_data_to_plot_ind['mean of median fluorescence (a.u.)'].tolist()
            errVec = stat_data_to_plot_ind['std of median fluorescence (a.u.)'].tolist()
            ax[row_ax_ID,col_ax_ID].errorbar(x=sortIDVec, y=meanVec, yerr = errVec, fmt='none',
              zorder=6, color=color_mapper[induction]
              )
    
        # Remove legend from plot
        ax[row_ax_ID,col_ax_ID].get_legend().remove()
        
        # Plot horizontal lines
        hline_data_to_plot = hline_data[hline_data['Post-induction (hrs)'] == pi_time]
        
        for hline_data_to_plot_row in hline_data_to_plot.itertuples():
            hline_mean = hline_data_to_plot_row[2]
            hline_sd = hline_data_to_plot_row[3]
            hline_color = hline_data_to_plot_row[5]
            hline_shade = hline_data_to_plot_row[6]
            
    #        if hline_data_to_plot_row[4] == True:
    #            ax[ax_ID,1].axhspan(hline_mean - hline_sd, hline_mean + hline_sd, facecolor=hline_shade, zorder=0)
    #        else:
    #            ax[ax_ID,1].axhline(hline_mean, color=hline_shade, zorder=0)
            ax[row_ax_ID,col_ax_ID].axhline(hline_mean, color=hline_color, dashes=[2,2], linestyle = '--', zorder=0)

# TODO: Customize legend for fluorescence
# frameon=False removes the frame where the labels were in
# Get the hue category label and removes it    
handles, labels = ax[0,1].get_legend_handles_labels()
labels = [
        '',
        'no induction',
        '1 mM arabinose (N-lobe)',
        '25 μM DAPG (C-lobe)',
        '1 mM arabinose + 25 μM DAPG (N+C lobes)'
        ]
ax[1,1].legend(
        handles=handles[1:], labels=labels[1:],
      bbox_to_anchor=(1, -0.05),frameon=False, ncol=2,
      prop={'size': 9})

# Control group labels and ticks
for i in [0,1]:
    ax[i,0].set_xlim(0.5, ctrl_len+0.5)
    ax[i,0].set_xticks(range(1,ctrl_len+1))

for i in [0,1]:
    ax[i,1].set_yticklabels('')
    ax[i,1].set_ylabel('')
    ax[i,1].set_xlim(0, candidate_len+1)
    ax[i,1].get_xaxis().set_ticks([])

# Customize outer labels
ax[0,0].set_ylabel('5 hours')
ax[1,0].set_ylabel('24 hours')
ax[1,0].set_xticklabels(ctrl_tick_labels, rotation=40, fontsize=10, ha="right")

# TODO: adjust label of y axis
ax[1,0].text(-5, 1500, 'median fluorescence (a.u.)', fontsize=12, rotation=90)

# TODO: adjust plot
plt.gcf().subplots_adjust(bottom=0.3, hspace=0.1, wspace=0.05)

'''Save Figure'''
outputDir = os.path.join(dataRootDir,dataFolderDir,figName)
fig.savefig(outputDir)