# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 19:36:40 2019

@author: Trevor Ho
"""
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
# matplotlib.rcParams['ps.fonttype'] = 42
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
dataFolderDir='FC033'
statDataFN = 'IBM_FC033R1,4,5_median&metadata&PRData_InductionRelabelled_Stats.csv'

hlineInfoFN = 'IBM_FC033R1,4,5_median&metadata&PRData_hlineInfo.csv'
figName = 'IBM_FC033R1,4,5_Result.pdf'

# Read data
statDataFP = os.path.join(dataRootDir,dataFolderDir,statDataFN)
statData = pd.read_csv(statDataFP,index_col=0)
hlineInfoFP = os.path.join(dataRootDir,dataFolderDir,hlineInfoFN)
hlineInfo = pd.read_csv(hlineInfoFP,index_col=0)

# TODO: Set induction and induction time information
inductions = ['no induction','1 mM arabinose','25 μM DAPG','1 mM arabinose + 25 μM DAPG']
pi_times = [5,24]

#%%
## Order the data from highest to lowest fluorescence based on induced fluorescence of PI24

# TODO: Update control list
control_list = ['IBMc186 + IBMc101',
                'IBMc120 + IBMc101',
                'IBMc120 + IBMc071'
                ]
# TODO: Set x tick labels of controls
ctrl_tick_labels = ['background','reporter','ECF20']

# Take out control data for ordering and focus on PI24 when induced
# TODO: Inspect data and set up criteria for ordering
ordered_df = statData[(~statData['SampleID'].isin(control_list)) & (statData['Post-induction (hrs)']==24) & \
               (statData['Induction']=='1 mM arabinose + 25 μM DAPG')]

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
statData = pd.merge(statData,ordered_df)

# TODO: Set value for plots
pi_dict = {0:5, 1:24}

#%% Plot

# Remove unnecessary data
candidate_stat_data = statData[(~statData['SampleID'].isin(control_list))]
ctrl_stat_data = statData[statData['SampleID'].isin(control_list)]

# Create Canvas
sns.set(style='ticks')
fig, ax = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(8, 6), dpi = 200)

# TODO: Decide color for inductions
color_mapper = {
        # '+ 10 μM 4-HT':'#F96495',
        # '+ DMSO':'#1D1D1B' ,
        'no induction':'#1D1D1B',
        '1 mM arabinose':'#55A0FB',
        '25 μM DAPG':'#099963',
        '1 mM arabinose + 25 μM DAPG':'#F96495' ,
        }

# Calculate width ratio by number of elements
ctrl_len = len(control_list)
candidate_len = len(candidate_stat_data)/len(inductions)/len(pi_dict)

# Grid spec
gs = gridspec.GridSpec(2, 2, width_ratios=[ctrl_len+2, candidate_len+4]) 

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
        
        # TODO: set y axis limits
        if row_ax_ID == 0:
            # ax[row_ax_ID,col_ax_ID].set_ylim(1e2, 1.5e3)
            ax[row_ax_ID,col_ax_ID].set_ylim(1e2, 1e5)
        elif row_ax_ID == 1:
            # ax[row_ax_ID,col_ax_ID].set_ylim(1e2, 1e4)
            ax[row_ax_ID,col_ax_ID].set_ylim(1e2, 1e5)
       
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
        hlineInfo_to_plot = hlineInfo[hlineInfo['Post-induction (hrs)'] == pi_time]
        
        for hlineInfo_to_plot_row in hlineInfo_to_plot.itertuples():
            hline_mean = hlineInfo_to_plot_row[2]
            hline_sd = hlineInfo_to_plot_row[3]
            hline_color = hlineInfo_to_plot_row[5]
            hline_shade = hlineInfo_to_plot_row[6]
            
            if hlineInfo_to_plot_row[4] == True:
                ax[row_ax_ID,1].axhspan(hline_mean - hline_sd, hline_mean + hline_sd, facecolor=hline_shade, zorder=0)
                ax[row_ax_ID,0].axhspan(hline_mean - hline_sd, hline_mean + hline_sd, facecolor=hline_shade, zorder=0)
            else:
                ax[row_ax_ID,1].axhline(hline_mean, color=hline_shade, zorder=0)
                ax[row_ax_ID,0].axhline(hline_mean, color=hline_shade, zorder=0)
            ax[row_ax_ID,col_ax_ID].axhline(hline_mean, color=hline_color, dashes=[2,2], linestyle = '--', zorder=1)

# TODO: Customize legend for fluorescence
# frameon=False removes the frame where the labels were in
# Get the hue category label and removes it    
handles, labels = ax[0,1].get_legend_handles_labels()
labels = [
        '',
        # '+ DMSO',
        # '+ 10 μM 4-HT'
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
# ax[1,0].text(-3, 1500, 'median fluorescence (a.u.)', fontsize=12, rotation=90)
ax[1,0].text(-7, 1500, 'median fluorescence (a.u.)', fontsize=12, rotation=90)

# TODO: adjust plot
# plt.gcf().subplots_adjust(bottom=0.3, hspace=0.1, wspace=0.05)
plt.gcf().subplots_adjust(bottom=0.3, hspace=0.1, wspace=0.05)

'''Save Figure'''
outputDir = os.path.join(dataRootDir,dataFolderDir,figName)
fig.savefig(outputDir)