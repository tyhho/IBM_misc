# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 18:30:52 2019

@author: Trevor Ho
"""

import os
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import gridspec

import seaborn as sns
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42


# Specify root dir
dataRootDir = r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'

# TODO: Specify CSV file containing the fluorescence data to plot
# TODO: Specify CSV file containing the information for plotting horizontal ref lines
dataRootDir=r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir='FC021'
dataFN = 'IBM_FC021R2-7_Finalized_FCmedian&metadata&PRData.csv'
hlineDataFN = 'IBM_FC021R2-7_FCmedian&metadata&PRData_CtrlStats.csv'

# TODO: Specify filename to be saved. Figure will be saved under path of dataRootDir/dataFolderDir
figName = 'IBM_FC021R2-7_Finalized_Result.pdf'

# Read data
csvDir = os.path.join(dataRootDir,dataFolderDir,dataFN)
data = pd.read_csv(csvDir,index_col=0)

# Read hline data
hlineDir = os.path.join(dataRootDir,dataFolderDir,hlineDataFN)
hlineData = pd.read_csv(hlineDir,index_col=0)

pi_dict = {0:5, 1:24}
sample_dict = {
        0:'IBMs0164',
        1:'IBMs0161',
        2:'IBMs0165'
        }

#%%
# Create canvas
sns.set(style='ticks')
fig, ax = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(8, 4.5), dpi = 200)
# Grid spec

gs = gridspec.GridSpec(2, 3) 

palette1 = sns.color_palette("hls", 8)

for ax_row_ID, pi_time in pi_dict.items():
    for ax_col_ID, sampleID in sample_dict.items():
        
        ax[ax_row_ID,ax_col_ID] = plt.subplot(gs[ax_row_ID,ax_col_ID])

        data2plot = data[(data['SampleID']==sampleID) & (data['Post-induction (hrs)']==pi_time)]
        sns.lineplot(x = 'caffeine (μM)', y = 'median fluorescence (a.u.)', data=data2plot,
                          hue='arabinose (mM)', palette=palette1,
                          err_style = 'bars', zorder = 9,
                          ax=ax[ax_row_ID,ax_col_ID])
        
        # Remove legend from plot
        ax[ax_row_ID,ax_col_ID].get_legend().remove()
        
        
        ax[ax_row_ID,ax_col_ID].set_xscale('symlog', linthreshx=0.1, linscalex=0.35
          , subsx=[2, 3, 4, 5, 6, 7, 8, 9]
          )
        ax[ax_row_ID,ax_col_ID].set(yscale="log")
        ax[ax_row_ID,ax_col_ID].set_xlim(-0.05, 2e2)
        
        if ax_row_ID == 0:
            ax[ax_row_ID,ax_col_ID].set_ylim(1e2, 5e3)
        elif ax_row_ID == 1:
            ax[ax_row_ID,ax_col_ID].set_ylim(1e2, 5e4)
        ax[ax_row_ID,ax_col_ID].spines['bottom'].set_bounds(0.06, 2e2) # values likely need to be adjusted in a new plot
        ax[ax_row_ID,ax_col_ID].set_xlabel('')
        ax[ax_row_ID,ax_col_ID].set_ylabel('')
        
        # Draw the left part of the broken x axis by using twiny()
        # Source: https://matplotlib.org/gallery/ticks_and_spines/multiple_yaxis_with_spines.html
        def make_patch_spines_invisible(ax):
            ax.set_frame_on(True)
            ax.patch.set_visible(False)
            for sp in ax.spines.values():
                sp.set_visible(False)
        
        lxAxis = ax[ax_row_ID,ax_col_ID].twiny()
        make_patch_spines_invisible(lxAxis)
        lxAxis.spines["bottom"].set_visible(True)
        lxAxis.spines["bottom"].set_bounds(-0.05, 0.08) # values likely need to be adjusted in a new plot
        lxAxis.tick_params(axis='both', which='both',top=False,bottom=False,labeltop=False)

        # Draw the two slashes
        # Source: https://matplotlib.org/examples/pylab_examples/broken_axis.html
        d = .02  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        kwargs = dict(transform=ax[ax_row_ID,ax_col_ID].transAxes, color='k', clip_on=False)
        ax[ax_row_ID,ax_col_ID].plot((0.07,0.09), (-d,+d), **kwargs, zorder=20)
        ax[ax_row_ID,ax_col_ID].plot((0.10,0.12), (-d,+d), **kwargs, zorder=21)
        
        # Plot horizontal lines
        hline_data_to_plot = hlineData[hlineData['Post-induction (hrs)'] == pi_time]
        
        for hline_data_to_plot_row in hline_data_to_plot.itertuples():
            hline_mean = hline_data_to_plot_row[2]
            hline_sd = hline_data_to_plot_row[3]
            hline_color = hline_data_to_plot_row[5]
            hline_shade = hline_data_to_plot_row[6]

            if hline_data_to_plot_row[4] == True:
                ax[ax_row_ID,ax_col_ID].axhspan(hline_mean - hline_sd, hline_mean + hline_sd, facecolor=hline_shade, zorder=0)
            else:
                ax[ax_row_ID,ax_col_ID].axhline(hline_mean, color=hline_shade, zorder=0)
            
            ax[ax_row_ID,ax_col_ID].axhline(hline_mean, color=hline_color, dashes=[2,2], linestyle = '--', zorder=0)    

# Remove xtick labels and ytick labels for non outer-frames
for ax_row_ID in [0,1]:
    for ax_col_ID in [0,1,2]:
        if ax_row_ID != 1:
            ax[ax_row_ID,ax_col_ID].set_xticklabels('')
#        else:
#            ax[ax_row_ID,ax_col_ID].tick_params(axis='x', which='major', labelsize=9)

        if ax_col_ID != 0:
            ax[ax_row_ID,ax_col_ID].set_yticklabels('')
        

ax[1,1].set_xlabel('caffeine (μM)')
ax[0,0].set_ylabel('5 hours')
ax[1,0].set_ylabel('24 hours')

# TODO: Customize legend for fluorescence
# frameon=False removes the frame where the labels were in
# Get the hue category label and removes it    
handles, labels = ax[0,0].get_legend_handles_labels()
ax[1,2].legend(
#        handles=handles[1:], labels=labels[1:],
      bbox_to_anchor=(1.005, 1),frameon=True,
      prop={'size': 8})

fig.tight_layout()

'''Save Figure'''
outputDir = os.path.join(dataRootDir,dataFolderDir,figName)
fig.savefig(outputDir)