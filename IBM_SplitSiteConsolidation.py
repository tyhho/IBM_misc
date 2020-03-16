# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 17:24:54 2019

@author: Trevor Ho

This script consolidate the split sites and calculate mean of median fluorescences.
Then plot everything on a graph
"""

import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter, LogLocator, MultipleLocator
from matplotlib import gridspec
import numpy as np

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42

sns.reset_defaults

# Specify root dir
dataRootDir = r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'

# TODO: Specify location of CSV file containing identified split sites
ssFolderDir1 = 'BM010'
ssFolderDir2 = 'Sequencing Results'
ssDataFN = 'IBM_BM010_IdentifiedSplitSites.csv'

# TODO: Specify location of CSV file containing the median fluorescence values of the strains
fluoFolderDir='FC033'
fluoDataFN = 'IBM_FC033R1,4,5_median&metadata&PRData_InductionRelabelled.csv'

# TODO: Specify location of CSV file containing the median fluorescence values of control strains
ctrlFluoDataFN = 'IBM_FC033R1,4,5_median&metadata&PRData_InductionRelabelled_Stats.csv'

# TODO: Specify filename to be saved. Figure will be saved under path of dataRootDir\ssFolderDir1
figName = 'IBM_BM010_InsertionMap_Raw.pdf'

# TODO: Specify filename of csv with merged data to be saved (final data saves under BM005 directly)
mergedDataFN = 'IBM_BM10_PooledResults.csv'

hlineInfoFN = 'IBM_FC033R1,4,5_median&metadata&PRData_hlineInfo.csv'

# Preparations

# TODO: Set induction and induction time information
# inductions = ['+ DMSO','+ 10 μM 4-HT']
inductions = ['no induction','1 mM arabinose','25 μM DAPG','1 mM arabinose + 25 μM DAPG']
pi_times = [5,24]

# TODO: Update control list
# control_list = ['IBMc101',
#         'IBMc307']
control_list = ['IBMc186 + IBMc101',
                'IBMc120 + IBMc101',
                'IBMc120 + IBMc071'
                ]

# TODO: Set x tick labels of controls
ctrl_tick_labels = ['background','M86']

# TODO: Define how many amino acids are being plotted
start_aa = 1
end_aa = 193

# TODO: Input transposition window
# n_trans_border = 7
# c_trans_border = 148
n_trans_border = 6
c_trans_border = 188

# TODO: Decide color for inductions
# color_mapper = {
#         '+ 10 μM 4-HT':'#F96495',
#         '+ DMSO':'#1D1D1B' 
#         }
color_mapper = {
        'no induction':'#1D1D1B',
        '1 mM arabinose':'#55A0FB',
        '25 μM DAPG':'#099963',
        '1 mM arabinose + 25 μM DAPG':'#F96495' ,
        }


#%%
'''Read data and map split sites to function'''
# Read split site csv file
ssDataDir = os.path.join(dataRootDir,ssFolderDir1,ssFolderDir2,ssDataFN)
ssData = pd.read_csv(ssDataDir,index_col=0)
ssData = ssData[['aa_ss_middle']]
ssData['SampleID']=ssData.index

# Read fluorescence csv file
fluoDataDir = os.path.join(dataRootDir,fluoFolderDir,fluoDataFN)
fluoData = pd.read_csv(fluoDataDir,index_col=0)

# Map split site to fluorescence
mergedData = fluoData.merge(ssData,sort=False)
# Generate a stripped mergedData set for histogram plot
mergedDataSingleCon=mergedData[(mergedData['Run']==4) \
                               & (mergedData['Induction']=='+ DMSO') \
                               & (mergedData['Post-induction (hrs)']==5)]
mergedDataSingleCon=mergedDataSingleCon['aa_ss_middle']

'''Loop through the amino acid split sites and consolidate median fluorescence'''

# Deduplicate the list of amino acid split sites and sort to a list
unique_ss_list = ssData['aa_ss_middle'].drop_duplicates().sort_values().tolist()

unique_ss_data = pd.DataFrame(columns=[])
for pi_time in pi_times:
    for induction in inductions:
        for ss in unique_ss_list:
            subdata = mergedData[(mergedData['Induction']==induction) \
                     & (mergedData['Post-induction (hrs)']==pi_time) \
                     & (mergedData['aa_ss_middle']==ss) 
                     ]
            medianFluos = subdata['median fluorescence (a.u.)']
            mean_medianFluo = medianFluos.mean()
            sd_medianFluo = medianFluos.std()
            count = len(subdata)/3  # divided by 3 because each subdata contains data from 3 runs
            rowInfo = pd.DataFrame([[ss,mean_medianFluo,sd_medianFluo,induction,pi_time,count]],
                                   columns=['position (aa)','mean of median fluorescence (a.u.)',
                                            'std of median fluorescence (a.u.)',
                                            'Induction','Post-induction (hrs)','count'])
            unique_ss_data = unique_ss_data.append(rowInfo,sort=False)

# Export the data to keep it somewhere safe
mergedDataDir = os.path.join(dataRootDir,ssFolderDir1,mergedDataFN)
mergedData.to_csv(mergedDataDir)
            
#%% Calculate stats for hlines
pi_dict = {1:5, 2:24}

hlineInfoFP = os.path.join(dataRootDir,fluoFolderDir,hlineInfoFN)
hlineInfo = pd.read_csv(hlineInfoFP,index_col=0)


# Import CSV containing control data for plotting
ctrlFluoDataDir = os.path.join(dataRootDir,fluoFolderDir,ctrlFluoDataFN)
ctrlFluoData = pd.read_csv(ctrlFluoDataDir,index_col=0)
ctrlFluoData = ctrlFluoData[ctrlFluoData['SampleID'].isin(control_list)]
ctrlFluoData['SortIndex'] = ctrlFluoData.apply(lambda data: control_list.index(data['SampleID'])+1, axis=1)
#%%
# Create canvas
sns.set(style='ticks')
fig, ax = plt.subplots(5, 2, sharex=True, sharey=True, figsize=(4, 6), dpi = 200)
# Grid spec
ctrl_len = len(control_list)
gs = gridspec.GridSpec(5, 2, width_ratios=[ctrl_len+10, end_aa - start_aa], height_ratios = [0.6,1,1,0.2,0.2]) 

# Delete top-left plot
fig.delaxes(ax[0,0])
fig.delaxes(ax[3,0])

# Plot the counts of strains being mapped to the same split site
ax[0,1] = plt.subplot(gs[0,1])
sns.distplot(mergedDataSingleCon, kde=False, color="k", ax=ax[0,1], bins=end_aa, hist_kws=dict(edgecolor='none'))
ax[0,1].set_xlabel('')
ax[0,1].set_ylabel('sequencing \n count')
ax[0,1].tick_params(axis='y', which='minor', left=True)
minorLocator = AutoMinorLocator(5)
ax[0,1].yaxis.set_minor_locator(minorLocator)
ax[0,1].axvline(n_trans_border, color='k', dashes=[2,2], linestyle = '--', zorder=5)
ax[0,1].axvline(c_trans_border, color='k', dashes=[2,2], linestyle = '--', zorder=5)
ax[0,1].set_xlim(start_aa,end_aa)
ax[0,1].set_xticklabels('')

# Plot the graphs

for col_ax_ID in [0,1]:
    for row_ax_ID, pi_time in pi_dict.items():
        
        if col_ax_ID ==0:
            # Create df to isolate data of interest
            subdata = ctrlFluoData[ctrlFluoData['Post-induction (hrs)']==pi_time]
            xaxis = 'SortIndex'
        elif col_ax_ID == 1:
            # Create df to isolate data of split sites for each time point
            subdata = unique_ss_data[unique_ss_data['Post-induction (hrs)']==pi_time]
            xaxis = 'position (aa)'
        
        ax[row_ax_ID,col_ax_ID] = plt.subplot(gs[row_ax_ID,col_ax_ID])
        sns.scatterplot(x=xaxis, y='mean of median fluorescence (a.u.)', 
                           hue='Induction',
                           hue_order=inductions,
                           palette=color_mapper,
                           data=subdata, ax=ax[row_ax_ID,col_ax_ID],zorder=10,
                           size = 3,
                           edgecolor=None)
            
        # Plot error bars
        for induction in inductions:
            subdata_ind = subdata[subdata['Induction']==induction]
            xVec = subdata_ind[xaxis].tolist()
            meanVec = subdata_ind['mean of median fluorescence (a.u.)'].tolist()
            errVec = subdata_ind['std of median fluorescence (a.u.)'].tolist()
            ax[row_ax_ID,1].errorbar(x=xVec, y=meanVec, yerr = errVec, fmt='none',
              zorder=6, color=color_mapper[induction]
              )

        # Remove legend from plot
        ax[row_ax_ID,col_ax_ID].get_legend().remove() 
                
        # TODO: Check range of fluorescence
        
        # ax[row_ax_ID,col_ax_ID].set_ylim(1e2, 3e5)
        
        ax[row_ax_ID,col_ax_ID].set(yscale="log")
        
        ax[row_ax_ID,col_ax_ID].yaxis.set_ticks(np.logspace(2, 5, num=4))
        locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
        ax[row_ax_ID,col_ax_ID].yaxis.set_minor_locator(locmin)       

        
        # Plot horizontal lines
        hline_data_to_plot = hlineInfo[hlineInfo['Post-induction (hrs)'] == pi_time]
        
        for hline_data_to_plot_row in hline_data_to_plot.itertuples():
            hline_mean = hline_data_to_plot_row[2]
            hline_sd = hline_data_to_plot_row[3]
            hline_color = hline_data_to_plot_row[5]
            hline_shade = hline_data_to_plot_row[6]
            
#            if hline_data_to_plot_row[4] == True:
#                ax[ax_ID,1].axhspan(hline_mean - hline_sd, hline_mean + hline_sd, facecolor=hline_shade, zorder=0)
#            else:
#            ax[row_ax_ID,1].axhline(hline_mean, color=hline_shade, zorder=0)
            ax[row_ax_ID,col_ax_ID].axhline(hline_mean, color=hline_color, dashes=[2,2], linestyle = '--', zorder=0)    
        
    #    # Annotate split sites to match 3D structures
    #    ss_color_list = ['#f9897a',
    #                     '#f9a4ed',
    #                     '#f9a4ed',
    #                     '#ffcd70',
    #                     '#7ffaff',
    #                     '#f2f90e',
    #                     '#91a9ea',
    #                     '#025ff4',
    #                     '#025ff4',
    #                     '#77fc5f',
    #                     '#77fc5f',
    #                     '#77fc5f'               
    #            ]
    #    ss_color_map = dict(zip(unique_ss_list,list(zip(meanVec,ss_color_list))))
    #    for ss, color_info in ss_color_map.items():
    ##    ax[3].axvline(ss, color=color, linestyle = '-',zorder=0)
    #        ax[ax_ID].annotate('',xy=(ss+0.25, color_info[0]*1.2), xytext=(ss+1.5, color_info[0]*2),
    #              zorder=3,
    #              arrowprops=dict(facecolor=color_info[1], width=5, headwidth=8, **{'edgecolor':'None'}))
    

# Extra customizations to control fluorescence plots
for row_ax_ID in [1,2]:
    ax[row_ax_ID,0].set_xlim(0, ctrl_len+1)
    ax[row_ax_ID,0].set_xticks(range(1,ctrl_len+1))
    ax[row_ax_ID,0].set_xticklabels('')
    ax[row_ax_ID,0].set_xlabel('')

# ax[1,0].set_ylim(1e2, 1.5e3)
# ax[1,1].set_ylim(1e2, 1.5e3)
# ax[2,0].set_ylim(1e2, 1e4)
# ax[2,1].set_ylim(1e2, 1e4)

ax[1,0].set_ylim(1e2, 1e5)
ax[1,1].set_ylim(1e2, 1e5)
ax[2,0].set_ylim(1e2, 1e5)
ax[2,1].set_ylim(1e2, 1e5)

# Extra customizations to BM fluorescence plots
for row_ax_ID in [1,2]:
    ax[row_ax_ID,1].set_xlim(start_aa,end_aa)
    ax[row_ax_ID,1].set_xlabel('')
    ax[row_ax_ID,1].set_yticklabels('')
    ax[row_ax_ID,1].set_ylabel('')
ax[1,1].set_xticklabels('')


# TODO: Customize legend for fluorescence
# frameon=False removes the frame where the labels were in
# Get the hue category label and removes it    
handles, labels = ax[2,1].get_legend_handles_labels()
# labels = [
#         '',
#         '+ DMSO',
#         '+ 10 μM 4-HT'
#         ]
labels = [
        '',
        # '+ DMSO',
        # '+ 10 μM 4-HT'
        'no induction',
        '1 mM arabinose (N-lobe)',
        '25 μM DAPG (C-lobe)',
        '1 mM arabinose + 25 μM DAPG (N+C lobes)'        
        ]
ax[2,1].legend(
        handles=handles[1:], labels=labels[1:],
      bbox_to_anchor=(1.005, 1),frameon=False, ncol=1,
      prop={'size': 12})

# Draw transposition window for all plots on the 2nd column
for ax_ID in [0,1,2]:
    ax[ax_ID,1].axvline(n_trans_border, color='k', dashes=[2,2], linestyle = '--', zorder=5)
    ax[ax_ID,1].axvline(c_trans_border, color='k', dashes=[2,2], linestyle = '--', zorder=5)  

# Customize outer labels
ax[1,0].set_ylabel('5 hours')
ax[2,0].set_ylabel('24 hours')
ax[2,1].set_xlabel('position (aa)')
ax[2,0].set_xticklabels(ctrl_tick_labels, rotation=40, fontsize=10, ha="right")


# TODO: Extra x ticks

minorLocator = AutoMinorLocator(10)

for row_ax_ID in [0,1,2]:
    
    xticks = ax[row_ax_ID,1].get_xticks() 
    xticks = np.append(xticks,start_aa)
    xticks = np.append(xticks,end_aa)
    
    xtickslist=xticks.tolist()
    xtickslist[-2]=str(start_aa)
    xtickslist[-1]=str(end_aa)
    ax[row_ax_ID,1].set_xticks(xticks)
    ax[2,1].set_xticklabels(xtickslist)
    ax[2,1].xaxis.set_major_formatter(FormatStrFormatter('%g'))
    
    # Plot minor ticks on x axis (optional depending on graph)
    ax[row_ax_ID,1].xaxis.set_minor_locator(minorLocator)
    
    ax[row_ax_ID,1].set_xlim(start_aa,end_aa)

#%%

## TODO: Extra y ticks
#for col_ax_ID in [0,1]:
#    for row_ax_ID in [1,2]:
#        
# Plot the 2nd structure of the protein being split (if available)
# All code below copied and modified directly on top of code from Biotite

# Code source: Patrick Kunzmann
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import biotite
import biotite.structure as struc
import biotite.structure.io.mmtf as mmtf
import biotite.sequence as seq
import biotite.sequence.graphics as graphics
import biotite.database.rcsb as rcsb
import biotite.sequence.io.genbank as gb
import biotite.database.entrez as entrez


# Create 'FeaturePlotter' subclasses
# for drawing the scondary structure features

class HelixPlotter(graphics.FeaturePlotter):

    def __init__(self):
        pass

    # Check whether this class is applicable for drawing a feature
    def matches(self, feature):
        if feature.key == "SecStr":
            if "sec_str_type" in feature.qual:
                if feature.qual["sec_str_type"] == "helix":
                    return True
        return False
    
    # The drawing function itself
    def draw(self, axes, feature, bbox, loc, style_param):
        # Approx. 1 turn per 3.6 residues to resemble natural helix
        n_turns = np.ceil((loc.last - loc.first + 1) / 3.6)
        x_val = np.linspace(0, n_turns * 2*np.pi, 100)
        # Curve ranges from 0.3 to 0.7
        y_val = (-0.4*np.sin(x_val) + 1) / 2
        
        # Transform values for correct location in feature map
        x_val *= bbox.width / (n_turns * 2*np.pi)
        x_val += bbox.x0
        y_val *= bbox.height
        y_val += bbox.y0
        
        # Draw white background to overlay the guiding line
        background = Rectangle(
            bbox.p0, bbox.width, bbox.height, color="white", linewidth=0
        )
        axes.add_patch(background)
        axes.plot(
            x_val, y_val, linewidth=1, color=biotite.colors["dimgreen"]
        )


class SheetPlotter(graphics.FeaturePlotter):

    def __init__(self, head_width=0.8, tail_width=0.5):
        self._head_width = head_width
        self._tail_width = tail_width


    def matches(self, feature):
        if feature.key == "SecStr":
            if "sec_str_type" in feature.qual:
                if feature.qual["sec_str_type"] == "sheet":
                    return True
        return False
    
    def draw(self, axes, feature, bbox, loc, style_param):
        x = bbox.x0
        y = bbox.y0 + bbox.height/2
        dx = bbox.width
        dy = 0
        
        if  loc.defect & seq.Location.Defect.MISS_RIGHT:
            # If the feature extends into the prevoius or next line
            # do not draw an arrow head
            draw_head = False
        else:
            draw_head = True
        
        axes.add_patch(biotite.AdaptiveFancyArrow(
            x, y, dx, dy,
            self._tail_width*bbox.height, self._head_width*bbox.height,
            # Create head with 90 degrees tip
            # -> head width/length ratio = 1/2
            head_ratio=0.5, draw_head=draw_head,
            color=biotite.colors["orange"], linewidth=0
        ))


#%%

ax[4,1] = plt.subplot(gs[4,1])

# Fetch GenBank files of the TK's first chain and extract annotatation
file_name = entrez.fetch("6FRH_A", biotite.temp_dir(), "gb", "protein", "gb")
gb_file = gb.GenBankFile()
gb_file.read(file_name)
annotation = gb.get_annotation(gb_file, include_only=["SecStr"])
# Length of the sequence
_, length, _, _, _, _ = gb.get_locus(gb_file)

graphics.plot_feature_map(
    ax[4,1], annotation, symbols_per_line=154,
    show_numbers=False, show_line_position=False,
    # 'loc_range' takes exclusive stop -> length+1 is required
    loc_range=(7,161),
    feature_plotters=[HelixPlotter(), SheetPlotter()]
)

ax[4,1].set_xlim(start_aa,end_aa)

# Make room for plot at the bottom
plt.gcf().subplots_adjust(hspace=0.15, wspace=0.1)

'''Save Figure'''
outputDir = os.path.join(dataRootDir,ssFolderDir1,figName)
fig.savefig(outputDir)