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
ssFolderDir1 = 'BM005'
ssFolderDir2 = 'BM005_Sequencing Results'
ssDataFN = 'IBM_BM005_IdentifiedSplitSites.csv'

# TODO: Specify location of CSV file containing the median fluorescence values of the strains
fluoFolderDir='FC023'
fluoDataFN = 'IBM_FC023R3-5_FCmedian&metadata&PRData_InductionProcessed_FoldChangeFiltered.csv'

# TODO: Specify location of CSV file containing the median fluorescence values of control strains
ctrlFluoDataFN = 'IBM_FC023R3-5_FCmedianStats&metadata_FoldChangeFiltered.csv'

# TODO: Specify filename to be saved. Figure will be saved under path of dataRootDir\ssFolderDir1
figName = 'IBM_BM005_BisectionMap_Raw.pdf'

# TODO: Specify filename of csv with merged data to be saved (final data saves under BM005 directly)
mergedDataFN = 'IBM_BM005_PooledResults.csv'


# Preparations

# TODO: Set induction and induction time information
inductions = ['no induction','1 mM arabinose','25 μM DAPG','1 mM arabinose + 25 μM DAPG']
pi_times = [5,24]

# TODO: Update control list
control_list = ['IBMc252 + IBMc249',
        'IBMc259 + IBMc249',
        'IBMc259 + IBMc226']

# TODO: Set x tick labels of controls
ctrl_tick_labels = ['background','reporter','TetR + reporter']

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

# TODO: Define how many amino acids are being plotted
start_aa = 1
end_aa = 207

# TODO: Input transposition window
n_trans_border = 21
c_trans_border = 188

# TODO: Decide color for inductions
color_mapper = {
        'no induction':'#F96495',
        '1 mM arabinose':'#55A0FB',
        '25 μM DAPG':'#099963',
        '1 mM arabinose + 25 μM DAPG':'#1D1D1B'
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
mergedDataSingleCon=mergedData[(mergedData['Run']==3) \
                               & (mergedData['Induction']=='1 mM arabinose') \
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
#mergedData.to_csv(mergedDataDir)
            
#%% Calculate stats for hlines
pi_dict = {1:5, 2:24}
hline_data = pd.DataFrame(columns=[])

for _, pi_time in pi_dict.items():
    for queryLine, color in hline_input.items():
        hRefSubData = fluoData.query(queryLine)
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


# Import CSV containing control data for plotting
ctrlFluoDataDir = os.path.join(dataRootDir,fluoFolderDir,ctrlFluoDataFN)
ctrlFluoData = pd.read_csv(ctrlFluoDataDir,index_col=0)
ctrlFluoData = ctrlFluoData[ctrlFluoData['SampleID'].isin(control_list)]
ctrlFluoData['SortIndex'] = ctrlFluoData.apply(lambda data: control_list.index(data['SampleID'])+1, axis=1)
#%%
# Create canvas
sns.set(style='ticks')
fig, ax = plt.subplots(5, 2, sharex=True, sharey=True, figsize=(6, 6), dpi = 200)
# Grid spec
ctrl_len = len(control_list)
gs = gridspec.GridSpec(5, 2, width_ratios=[ctrl_len+10, end_aa], height_ratios = [0.6,1,1,0.2,0.2]) 

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
        ax[row_ax_ID,col_ax_ID].set_ylim(1e2, 3e5)
        ax[row_ax_ID,col_ax_ID].set(yscale="log")
        ax[row_ax_ID,col_ax_ID].yaxis.set_ticks(np.logspace(2, 5, num=4))
        locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
        ax[row_ax_ID,col_ax_ID].yaxis.set_minor_locator(locmin)       

        
        # Plot horizontal lines
        hline_data_to_plot = hline_data[hline_data['Post-induction (hrs)'] == pi_time]
        
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
labels = [
        '',
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

## TODO: Extra y ticks
#for col_ax_ID in [0,1]:
#    for row_ax_ID in [1,2]:
#        
#%% Plot the 2nd structure of the protein being split (if available)
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

length = 207
# Dictionary to convert 'secStructList' codes to DSSP values
# https://github.com/rcsb/mmtf/blob/master/spec.md#secstructlist
sec_struct_codes = {0 : "I",
                    1 : "S",
                    2 : "H",
                    3 : "E",
                    4 : "G",
                    5 : "B",
                    6 : "T",
                    7 : "C"}
# Converter for the DSSP secondary structure elements
# to the classical ones
dssp_to_abc = {"I" : "c",
               "S" : "c",
               "H" : "a",
               "E" : "b",
               "G" : "c",
               "B" : "b",
               "T" : "c",
               "C" : "c"}


# Fetch and load structure
file_name = rcsb.fetch("4AC0", "mmtf", biotite.temp_dir())
mmtf_file = mmtf.MMTFFile()
mmtf_file.read(file_name)
array = mmtf.get_structure(mmtf_file, model=1)

# homodimer
tetR_dimer = array[struc.filter_amino_acids(array)]
# monomer
tetR_mono = tetR_dimer[tetR_dimer.chain_id == "A"]

# The chain ID corresponding to each residue
chain_id_per_res = array.chain_id[struc.get_residue_starts(tetR_dimer)]
sse = mmtf_file["secStructList"]
sse = sse[sse != -1]
sse = sse[chain_id_per_res == "A"]
sse = np.array([sec_struct_codes[code] for code in sse if code != -1],
               dtype="U1")
sse = np.array([dssp_to_abc[e] for e in sse], dtype="U1")

# convert secondary structure array to annotation and visualize it
# Visualize seconday structure array
# Sine the residues may not start at 1,
# provide the actual first residue ID


first_id = tetR_mono.res_id[0]
def _add_sec_str(annotation, first, last, str_type):
    if str_type == "a":
        str_type = "helix"
    elif str_type == "b":
        str_type = "sheet"
    else:
        # coil
        return
    feature = seq.Feature(
        "SecStr", [seq.Location(first, last)], {"sec_str_type" : str_type}
    )
    annotation.add_feature(feature)

# Find the intervals for each secondary structure element
# and add to annotation
annotation = seq.Annotation()
curr_sse = None
curr_start = None
for i in range(len(sse)):
    if curr_start is None:
        curr_start = i
        curr_sse = sse[i]
    else:
        if sse[i] != sse[i-1]:
            _add_sec_str(
                annotation, curr_start+first_id, i-1+first_id, curr_sse
            )
            curr_start = i
            curr_sse = sse[i]
# Add last secondary structure element to annotation
_add_sec_str(annotation, curr_start+first_id, i-1+first_id, curr_sse)

#fig = plt.figure(figsize=(8.0, 1.0))
#ax = fig.add_subplot(111)
ax[4,1] = plt.subplot(gs[4,1])
graphics.plot_feature_map(
    ax[4,1], annotation, multi_line=False, loc_range=(1,length+1),
    show_numbers=False, show_line_position=False,
    feature_plotters=[HelixPlotter(), SheetPlotter()]
)

ax[4,1].set_xlim(start_aa,end_aa)

#fig.tight_layout()


# Make room for plot at the bottom
plt.gcf().subplots_adjust(hspace=0.15, wspace=0.1)

'''Save Figure'''
outputDir = os.path.join(dataRootDir,ssFolderDir1,figName)
fig.savefig(outputDir)
