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
from matplotlib.ticker import AutoMinorLocator

sns.reset_defaults

# Specify root dir
dataRootDir = r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'

# TODO: Specify location of CSV file containing identified split sites
ssFolderDir1 = 'BM005'
ssFolderDir2 = 'BM005_Sequencing Results'
ss_csv = 'IBM_BM005_IdentifiedSplitSites.csv'

# TODO: Specify location of CSV file containing the median fluorescence values of the strains
fluoFolderDir='FC023'
fluo_csv = 'IBM_FC023R3-5_FCmedian&metadata&PRData_InductionProcessed.csv'
# TODO: Specify location of CSV file containing the median fluorescence values of control strains
control_fluo_csv = 'IBM_FC023R3-5_FCmedianStats&metadata.csv'

# TODO: Specify filename to be saved. Figure will be saved under path of dataRootDir\ssFolderDir1
figName = 'IBM_BM005_BisectionMap_Raw.eps'

'''Read data and map split sites to function'''
# Read split site csv file
ss_csvDir = os.path.join(dataRootDir,ssFolderDir1,ssFolderDir2,ss_csv)
ss_data = pd.read_csv(ss_csvDir,index_col=0)
ss_data = ss_data[['aa_ss_middle']]
ss_data['SampleID']=ss_data.index

# Read fluorescence csv file
fluo_csvDir = os.path.join(dataRootDir,fluoFolderDir,fluo_csv)
fluo_data = pd.read_csv(fluo_csvDir,index_col=0)

# Map split site to fluorescence
mergedData = fluo_data.merge(ss_data,sort=False)
# Generate a stripped mergedData set for histogram plot
mergedDataSingleCon=mergedData[(mergedData['Run']==3) \
                               & (mergedData['Induction']=='1 mM arabinose') \
                               & (mergedData['Post-induction (hrs)']==5)]
mergedDataSingleCon=mergedDataSingleCon['aa_ss_middle']

'''Loop through the amino acid split sites and consolidate median fluorescence'''
# Preparations
inductions = ['no induction','1 mM arabinose','25 μM DAPG','1 mM arabinose + 25 μM DAPG']
pi_times = [5,24]

# Deduplicate the list of amino acid split sites and sort to a list
unique_ss_list = ss_data['aa_ss_middle'].drop_duplicates().sort_values().tolist()

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


#%%

'''Plot the graphs'''

# TODO: Define how many amino acids are being plotted
start_aa = 1
end_aa = 200

# Import CSV containing control data
control_fluo_csvDir = os.path.join(dataRootDir,fluoFolderDir,control_fluo_csv)
control_fluo_data = pd.read_csv(control_fluo_csvDir,index_col=0)

# TODO: Check control data to be used
# Keep control data only
control_list = ['IBMc252 + IBMc249',
    'IBMc259 + IBMc249']
control_fluo_data = control_fluo_data[control_fluo_data['SampleID'].isin(control_list)]

# Create canvas
sns.set(style='ticks')
fig, ax = plt.subplots(3, 1, sharex=True, figsize=(6, 8), dpi = 100)

# Plot the counts of strains being mapped to the same split site
sns.distplot(mergedDataSingleCon, kde=False, color="k", ax=ax[0], bins=end_aa)
ax[0].set_xlabel('')
ax[0].set_ylabel('sequencing \n count')
ax[0].tick_params(axis='y', which='minor', left=True)
minorLocator = AutoMinorLocator(5)
ax[0].yaxis.set_minor_locator(minorLocator)

# Plot fluorescence of split sites

pi_dict = {1:5, 2:24}
for ax_ID, pi_time in pi_dict.items():
    
    # Create df to isolate data of interest
    subdata = unique_ss_data[unique_ss_data['Post-induction (hrs)']==pi_time]
    
    sns.scatterplot(x='position (aa)', y='mean of median fluorescence (a.u.)', 
                       hue_order=['no induction','1 mM arabinose','25 μM DAPG','1 mM arabinose + 25 μM DAPG'],
                       hue='Induction',
                       #palette=['#808080','r'],
                       data=subdata, ax=ax[ax_ID],zorder=10,
                       edgecolor=None)
    ax[ax_ID].set(yscale="log")

    # Plot error bars
    errbar_colors = {'1 mM arabinose': '#808080',
                     '1 mM arabinose + 0.1 mM caffeine': 'r'}
    for induction in inductions:
        sub_unique_ss_data = subdata[subdata['Induction']==induction]
        ssVec = sub_unique_ss_data['position (aa)'].tolist()
        meanVec = sub_unique_ss_data['mean of median fluorescence (a.u.)'].tolist()
        errVec = sub_unique_ss_data['std of median fluorescence (a.u.)'].tolist()
        ax[ax_ID].errorbar(x=ssVec, y=meanVec, yerr = errVec, fmt='none',
          zorder=6
#          , color=errbar_colors[induction]
          )
    
    # TODO: Check range of identified split sites & range of fluorescence
    ax[ax_ID].set_ylim(1e0, 6e5)
    ax[ax_ID].set_xlim(start_aa,end_aa)

    # TODO: Check transposition window location
    n_trans_border = 6
    ax[ax_ID].axvline(n_trans_border, color='k', dashes=[2,2], linestyle = '--', zorder=5)
    # also add line for count plot
    ax[0].axvline(n_trans_border, color='k', dashes=[2,2], linestyle = '--', zorder=5)
    
    # Plot minor ticks on x axis (optional depending on graph)
    minorLocator = AutoMinorLocator(5)
    ax[ax_ID].xaxis.set_minor_locator(minorLocator)

    # Calculate values of control data for horizontal line plotting
    sub_control_data = control_fluo_data[control_fluo_data['Post-induction (hrs)']==pi_time]
    negdata = sub_control_data[(sub_control_data['SampleID']==control_list[0])]
    negFluos = negdata['mean of median fluorescence (a.u.)']
    neg_mean = negFluos.mean()
    neg_sd = negFluos.std()
    
    posdata = sub_control_data[sub_control_data['SampleID']==control_list[1]]
    posFluos = posdata['mean of median fluorescence (a.u.)']
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
    
# Customize y labels of the two different graphs
ax[1].set_ylabel('5 hours')
ax[2].set_ylabel('24 hours')

# Customize legend for fluorescence
# frameon=False removes the frame where the labels were in
# Get the hue category label and removes it    
handles, labels = ax[1].get_legend_handles_labels()
ax[1].legend(handles=handles[1:], labels=labels[1:],
  loc='upper right',frameon=False)

#%%
# Make room for plot at the bottom
plt.gcf().subplots_adjust(bottom=0.2)

'''Save Figure'''
outputDir = os.path.join(dataRootDir,ssFolderDir1,figName)
fig.savefig(outputDir)