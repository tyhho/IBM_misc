# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 15:16:30 2020

@author: User
"""
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
from matplotlib import gridspec
from matplotlib import pyplot as plt



# TODO: Specify folder location
    # Each folder must contain only fcs files that end with well location
root_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), os.pardir))
data_folder = 'EXP012'
data_fn = 'IBM_EXP012_annotated_cyto_&_PR_data.json'
plot_inst_fn = 'IBM_EXP012R1_PlotOrder.csv'

fig_name = 'IBM_EXP012R1_Result.pdf'

data_file_path = os.path.join(root_path, data_folder, data_fn)
plot_inst_path = os.path.join(root_path, data_folder, plot_inst_fn)
fig_path = os.path.join(root_path, data_folder, fig_name)

imported_data = pd.read_json(data_file_path)
plot_data = imported_data[['SampleID', 'fluorescence (a.u.)', '4-HT (μM)']]
plot_inst = pd.read_csv(plot_inst_path)

plot_data = plot_data.merge(plot_inst)
plot_data.sort_values('order')

# data = plot_data['fluorescence (a.u.)']
data = np.array(plot_data['fluorescence (a.u.)'].values)


#%% Plot graphs

sns.set(style='ticks')
fig, ax = plt.subplots(1, 14, figsize=(4,2),
                       sharey=True, dpi = 600)
gs = gridspec.GridSpec(1, 14, wspace=0)

# fluo_values = data[10]

for i, fluo_values in enumerate(data):

    fluo_values = np.array(fluo_values)
    log_min, log_max = np.log10(fluo_values.min()), np.log10(fluo_values.max())
    bins = np.logspace(log_min, log_max, 100)
    
    ax[i] = plt.subplot(gs[i])
    sns.distplot(fluo_values, bins=bins, hist=True, kde=False, ax=ax[i],
                 color='k',
                 hist_kws=dict(edgecolor="k", linewidth=0),
                 vertical=True
                 )
    ax[i].set(yscale="log")
    
    xlim_min, xlim_max = ax[i].get_xlim()
    ax[i].set_xlim(xlim_max, xlim_min)
    ax[i].set_ylim(1e2, 5e4)
    ax[i].get_xaxis().set_visible(False)
    
    if i !=0:
        ax[i].axes.get_yaxis().set_visible(False)
    
# fig.tight_layout()

fig.savefig(fig_path)
