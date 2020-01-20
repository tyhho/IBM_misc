# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 18:04:33 2019

@author: Trevor Ho
"""

import os
import pandas as pd


# TODO: Specify folder location
    # csv file containing data to be plotted should be in this folder

dataRootDir=r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir='FC029'
dataFN = 'IBM_FC029R4,5,7_median&metadata&PRData.csv'
hlineInfoFN = dataFN.split('.csv')[0] + '_hlineInfo.csv'

# Read data
dataFP = os.path.join(dataRootDir,dataFolderDir,dataFN)
data = pd.read_csv(dataFP,index_col=0)

# TODO: Input information for constructing hlines

# hline_input instruction:
    # key = str for df.query. Note that str inside str are double-quoted
    # value = list of colors to use, 1st color = line color, 2nd color = shade color
hline_input = {
        'SampleID == "IBMc307"': ['#F96495','#EDA3BA'],
        'SampleID == "IBMc101"': ['k','#E8E8E8']
        }

#%% Calculate stats for hlines
pi_dict = {0:5, 1:24}
hlineInfo = pd.DataFrame(columns=[])

for key, pi_time in pi_dict.items():
    for queryLine, color in hline_input.items():
        hRefSubData = data.query(queryLine)
        hRefSubData = hRefSubData[hRefSubData['Post-induction (hrs)']==pi_time]
        hRefSubData = hRefSubData['median fluorescence (a.u.)']        
        hline_info_row_dict = {
                'Post-induction (hrs)': [pi_time],
                'mean of median fluorescence (a.u.)': [hRefSubData.mean()],
                'std of median fluorescence (a.u.)': [hRefSubData.std()],
                'span for std': [(hRefSubData.std() > hRefSubData.mean()*0.1)],
                'line color': [color[0]],
                'shade color': [color[1]]
                }
        hlineInfo = hlineInfo.append(pd.DataFrame.from_dict(hline_info_row_dict))

#%%
hlineInfoFP = os.path.join(dataRootDir,dataFolderDir,hlineInfoFN)
hlineInfo.to_csv(hlineInfoFP)