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
dataFolderDir='FC018'
dataFN = 'IBM_FC018R5-7_FCmedian&metadata&PRData.csv'
hlineDataFN = 'IBM_FC018R5-7_FCmedian&metadata&PRData_CtrlStats.csv'

# Read data
csvDir = os.path.join(dataRootDir,dataFolderDir,dataFN)
data = pd.read_csv(csvDir,index_col=0)


# TODO: Input information for constructing hlines

# hline_input instruction:
    # key = str for df.query. Note that str inside str are double-quoted
    # value = list of colors to use, 1st color = line color, 2nd color = shade color
hline_input = {
#        'SampleID == "IBMc186 + IBMc101"': ['k','#E8E8E8'],
        'SampleID == "IBMc120 + IBMc101"': ['k','#E8E8E8']
        }

#%% Calculate stats for hlines
pi_dict = {0:5, 1:24}
hlineData = pd.DataFrame(columns=[])

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
        hlineData = hlineData.append(pd.DataFrame.from_dict(hline_data_row_dict))

#%%
hlineDataDir = os.path.join(dataRootDir,dataFolderDir,hlineDataFN)
hlineData.to_csv(hlineDataDir)