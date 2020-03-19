# -*- coding: utf-8 -*-

"""
@author: Trevor Ho
"""

import os
import pandas as pd

# TODO: Specify folder location

dataRootDir = r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir = 'FC034'
dataFN = 'IBM_FC034R1-4_median&metadata&PRData.csv'
outputFN = dataFN.split('.csv')[0] + '_InductionRelabelled.csv'

dataFP = os.path.join(dataRootDir,dataFolderDir,dataFN)
data = pd.read_csv(dataFP,index_col=0)

#%%
def processInduction(ara,dapg):
    if ara == 0 and dapg == 0:
        induction = 'no induction'
    elif ara == 1 and dapg == 0:
        induction = '1 mM arabinose'
    elif ara == 0 and dapg == 25:
        induction = '25 μM DAPG'
    elif ara == 1 and dapg == 25:
        induction = '1 mM arabinose + 25 μM DAPG'
    return induction

# def processInduction(ht):
#     induction = ''
#     if ht == 0:
#         induction = '+ DMSO'
#     elif ht == 10:
#         induction = '+ 10 μM 4-HT'
#     return induction

data['Induction'] = data.apply(lambda data: processInduction(data['arabinose (mM)'], data['DAPG (μM)']), axis=1)
data = data.drop(['arabinose (mM)','DAPG (μM)'], axis=1)

# data['Induction'] = data.apply(lambda data: processInduction(data['4-HT (μM)']), axis=1)
# data = data.drop(['arabinose (mM)','4-HT (μM)'], axis=1)

outputDataFP = os.path.join(dataRootDir,dataFolderDir,outputFN)
data.to_csv(outputDataFP)


