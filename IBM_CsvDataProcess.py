# -*- coding: utf-8 -*-

"""
@author: Trevor Ho
"""

import os
import pandas as pd

# TODO: Specify folder location

dataRootDir = r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir = 'FC023'
datafilename = 'IBM_FC023R3-5_FCmedian&metadata&PRData.csv'
csvDir = os.path.join(dataRootDir,dataFolderDir,datafilename)
data = pd.read_csv(csvDir,index_col=0)


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

data['Induction'] = data.apply(lambda data: processInduction(data['arabinose (mM)'], data['DAPG (μM)']), axis=1)

data = data.drop(['arabinose (mM)','DAPG (μM)'], axis=1)

outputFilename = 'IBM_FC023R3-5_FCmedian&metadata&PRData_InductionProcessed.csv'
outputDir = os.path.join(dataRootDir,dataFolderDir,outputFilename)
data.to_csv(outputDir)