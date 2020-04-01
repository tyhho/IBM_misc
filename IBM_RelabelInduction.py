# -*- coding: utf-8 -*-

"""
@author: Trevor Ho
"""

import os
import pandas as pd

#%%

# def processInduction(ara,dapg):
#     if ara == 0 and dapg == 0:
#         induction = 'no induction'
#     elif ara == 1 and dapg == 0:
#         induction = '1 mM arabinose'
#     elif ara == 0 and dapg == 25:
#         induction = '25 μM DAPG'
#     elif ara == 1 and dapg == 25:
#         induction = '1 mM arabinose + 25 μM DAPG'
#     return induction

def processInduction(ara, ahl):
    if ara == 0 and ahl == 0:
        induction = 'no induction'
    elif ara == 1 and ahl == 0:
        induction = '1 mM arabinose'
    elif ara == 0 and ahl == 1:
        induction = '1 μM AHL'
    elif ara == 1 and ahl == 1:
        induction = '1 mM arabinose + 1 μM AHL'
    return induction


# def processInduction(ht):
#     induction = ''
#     if ht == 0:
#         induction = '+ DMSO'
#     elif ht == 10:
#         induction = '+ 10 μM 4-HT'
#     return induction

#%%

# TODO: Specify folder location
root_path = '..'
data_folder = 'FC007'
json_data_fn = 'IBM_FC007R1-5,7_cyto_&_PR_data.json'
json_data_path = os.path.join(root_path, data_folder, json_data_fn) 
json_data = pd.read_json(json_data_path)

json_output_fn = json_data_fn.split('.json')[0] + '_InductionRelabelled.json'
json_output_path = os.path.join(root_path, data_folder, json_output_fn) 

# json_data['Induction'] = json_data.apply(
#     lambda json_data: processInduction(json_data['arabinose (mM)'], json_data['DAPG (μM)'])
#     , axis=1)
# json_data = json_data.drop(['arabinose (mM)','DAPG (μM)'], axis=1)

json_data['Induction'] = json_data.apply(
    lambda json_data: processInduction(json_data['arabinose (mM)'], json_data['AHL (μM)'])
    , axis=1)
json_data = json_data.drop(['arabinose (mM)','AHL (μM)'], axis=1)


# data['Induction'] = data.apply(lambda data: processInduction(data['4-HT (μM)']), axis=1)
# data = data.drop(['arabinose (mM)','4-HT (μM)'], axis=1)

json_data.to_json(json_output_path)

csv_data = json_data.drop(['FC_fluorescence (a.u.)'], axis=1)
csv_output_path = json_output_path.split('.json')[0] + '.csv'
csv_data.to_csv(csv_output_path)