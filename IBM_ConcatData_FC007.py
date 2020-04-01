# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 01:17:32 2020

@author: Trevor Ho
"""

import os
import pandas as pd

root_path = '..'
data_folder = 'FC007'

json_data_fns = [
    'IBM_FC007R1-3_cyto_&_PR_data.json',
    'IBM_FC007R4,5,7_cyto_&_PR_data.json'
    ]

json_output_fn = 'IBM_FC007R1-5,7_cyto_&_PR_data.json'

json_data_paths = [
    os.path.join(root_path, data_folder, json_data_fn)
    for json_data_fn in json_data_fns 
    ]

json_data = [
    pd.read_json(json_data_path) for json_data_path in json_data_paths
    ]

concat_data = pd.DataFrame()

for data in json_data:
    concat_data = concat_data.append(data, ignore_index=True)
    
json_output_path = os.path.join(root_path, data_folder, json_output_fn)
concat_data.to_json(json_output_path)

csv_data = concat_data.drop(['FC_fluorescence (a.u.)'], axis=1)
csv_output_path = json_output_path.split('.json')[0] + '.csv'
csv_data.to_csv(csv_output_path)