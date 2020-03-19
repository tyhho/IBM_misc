# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 20:13:37 2020

@author: Trevor Ho
"""

import os
import pandas as pd

# Specify root dir
dataRootDir = r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir1 = 'BM011'

ss_file_path = os.path.join(dataRootDir,dataFolderDir1,'SrpR_structure_model','SrpR_Jpred_sec_struct.txt')

f=open(ss_file_path, "r")

predicted_sec_struct_str = f.read()



aa_positions = []
predicted_sec_structs = []

for aa_pos, predicted_sec_struct in enumerate(predicted_sec_struct_str, start = 1):
    aa_positions.append(aa_pos)
    if predicted_sec_struct == "E":
        predicted_sec_struct = "S"
    elif predicted_sec_struct == "-":
        predicted_sec_struct = "L"
    predicted_sec_structs.append(predicted_sec_struct)

ss_df = pd.DataFrame(
    {'aa': aa_positions,
     'ss': predicted_sec_structs
     }
    )

ss_outputfile_path = ss_file_path.split(".txt")[0] + ".csv"
ss_df.to_csv(ss_outputfile_path, index = False)