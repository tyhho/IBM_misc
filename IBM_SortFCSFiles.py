# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 21:56:24 2019

@author: Trevor Ho

This script sort FCS files into the folder according to their filenames
"""

import os
import shutil
import fnmatch

# TODO: Specify folder location
    # Each folder must contain only fcs files that end with well location
dataRootDir = r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir = 'FC018'

sourceFolder = 'FC018_awaitng sort FCS'

fn_general_pattern = '*_IBM_FC018*.fcs'

sourceFolderDir = os.path.join(dataRootDir,dataFolderDir,sourceFolder)
allFiles = os.listdir(sourceFolderDir)

# Get all files that match the search criteria

matchedFiles = fnmatch.filter(allFiles,fn_general_pattern)

fn_keylist = []
for filename in matchedFiles:
    
    fn_keyInfo = filename.split('_IBM_')[1].split('_Experiment')[0]
    sourceFileDir = os.path.join(dataRootDir,dataFolderDir,sourceFolder,filename)
    #fn_keylist.append(fn_keyInfo)
    
    destFolder ='IBM_' + fn_keyInfo + '_FCS'
    destFolderDir = os.path.join(dataRootDir,dataFolderDir,sourceFolder,destFolder)
    
    shutil.move(sourceFileDir,destFolderDir)

#
#
#print(fn_keylist)
