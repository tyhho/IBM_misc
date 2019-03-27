# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 16:46:43 2019

@author: Trevor Ho

Split Site Identification

This script identifies the split sites from sequencing results of the PCR
products using the Smith-Waterman local alignment algorithmn.

Summary of process
1. Import and read a .seq file into a sequence
2. Align the sequence to a signature sequence of the substitution insert
3. Extract the target sequence that was adjacent to the substitution insert
4. (Optional) Remove bases from the extracted sequence if they came from
    duplication by transposon
5. Align the extracted sequence to the full sequence of the target CDS sequence
6. Locate the location of the split site in DNA and amino acid sequence using
    information from 5
7. Export the data as a .csv file

"""

import os
from Bio import SeqIO, Align
from Bio.Seq import Seq
import pandas as pd

# TODO: Provide the folder path that leads to the sequence

# Specify folder location
    # all .seq files should be in this folder
    # this is also the folder where the output will be deposited
dataRootDir = r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir = 'BM004'
dataFolder2Dir = 'Sequencing Results'
#seqFile = '450185201_IBMs0161_IBMo0503_A01.ab1' # for debugging

# TODO: Provide the CDS sequence of the protein to align to
targetCDS_sequence = 'tgcatctcgggagatagtttgatcagcttggcgagcacagggaaaagagttcctattaaggatttgttaggcgaaaaagattttgaaatatgggcaattaatgaacagacgatgaagctggaatcagctaaagttagtcgtgtattttgtaccggcaaaaagctagtctatactctaaaaactcgactaggtagaactatcaaggcaacagcaaatcatagatttttaactattgatggttggaaaagattagatgagctatctttaaaagagcatattgctctaccccgtaaactagaaagctcctctttacaattggcaccagaaatagaaaagttgcctcagagtgatatttactgggaccccatcgtttctattacggagactggagtcgaagaggtttttgatttgactgtgccaggactacgtaactttgtcgccaatgacatcattgtacataac'

# TODO: Provide the signature sequence. Usually, 50 bp are more than enough
signature_sequence= 'CCCAGGTTACCGTTAGCCATGGTGGAGGTGGTTCAGGTGGCGGAGGTTCA'

#%% Execution
ss_output = pd.DataFrame(columns = [])  # set up empty df to store all data

# Parse all files in the same folder

folderDir = os.path.join(dataRootDir, dataFolderDir,dataFolder2Dir)

for file in os.listdir(folderDir):
     filename = os.fsdecode(file)
     if filename.endswith(".ab1"):
        ## Import the data & extract metadata from the sequenicng result
        seqDir = os.path.join(folderDir,filename)
        record = SeqIO.read(seqDir, "abi")
        seq_order,strainID,seq_Primer = record.id.split('_')
        
        ## Sequence alignment part
        
        # Set up aligner for alignment
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -1
        
        # Perform local alignment between the sequencing result and the signature sequence
        alignment = aligner.align(record.seq,signature_sequence)
        # Locate the indices where the signature sequence was aligned to
        sig_start,sig_end = alignment[0].path[0][0], alignment[0].path[1][0]    
        # Extract the bp that would correspond to the target CDS
        ex_seq = record.seq[sig_end:]   # Note: no need to +1 here because the indexing of strings starts from 0
        
        # Perform local alignment between the extracted sequence and the target CDS sequence
        targetCDS = Seq(targetCDS_sequence) # import the targetCDS sequence as a Seq.Seq object
        alignment = aligner.align(targetCDS,ex_seq)
        # Locate the start index of where the extracted sequence was aligned to, i.e. ID the split site
        dna_ss_minus1 = alignment[0].path[0][0] # Note: indexing of string starts from 0, so this becomes -1
        
        # Workout the rest of the split sites
        dna_ss_plus1 = dna_ss_minus1 + 1
        aa_ss_minus1 = dna_ss_minus1 / 3
        aa_ss_plus1 = aa_ss_minus1 + 1
        
        # Create a dictionary for import into pandas dataframe
        ss_data = {
                'order': seq_order,
                'primer': seq_Primer,
                'raw_seq': record.seq._data, # raw sequence
                'read_length': len(record.seq), # length of the sequence
                'dna_ss_minus1': dna_ss_minus1,
                'dna_ss_plus1':  dna_ss_minus1 + 1,
                'aa_ss_minus1': dna_ss_minus1 / 3,
                'aa_ss_plus1': dna_ss_minus1 / 3 + 1        
                }
        
        ss_dfRow = pd.DataFrame(ss_data,index=[strainID])
        ss_output = ss_output.append(ss_dfRow, sort=True)
#%% Export output to file
outputFilename = 'BM004_SplitSitesIDed.csv'
outputDir = os.path.join(folderDir,outputFilename)
ss_output.to_csv(outputDir)

#%% Reads the sequencing trace. Not actively looking into at the moment
#from matplotlib import pyplot as plt
#channels = ['DATA9', 'DATA10', 'DATA11', 'DATA12']
#from collections import defaultdict
#trace = defaultdict(list)
#for c in channels:
#    trace[c] = record.annotations['abif_raw'][c]
#
#plt.plot(trace['DATA9'], color='blue')
#plt.plot(trace['DATA10'], color='red')
#plt.plot(trace['DATA11'], color='green')
#plt.plot(trace['DATA12'], color='yellow')
#plt.ylim(500,3500)