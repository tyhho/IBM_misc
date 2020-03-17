# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 14:55:50 2020

@author: Trevor Ho

Provide functions to plot the 2nd structure of the protein being split (if available)
All code below copied and modified directly on top of code from Biotite

Also provide functions to convert csv files containing secondary structure
information into biotite.sequence.annotation objects for plotting

# Code source: Patrick Kunzmann
# License: BSD 3 clause

"""

import numpy as np
from matplotlib.patches import Rectangle
import biotite
import biotite.structure as struc
import biotite.sequence as seq
import biotite.sequence.io.genbank as gb
import biotite.database.entrez as entrez
import biotite.sequence.graphics as graphics

import pandas as pd

# Create 'FeaturePlotter' subclasses
# for drawing the scondary structure features

class HelixPlotter(graphics.FeaturePlotter):

    def __init__(self):
        pass

    # Check whether this class is applicable for drawing a feature
    def matches(self, feature):
        if feature.key == "SecStr":
            if "sec_str_type" in feature.qual:
                if feature.qual["sec_str_type"] == "helix":
                    return True
        return False
    
    # The drawing function itself
    def draw(self, axes, feature, bbox, loc, style_param):
        # Approx. 1 turn per 3.6 residues to resemble natural helix
        n_turns = np.ceil((loc.last - loc.first + 1) / 3.6)
        x_val = np.linspace(0, n_turns * 2*np.pi, 100)
        # Curve ranges from 0.3 to 0.7
        y_val = (-0.4*np.sin(x_val) + 1) / 2
        
        # Transform values for correct location in feature map
        x_val *= bbox.width / (n_turns * 2*np.pi)
        x_val += bbox.x0
        y_val *= bbox.height
        y_val += bbox.y0
        
        # Draw white background to overlay the guiding line
        background = Rectangle(
            bbox.p0, bbox.width, bbox.height, color="white", linewidth=0
        )
        axes.add_patch(background)
        axes.plot(
            x_val, y_val, linewidth=1, color=biotite.colors["dimgreen"]
        )


class SheetPlotter(graphics.FeaturePlotter):

    def __init__(self, head_width=0.8, tail_width=0.5):
        self._head_width = head_width
        self._tail_width = tail_width


    def matches(self, feature):
        if feature.key == "SecStr":
            if "sec_str_type" in feature.qual:
                if feature.qual["sec_str_type"] == "sheet":
                    return True
        return False
    
    def draw(self, axes, feature, bbox, loc, style_param):
        x = bbox.x0
        y = bbox.y0 + bbox.height/2
        dx = bbox.width
        dy = 0
        
        if  loc.defect & seq.Location.Defect.MISS_RIGHT:
            # If the feature extends into the prevoius or next line
            # do not draw an arrow head
            draw_head = False
        else:
            draw_head = True
        
        axes.add_patch(biotite.AdaptiveFancyArrow(
            x, y, dx, dy,
            self._tail_width*bbox.height, self._head_width*bbox.height,
            # Create head with 90 degrees tip
            # -> head width/length ratio = 1/2
            head_ratio=0.5, draw_head=draw_head,
            color=biotite.colors["orange"], linewidth=0
        ))

def ss_csv_to_annotation(csv_path=str):

    # Codes retained for debugging
    # dataRootDir=r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
    # dataFolderDir='BM010\ECF20_structure_model'
    # exported_ss = pd.read_csv(os.path.join(dataRootDir,dataFolderDir,'ECF20_ExPASy_sec_struct.csv'))
    
    exported_ss = pd.read_csv(csv_path)
    ss_segments = pd.DataFrame()
    
    # Take info for the first ss segment without knowing when it ends
    
    start_aa, last_ss = exported_ss.iloc[0]
    previous_ss_seg = last_ss
    
    seq_end_aa, _ = exported_ss.iloc[-1] # for recording the last segment of ss
    
    for _, aa, ss in exported_ss.itertuples():
    
        # Only when a new ss sgement is detected would an entry
        # for the previous ss segment be recorded
        
        if ss != last_ss:
            
            ss_unit_entry = pd.DataFrame({
                            'start_aa': [start_aa],
                            'end_aa': [aa - 1],
                            'sec_str_type': [previous_ss_seg]
                            })
            ss_segments = ss_segments.append(ss_unit_entry)
            
            previous_ss_seg = ss
            start_aa = aa
    
        if aa == seq_end_aa:
            ss_unit_entry = pd.DataFrame({
                            'start_aa': [start_aa],
                            'end_aa': [aa],
                            'sec_str_type': [previous_ss_seg]
                            })
            ss_segments = ss_segments.append(ss_unit_entry)
        
        last_ss = ss
    
    # At this point the df ss_segments also contains 'L' linkers
    # Create new df to store only those relevent for plotting
    
    ss_segments_plot = ss_segments.query('sec_str_type != "L"')
    
    annotation = seq.Annotation()
    for _, start_aa, end_aa, ss_type in ss_segments_plot.itertuples():
        if ss_type == "H":
            ss_type = "helix"
        elif ss_type == "S":
            ss_type = "sheet"
        
        feature = seq.Feature(
                "SecStr", [seq.Location(start_aa, end_aa)], {"sec_str_type" : ss_type}
            )
        annotation.add_feature(feature)

    return annotation

def fetch_gb_annotation(pdb_chain=str):
    
    # input line retained for debugging
    # pdb_chain = "6FRH_A"
    
    # Fetch GenBank files of the TK's first chain and extract annotatation
    file_name = entrez.fetch(pdb_chain, biotite.temp_dir(), "gb", "protein", "gb")
    gb_file = gb.GenBankFile()
    gb_file.read(file_name)
    annotation = gb.get_annotation(gb_file, include_only=["SecStr"])
    # Length of the sequence
    _, length, _, _, _, _ = gb.get_locus(gb_file)