# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 16:22:08 2020

@author: Trevor Ho
"""

# Fetch and load structure
# file_name = rcsb.fetch("1QGD", "mmtf", biotite.temp_dir())

import biotite.structure.io.pdb as pdb
import os
import biotite
import biotite.structure as struc
import biotite.sequence as seq
import biotite.sequence.graphics as graphics
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


dataRootDir=r'W:\Data storage & Projects\PhD Project_Trevor Ho\3_Intein-assisted Bisection Mapping'
dataFolderDir='BM010\ECF20_structure_model'

pdb_fn = "ECF20.pdb"
pdb_path = os.path.join(dataRootDir,dataFolderDir,pdb_fn)

pdb_file = pdb.PDBFile()
pdb_file.read(pdb_path)

array = pdb_file.get_structure(model=1)


#%%



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
            x_val, y_val, linewidth=2, color=biotite.colors["dimgreen"]
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

# Test our drawing functions with example annotation
annotation = seq.Annotation([
    seq.Feature("SecStr", [seq.Location(10, 40)], {"sec_str_type" : "helix"}),
    seq.Feature("SecStr", [seq.Location(60, 90)], {"sec_str_type" : "sheet"}),
])

#     def _add_sec_str(annotation, first, last, str_type):
#         if str_type == "a":
#             str_type = "helix"
#         elif str_type == "b":
#             str_type = "sheet"
#         else:
#             # coil
#             return
#         feature = seq.Feature(
#             "SecStr", [seq.Location(first, last)], {"sec_str_type" : str_type}
#         )
#         annotation.add_feature(feature)


fig = plt.figure(figsize=(8.0, 0.8))
ax = fig.add_subplot(111)
graphics.plot_feature_map(
    ax, annotation, multi_line=False, loc_range=(1,100),
    # Register our drawing functions
    feature_plotters=[HelixPlotter(), SheetPlotter()]
)

# The chain ID corresponding to each residue
# chain_id_per_res = array.chain_id[struc.get_residue_starts(array)]
# sse = pdb_file["secStructList"]

#%%
# sse = sse[sse != -1]
# sse = sse[chain_id_per_res == "A"]
# sse = np.array([sec_struct_codes[code] for code in sse if code != -1],
#                dtype="U1")
# sse = np.array([dssp_to_abc[e] for e in sse], dtype="U1")

# # Helper function to convert secondary structure array to annotation
# # and visualize it
# def visualize_secondary_structure(sse, first_id):

#     def _add_sec_str(annotation, first, last, str_type):
#         if str_type == "a":
#             str_type = "helix"
#         elif str_type == "b":
#             str_type = "sheet"
#         else:
#             # coil
#             return
#         feature = seq.Feature(
#             "SecStr", [seq.Location(first, last)], {"sec_str_type" : str_type}
#         )
#         annotation.add_feature(feature)

#     # Find the intervals for each secondary structure element
#     # and add to annotation
#     annotation = seq.Annotation()
#     curr_sse = None
#     curr_start = None
#     for i in range(len(sse)):
#         if curr_start is None:
#             curr_start = i
#             curr_sse = sse[i]
#         else:
#             if sse[i] != sse[i-1]:
#                 _add_sec_str(
#                     annotation, curr_start+first_id, i-1+first_id, curr_sse
#                 )
#                 curr_start = i
#                 curr_sse = sse[i]
#     # Add last secondary structure element to annotation
#     _add_sec_str(annotation, curr_start+first_id, i-1+first_id, curr_sse)

#     fig = plt.figure(figsize=(8.0, 3.0))
#     ax = fig.add_subplot(111)
#     graphics.plot_feature_map(
#         ax, annotation, symbols_per_line=150,
#         loc_range=(first_id, first_id+len(sse)),
#         show_numbers=True, show_line_position=True,
#         feature_plotters=[HelixPlotter(), SheetPlotter()]
#     )