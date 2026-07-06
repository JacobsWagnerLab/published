# -*- coding: utf-8 -*-
"""
Peptide sequence analysis and feature extraction (Figure 1A-B; Figures S1A, S10B).

Defines the panel of antimicrobial peptides/peptoids used in this study (amino-
acid sequences, source organism, predicted secondary structure, C-terminal amidation,
disulfide bridges, MIC, reported mechanism of action), computes their
physicochemical properties via `get_protein_df` (GRAVY index, net charge,
molecular weight, polar/hydrophobic/cationic fractions, etc.; see
protein_analysis_functions.py), and produces:
  - the AMP summary table          (Figure 1A)
  - the property scatterplots       (Figure 1B)
  - the physicochemical-properties table (Figure S1A)
  - the sequences / class-and-type table (Figure S10B)

@author: alessio fragasso
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from matplotlib.ticker import MaxNLocator
import warnings
warnings.filterwarnings('ignore')

from protein_analysis_functions import *


'''Settings'''
# Output folder for the figures/tables. Leave as "" to only display them
# (no files written). Set to a path to also save PNG + PDF.
save_path = ""   # e.g. r"path_to_output_folder"


def save_figure(basename):
    """Save the current figure as PNG + PDF into save_path, if one is set."""
    if not save_path:
        return
    os.makedirs(save_path, exist_ok=True)
    plt.savefig(os.path.join(save_path, basename + '.png'), format='png', dpi=300)
    plt.savefig(os.path.join(save_path, basename + '.pdf'), format='pdf')


'''Global font settings'''
plt.rcParams['font.family'] = 'sans-serif'       # Use sans-serif globally
plt.rcParams['font.sans-serif'] = ['Arial']      # Specify Arial as the sans-serif font
plt.rcParams['mathtext.fontset'] = 'custom'      # Use custom settings for math
plt.rcParams['mathtext.rm'] = 'Arial'            # Regular math text in Arial
plt.rcParams['mathtext.it'] = 'Arial:italic'     # Italic math text in Arial
plt.rcParams['mathtext.bf'] = 'Arial:bold'       # Bold math text in Arial
mpl.rcParams['pdf.fonttype'] = 42                # Store text as TrueType in PDFs


'''Peptide panel definition'''
label = ['LL-37',
         'Cecropin A',
         'Magainin-2',
         'Human beta-defensin 3',
         'Indolicidin',
         'PR-39',
         'Bactenecin-7',
         'DJK-5',
         'Melittin',
         'Tachyplesin-1',
         'Tur1A',
         'CRAMP',
         'Camel Bactenecin',
         'Protegrin-1',
         'Tur1B'
         ]
abbr = ['LL-37',
        'CecA',
        'Mag2',
        'HBD-3',
        'Indo',
        'PR-39',
        'Bac7',
        'DJK-5',
        'Mltt',
        'Tac1',
        'Tur1A',
        'CRAMP',
        'CamBac',
        'PG-1',
        'Tur1B'
        ]
source = [
    'Human',
    'Cecropia Moth',
    'Xenopus laevis',
    'Human',
    'Bovine',
    'Porcine',
    'Bovine',
    'Synthetic',
    'Honey Bee',
    'Horseshoe crab',
    'Cetacean',
    'Mouse',
    'Camel',
    'Porcine',
    'Cetacean'
    ]
aa_seq = ['LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES',    # LL-37
          'KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK',    # CecA
          'GIGKFLHSAKKFGKAFVGEIMNS',                    # Mag2
          'GIINTLQKYYCRVRGGRCAVLSCLPKEEQIGKCSTRGRKCCRRKK',    # HBD3
          'ILPWKWPWWPWRR',                                      # Indo
          'RRRPRPPYLPRPRPPPFFPPRLPPRIPPGFPPRFPPRFP',      # PR-39
          'RRIRPRPPRLPRPRPRPLPFPRPGPRPIPRPLPFP',          # Bac7
          'VQWRAIRVRVIR',                               # DJK-5
          'GIGAVLKVLTTGLPALISWIKRKRQQ',         # Mltt
          'KWCFRVCYRGICYRRCR',                  # Tac1
          'RRIRFRPPYLPRPGRRPRFPPPFPIPRIPRIP',      # Tur1A
          'GLLRKGGEKIGEKLKKIGQKIKNFFQKLVPQPEQ',     # CRAMP
          'RSIRRPRLPRPRVPRPYIPPRIPRPVLPPPRFPIPRFPRGR',     # CamBac
          'RGGRLCYCRRRFCVCVGR',                  # PG-1
          'RRIPFWPPNWPGPWLPPWSPPDFRIPRILRKR'       # Tur1B
        ]
structure = [
    'alpha-helix',
    'alpha-helix with hinge domain',
    'alpha-helix',
    'beta-sheet',
    'extended',
    'extended',
    'extended',
    'unknown',
    'alpha-helix',
    'beta-hairpin',
    'extended',
    'alpha-helix',
    'extended',
    'beta-hairpin',
    'extended'
]
structure_ref = [
    'PDB ID: 5NMN',    # LL-37
    'AF-P01507-F1',    # CecA
    'PDB ID: 2MAG',    # Mag2
    'PDB ID: 1KJ6',    # HBD3
    'PDB ID: 1G8C',    # Indo
    'AF-P80054-F1',    # PR-39
    'PDB ID: 5HAU',    # Bac7
    '',                # DJK-5
    'PDB ID: 2MLT',    # Mltt
    'PDB ID: 2RTV',    # Tac1
    'PDB ID: 6FKR',    # Tur1A
    '',                # CRAMP
    '',                # CamBac
    'PDB ID: 1PG1',    # PG-1
    ''                 # Tur1B
]
moa = ['Membrane permeabilization, cytoplasmic rigidification',
       'Membrane permeabilization',
       'Membrane permeabilization',
       'Membrane permeabilization',
       'DNA binding and synthesis inhibition',
       'Inhibition of DNA and protein synthesis',
       'Inhibition of DNA and protein synthesis; ribosome binding',
       'Binding and inhibition of (p)ppGpp',
       'Membrane permeabilization',
       'Binding to DNA minor groove',
       'Binding to ribosome',
       'Membrane permeabilization',
       'Inhibition of DNA and protein synthesis; ribosome binding',
       'Membrane permeabilization',
       'Membrane permeabilization'
      ]
MIC = [1,     ## LL37
       0.5,   ## CecA
       40,    ## Mag2
       20,    ## HBD3
       12.5,  ## Indo
       2,     ## PR-39
       0.6,   ## Bac7
       2.5,   ## DJK5
       4,     ## Melittin
       0.6,   ## Tac1
       1,     ## Tur1A
       4,     ## CRAMP
       0.5,   ## CamBac
       1,     ## PG-1
       8      ## Tur1B
       ]
amidated = [0,   ## LL37
            1,   ## CecA
            0,   ## Mag2
            0,   ## HBD3
            1,   ## Indo
            1,   ## PR-39
            0,   ## Bac7
            1,   ## DJK5
            1,   ## Melittin
            1,   ## Tac1
            0,   ## Tur1A
            0,   ## CRAMP
            0,   ## CamBac
            1,   ## PG-1
            0    ## Tur1B
            ]
amidation = ['No',    ## LL37
             'Yes',   ## CecA
             'No',    ## Mag2
             'No',    ## HBD3
             'Yes',   ## Indo
             'Yes',   ## PR-39
             'No',    ## Bac7
             'Yes',   ## DJK5
             'Yes',   ## Melittin
             'Yes',   ## Tac1
             'No',    ## Tur1A
             'No',    ## CRAMP
             'No',    ## CamBac
             'Yes',   ## PG-1
             'No'     ## Tur1B
             ]
disulfide_bridges = [0,
                     0,
                     0,
                     3,   ## HBD3
                     0,
                     0,
                     0,
                     0,
                     0,
                     2,   ## Tac1
                     0,
                     0,   ## CRAMP
                     0,   ## CamBac
                     2,   ## PG-1
                     0    ## Tur1B
                     ]


'''Build the per-peptide property table'''
protein_df = get_protein_df(label, aa_seq, amidated, disulfide_bridges, source)
protein_df[r'MIC ($\mu$M)'] = MIC
protein_df['Reported Mechanism of Action'] = moa
structure_ref = [x + ' [' + y + ']' for x, y in zip(structure, structure_ref)]
protein_df['Structure'] = structure
protein_df['C-terminal NH2'] = amidation
protein_df['Molecular weight (kDa)'] = protein_df['Molecular weight (Da)'] / 1000
protein_df['Abbr.'] = abbr

protein_df['Name (abbreviation)'] = [s + ' (' + t + ')' for (s, t) in zip(protein_df['Name'], protein_df['Abbr.'])]
protein_df = protein_df.sort_values(by='Name (abbreviation)', ascending=True)

# Palette of alternating row background colors for the tables
palette = [
    "#f2ede4",  # light cream/brown
    "#d8d1c4"   # light warm grey
    ]


'''AMP summary table (Figure 1A)'''
df = pd.DataFrame(protein_df)
df = df.round(2)
use_cols = ['Name (abbreviation)', 'Sequence', 'Source', 'Structure', 'Cys-Cys (#)', 'C-terminal NH2', r'MIC ($\mu$M)']

cell_text = df[use_cols].values.tolist()
column_labels = use_cols

fig, ax = plt.subplots(figsize=(150, 60))
ax.set_axis_off()

the_table = ax.table(
    cellText=cell_text,
    colLabels=column_labels,
    cellLoc='center',
    loc='center'
)

the_table.auto_set_font_size(False)
the_table.set_fontsize(200)

num_rows = len(cell_text)
num_cols = len(column_labels)

the_table.auto_set_column_width(list(range(num_cols)))

# Header row: thick black line only below the header
for c in range(num_cols):
    header_cell = the_table.get_celld()[(0, c)]
    header_cell.visible_edges = 'B'  # Only bottom line
    header_cell.set_edgecolor('black')
    header_cell.set_linewidth(2)
    header_cell.set_facecolor('#CCCCCC')  # Distinct background for header
    header_text = header_cell.get_text()
    header_text.set_weight('bold')
    header_text.set_fontname('Arial')
# Data rows: no lines, just alternating background colors
for i in range(num_rows):
    row_color = palette[i % len(palette)]
    for j in range(num_cols):
        data_cell = the_table.get_celld()[(i + 1, j)]
        data_cell.set_linewidth(0)
        data_cell.set_facecolor(row_color)  # Apply row color
        data_text = data_cell.get_text()
        data_text.set_fontname('Arial')

the_table.scale(1, 20)  # Second parameter sets row height
plt.tight_layout()
save_figure('AMP_table')
plt.show()


'''Scatterplot of physicochemical properties (Figure 1B)'''
ft = 38       # axis label font size
ft_2 = 34     # tick label font size
ft_leg = 24   # legend font size
sc = 600      # marker size

# Pairs of columns to plot
pairs = [
    ("Molecular weight (kDa)", "Net charge"),
    ("Lysine fraction",       "Arginine fraction"),
    ("Gravy index",           "Proline fraction"),
    ("Hydrophobic fraction",  "Polar fraction")
]

# Unique hue categories
unique_abbr = protein_df["Abbr."].unique()

# Marker mapping for each group
marker_list = ['o'] * 15
if len(unique_abbr) > len(marker_list):
    raise ValueError("Not enough markers for the number of unique groups.")
marker_dict = dict(zip(unique_abbr, marker_list[:len(unique_abbr)]))

# Distinct colors using the "Paired" palette
n_groups = len(unique_abbr)
colors = sns.color_palette("Paired", n_colors=n_groups)
colors[10] = '#FFDF00'   # brighten one color for visibility
color_dict = dict(zip(unique_abbr, colors))

# Figure with 4 subplots in one row
fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(34, 8))

for idx, (ax, (xx, yy)) in enumerate(zip(axes, pairs)):
    # Plot each group separately (open markers)
    for abbr_, group in protein_df.groupby("Abbr."):
        marker = marker_dict[abbr_]
        color = color_dict[abbr_]
        ax.scatter(
            group[xx],
            group[yy],
            s=500,
            marker=marker,
            facecolors="none",
            edgecolors=color,
            alpha=1,
            linewidths=5,
            label=abbr_
        )

    # Legend only on the last subplot (deduplicated)
    if idx == len(axes) - 1:
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), fontsize=ft_leg)

    # Axis labels and tick parameters
    ax.set_xlabel(xx, fontsize=ft)
    ax.set_ylabel(yy, fontsize=ft)
    ax.tick_params(axis="both", which="major", labelsize=ft_2, length=10)

    # Expand plot limits for whitespace
    x_min, x_max = protein_df[xx].min(), protein_df[xx].max()
    y_min, y_max = protein_df[yy].min(), protein_df[yy].max()
    margin_x = 0.2 * (x_max - x_min)
    margin_y = 0.2 * (y_max - y_min)
    ax.set_xlim(x_min - margin_x, x_max + margin_x)
    ax.set_ylim(y_min - margin_y, y_max + margin_y)

    # Integer y-ticks for "Net charge"
    if yy == "Net charge":
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # Thicker axes
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    ax.tick_params(width=2)

plt.tight_layout(pad=5.0, w_pad=10.0, h_pad=6.0)
save_figure('AMP_scatterplot')
plt.show()


'''Physicochemical-properties table (Figure S1A)'''
feature_cols = [
    'Name (abbreviation)',
    'Molecular weight (kDa)',
    'Net charge',
    'Lysine fraction',
    'Arginine fraction',
    'Gravy index',
    'Proline fraction',
    'Hydrophobic fraction',
    'Polar fraction'
]

df = protein_df[feature_cols].copy().round(2)
cell_text = df.values.tolist()
column_labels = df.columns.tolist()

fig, ax = plt.subplots(figsize=(100, 60))
ax.set_axis_off()

the_table = ax.table(
    cellText=cell_text,
    colLabels=column_labels,
    cellLoc='center',
    loc='center'
)

the_table.auto_set_font_size(False)
the_table.set_fontsize(160)

num_rows = len(cell_text)
num_cols = len(column_labels)

the_table.auto_set_column_width(list(range(num_cols)))

# Header
for c in range(num_cols):
    cell = the_table[(0, c)]
    cell.visible_edges = 'B'
    cell.set_edgecolor('black')
    cell.set_linewidth(2)
    cell.set_facecolor('#CCCCCC')
    cell.get_text().set_weight('bold')
    cell.get_text().set_fontname('Arial')

# Rows
for i in range(num_rows):
    row_color = palette[i % len(palette)]
    for j in range(num_cols):
        cell = the_table[(i + 1, j)]
        cell.set_facecolor(row_color)
        cell.set_linewidth(0)
        cell.get_text().set_fontname('Arial')

the_table.scale(1, 20)
plt.tight_layout()
save_figure('AMP_physicochem_table')
plt.show()


'''Sequences / class-and-type table (Figure S10B)'''
# Class and type per AMP 
class_map = {
    'LL-37':  ('I', 'A'),
    'CecA':   ('I', 'A'),
    'Mag2':   ('I', 'A'),
    'CRAMP':  ('I', 'A'),
    'Tac1':   ('I', 'A'),
    'Mltt':   ('I', 'A'),
    'PG-1':   ('I', 'A'),
    'HBD-3':  ('I', 'B'),
    'Indo':   ('I', 'B'),
    'Tur1B':  ('I', 'B'),
    'PR-39':  ('II', 'B'),
    'Bac7':   ('II', 'B'),
    'Tur1A':  ('II', 'B'),
    'CamBac': ('II', 'B'),
    'DJK-5':  ('II', 'B'),
}

protein_df['Class'] = protein_df['Abbr.'].map(lambda x: class_map.get(x, ('', ''))[0])
protein_df['Type']  = protein_df['Abbr.'].map(lambda x: class_map.get(x, ('', ''))[1])

order = list(class_map.keys())
protein_df['_sort'] = protein_df['Abbr.'].map({k: i for i, k in enumerate(order)})
protein_df = protein_df.sort_values('_sort').drop(columns='_sort')

feature_cols = ['Name (abbreviation)', 'Sequence', 'Class', 'Type']

df = protein_df[feature_cols].copy()
cell_text = df.values.tolist()
column_labels = df.columns.tolist()

fig, ax = plt.subplots(figsize=(120, 60))
ax.set_axis_off()
the_table = ax.table(
    cellText=cell_text,
    colLabels=column_labels,
    cellLoc='center',
    loc='center'
)
the_table.auto_set_font_size(False)
the_table.set_fontsize(160)
num_rows = len(cell_text)
num_cols = len(column_labels)
the_table.auto_set_column_width(list(range(num_cols)))

# Header
for c in range(num_cols):
    cell = the_table[(0, c)]
    cell.visible_edges = 'B'
    cell.set_edgecolor('black')
    cell.set_linewidth(2)
    cell.set_facecolor('#CCCCCC')
    cell.get_text().set_weight('bold')
    cell.get_text().set_fontname('Arial')

# Rows
for i in range(num_rows):
    row_color = palette[i % len(palette)]
    for j in range(num_cols):
        cell = the_table[(i + 1, j)]
        cell.set_facecolor(row_color)
        cell.set_linewidth(0)
        cell.get_text().set_fontname('Arial')

the_table.scale(1, 20)
plt.tight_layout()
save_figure('AMP_class_table')
plt.show()
