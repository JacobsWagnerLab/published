
# -*- coding: utf-8 -*-
"""
Post-processing and curation script for wide-trench single-cell data.

Differences from post_processing_curation_functions.py (standard PDMS-microwell single-cell data):
- Cells are tracked in 2D; rcm_x_offset is computed per position with optional axis flip (if original trenches were imaged upside down)
- Spatial filtering via xwin_by_pos dictionary to restrict analysis to valid trench regions (exclude partially imaged trenches, for example, or that were lost due to x-drift)
- clean_track_drop_one curation step removes frames where the cell center jumps > threshold pixels

Output:
- Curated cell-feature dataframe exported to output_2/ as '_cell_features_df_trenches'

@author: alessio fragasso
Fri Jan 16 10:59:18 2026
"""

from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import numpy as np
import os
from scipy import io
import skimage
import warnings
import re
import imageio.v3 as iio
from matplotlib.colors import Normalize
import matplotlib.cm as cm
warnings.filterwarnings('ignore')
from matplotlib.collections import LineCollection
from matplotlib import cm, colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

runfile(r"path_to_analysis_functions_library.py",
        wdir=r"path_to_folder")

runfile(r"path_to_post_processing_functions.py",
        wdir=r"path_to_folder")

runfile(r"path_to_plot_videos_functions.py",
        wdir=r"path_to_folder")


'''Global variables'''
fr = 1    # [min] time between frame
px_size = 0.065841  # um/px   (pixel size, 100x)


'''Import dataframes'''
master_folder_path = #r"path_to_strain_master_folder"

AMP_name = #'AMP_condition_folder_name'

exp_name = AMP_name + '_' + os.path.basename(master_folder_path)

rep_list = [s for s in os.listdir(master_folder_path+'/'+AMP_name) if '.DS_Store' not in s]

rep_path_list = [master_folder_path+'/'+AMP_name+'/'+rep for rep in rep_list]
exp_list = [AMP_name+'_'+rep+'_'+os.path.basename(master_folder_path) for rep in rep_list]
cell_features_all_df, cell_id_list = import_cell_df(
    rep_path_list, exp_list, fr, df_mark='cell_features_df', output_folder='output')

df = cell_features_all_df.copy()


'''Compute rcm_x_offset: align cell positions within each xy position.
flip_list controls orientation: 1 = already well oriented, 0 = needs to be flipped.'''
flip_list = np.ones(len(df['xy'].unique()))    # set per-position flip flags as needed

pos_list = df['xy'].unique()
for pos, do_flip in zip(pos_list, flip_list):
    m = df['xy'] == pos
    y = df.loc[m, 'rcm_x'].to_numpy()
    if do_flip:
        df.loc[m, 'rcm_x_offset'] = y - y.max()
    else:
        df.loc[m, 'rcm_x_offset'] = y.min() - y
    print(y.max(), ', xy='+pos)


'''Analyse and select trenches: define valid x-range per position.
Restrict each position to the trench region of interest.'''
xwin_by_pos = {
    ## 'pos_01': (x_min, x_max),
    ## 'pos_02': (x_min, x_max),
    }

df_tr_list = []
for pos in sorted(df['xy'].unique()):
    if len(xwin_by_pos)==0:
        x1, x2 = -2000, 4200    # default: keep all; replace with xwin_by_pos[pos] to restrict
    else: x1, x2 = xwin_by_pos[pos] 
    df_pos = df[df['xy'] == pos]
    df_pos_new = df_pos[(df_pos['rcm_y'] > x1) & (df_pos['rcm_y'] < x2)]
    df_tr_list.append(df_pos_new)

    df_tr = pd.concat(df_tr_list, ignore_index=True)
    sns.scatterplot(df_pos_new, x='rcm_y', y='rcm_x_offset', hue='xy')
    plt.show()

df_tr = df_tr.groupby('cell_id', group_keys=False).apply(
    add_gr_columns(window_size=15, min_len_for_savgol=15))

df = df_tr.copy()

df = (df.groupby('cell_id', group_keys=False)
        .apply(clean_track_drop_one, thr=50.0, max_iter=3))    # remove frames where rcm jumps > 50 px


'''Curation steps'''
df = curate_short_trajectories(df, 15)
df = curate_fast_jumps(df, gr_max=0.04, gr_min=-0.01, gr='norm_gr_savgol')
df = curate_area(df, max_area=3000, min_area=200)


'''Plot cell area and growth rate single-cell trajectories'''
fig = plt.figure(figsize=(25, 10))
fig.suptitle(exp_name+', N='+str(len(df['cell_id'].unique().tolist())))

x_var = 'time_min'

ax1 = fig.add_subplot(1, 3, 1)
sns.lineplot(data=df, x=x_var, y='area_smooth_savgol', estimator=None,
             lw=0.8, units='cell_id', alpha=0.3)
plt.xlabel('Time (min)', fontsize=14, fontweight='bold')
plt.ylabel('Area (px)', fontsize=14, fontweight='bold')

ax2 = fig.add_subplot(1, 3, 2)
sns.lineplot(df, x=x_var, y='norm_gr_savgol', estimator=None,
             lw=0.8, alpha=0.3, units='cell_id')
plt.grid(True)
plt.xlabel('Time (min)', fontsize=14, fontweight='bold')
plt.ylabel('Normalized growth rate (1/min)', fontsize=14, fontweight='bold')
plt.ylim(-0.05, 0.05)

ax3 = fig.add_subplot(1, 3, 3)
sns.lineplot(df, x=x_var, y='rcm_x_offset', estimator=None,
             lw=0.8, alpha=0.3, units='cell_id')
plt.xlabel('Time (min)', fontsize=14, fontweight='bold')
plt.ylabel('rcm y-coord', fontsize=14, fontweight='bold')
ax3.invert_yaxis()

plt.show()


'''Export curated dataframe'''
df = df.groupby('cell_id', group_keys=False).apply(add_gr_columns)    # recalculate gr after curation
export_df(df, rep_path_list, exp_list,
          output_folder='output_2', df_name='_cell_features_df_trenches')
