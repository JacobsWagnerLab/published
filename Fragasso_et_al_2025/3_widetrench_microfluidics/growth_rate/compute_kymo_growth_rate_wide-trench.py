# -*- coding: utf-8 -*-
"""
Compute and export spatial kymographs of normalized growth rate for wide-trench data.

Steps:
1) Import curated cell-feature dataframes for one or more AMP conditions
2) Time-align traces to the frame of AMP injection (time_min_offs = 0)
3) Bin cells spatially along the trench axis (rcm_x_offset_neg) in dx-pixel bins
4) Optionally prune isolated kymograph pixels using a neighborhood-count filter
5) Plot and export the kymograph as a PDF

@author: alessio fragasso
Wed Sep 27 12:22:30 2023
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
import warnings
import re
import imageio.v3 as iio
from matplotlib.colors import Normalize, LinearSegmentedColormap
import matplotlib.cm as cm
warnings.filterwarnings('ignore')
from matplotlib.collections import LineCollection
from matplotlib import cm, colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pathlib import Path
from scipy import ndimage as ndi

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

AMP_list = [
    ## 'AMP_condition_1',
    ## 'AMP_condition_2',
    ## 'AMP_condition_3',
    ]

df_all = pd.DataFrame()
for AMP_name in AMP_list:
    exp_name = AMP_name + '_' + os.path.basename(master_folder_path)
    rep_list = [s for s in os.listdir(master_folder_path+'/'+AMP_name) if '.DS_Store' not in s]
    rep_path_list = [master_folder_path+'/'+AMP_name+'/'+rep for rep in rep_list]
    exp_list = [AMP_name+'_'+rep+'_'+os.path.basename(master_folder_path) for rep in rep_list]
    cell_features_all_df, cell_id_list = import_cell_df(
        rep_path_list, exp_list, fr, df_mark='df_tr_filtered', output_folder='output_tr_filtered')
    cell_features_all_df['AMP'] = AMP_name
    df_all = pd.concat([df_all, cell_features_all_df])


'''Plot kymographs'''
prune = True    # enable neighborhood-based pruning of isolated kymograph pixels

base_dir = Path(#r"path_to_figure_output_folder"
                )
plot_dir = base_dir / "kymo_plots"
df_dir = base_dir / "kymo_dfs"
plot_dir.mkdir(parents=True, exist_ok=True)
df_dir.mkdir(parents=True, exist_ok=True)

for AMP in AMP_list:
    df_tr = df_all[df_all['AMP'] == AMP].copy()

    AMP_label = AMP.split('_')[1] + ', ' + AMP.split('_')[2]

    df_tr['time_min_offs'] = df_tr['time_min'] - fr * df_tr['fr_inj']
    df_tr['rcm_x_offset_neg'] = -df_tr['rcm_x_offset']

    '''Spatial binning along trench axis'''
    dx = 12    # bin width [pixels]; adjust as needed
    rcm_min_global = df_tr['rcm_x_offset_neg'].min()
    rcm_max_global = df_tr['rcm_x_offset_neg'].max()

    start = np.floor(rcm_min_global / dx) * dx
    stop = np.ceil(rcm_max_global / dx) * dx
    rcm_edges = np.arange(start, stop + dx, dx)

    df_tr['rcm_bin'] = pd.cut(df_tr['rcm_x_offset_neg'], bins=rcm_edges, include_lowest=True, right=False)
    rcm_centers = (rcm_edges[:-1] + rcm_edges[1:]) / 2
    df_tr['rcm_center'] = rcm_centers[df_tr['rcm_bin'].cat.codes]

    df_kymo_raw = (
        df_tr
        .groupby(['time_min_offs', 'rcm_bin'], observed=True)['norm_gr_savgol']
        .median()
        .reset_index()
    )

    Z_raw = (
        df_kymo_raw
        .pivot(index='time_min_offs', columns='rcm_bin', values='norm_gr_savgol')
        .sort_index(axis=0)
        .sort_index(axis=1)
    )

    occupied = Z_raw.notna().to_numpy()

    '''Prune isolated pixels: remove kymograph pixels with fewer than min_neighbors
    occupied neighbors in a (2*dt+1) x (2*dy+1) window, iterated n_iter times.
    The center column (same y-bin) is excluded from the kernel so isolated single-trench
    streaks do not self-validate across time.'''
    if prune:
        dt = 1             # half-window in time bins
        dy = 2             # half-window in spatial bins
        kernel = np.ones((2*dt + 1, 2*dy + 1), dtype=int)
        kernel[:, dy] = 0  # exclude same y-bin from neighbor count
        min_neighbors = 5
        n_iter = 5

        keep_pixels = occupied.copy()
        for _ in range(n_iter):
            neighbor_count = ndi.convolve(keep_pixels.astype(int), kernel,
                                          mode='constant', cval=0)
            new_keep_pixels = keep_pixels & (neighbor_count >= min_neighbors)
            if np.array_equal(new_keep_pixels, keep_pixels):
                break
            keep_pixels = new_keep_pixels

        Z = Z_raw.where(keep_pixels)

        keep_rc, keep_cc = np.where(keep_pixels)
        keep_pairs = pd.MultiIndex.from_arrays(
            [Z_raw.index.to_numpy()[keep_rc], Z_raw.columns.to_numpy()[keep_cc]],
            names=['time_min_offs', 'rcm_bin'])
        df_pairs = pd.MultiIndex.from_frame(df_tr[['time_min_offs', 'rcm_bin']])
        df_tr_surv = df_tr[df_pairs.isin(keep_pairs)].copy()
    else:
        Z = Z_raw.copy()
        df_tr_surv = df_tr.copy()

    t = Z.index.to_numpy()
    Z_mat = Z.to_numpy()

    N_cells_final = df_tr_surv['cell_id'].nunique()
    N_trenches = int(np.nansum(df_tr_surv['N_trenches'].unique()))

    '''Plot kymograph'''
    vmin, vmax = -0.01, 0.025
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.get_cmap('berlin')

    y0_um = 0.0
    y1_um = (rcm_edges[-1] - rcm_edges[0]) * px_size

    fig, ax = plt.subplots(figsize=(6, 4))
    fig.patch.set_facecolor('white')
    ax.set_facecolor((0.85, 0.85, 0.85))

    im = ax.imshow(
        Z_mat.T,
        origin='lower',
        aspect='auto',
        cmap=cmap,
        norm=norm,
        extent=[t.min(), t.max(), y0_um, y1_um],
        interpolation='nearest'
    )

    ax.set_xlabel('Time since start of AMP injection (min)', fontsize=16)
    ax.set_ylabel(r'Center of mass y-coordinate ($\mu$m)', fontsize=16)
    ax.set_title(f'{AMP_label}, N_cells = {N_cells_final}, N_trenches = {N_trenches}')
    cb = fig.colorbar(im, ax=ax)
    cb.set_label('Normalized growth rate (1/min)', fontsize=16)
    cb.ax.yaxis.set_tick_params(labelsize=16)
    ax.axvline(x=0, linestyle=':', c='black')
    ax.set_ylim(-5, 50)
    plt.tight_layout()

    safe_label = AMP_label.replace(',', '').replace(' ', '_').replace('/', '-')
    plt.savefig(plot_dir / f'{safe_label}.pdf', format='pdf')
    plt.show()
