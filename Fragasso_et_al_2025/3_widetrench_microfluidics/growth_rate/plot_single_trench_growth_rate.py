
# -*- coding: utf-8 -*-
"""
Generate frame-by-frame videos of a single wide-trench position showing:
  panel 1 — phase contrast image
  panel 2 — ribosome fluorescence (fluor1)
  panel 3 — nucleoid fluorescence (fluor2)
  panel 4 — scatter of instantaneous growth rate at each cell centre
  panel 5 — single-cell growth-rate trajectories vs. time (kymograph-style)

Frames are saved as PDF or PNG and can be assembled into a video externally.

@author: alessio fragasso
Wed Sep 27 12:22:30 2023
"""

from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
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
from skimage import measure
from skimage.measure import regionprops

runfile(r"path_to_analysis_functions_library.py",
        wdir=r"path_to_folder")

runfile(r"path_to_post_processing_functions.py",
        wdir=r"path_to_folder")

runfile(r"path_to_plot_videos_functions.py",
        wdir=r"path_to_folder")


'''Custom colormaps'''
black_to_green = LinearSegmentedColormap.from_list('BlackToGreen', ['black', 'lime'])
black_to_magenta = LinearSegmentedColormap.from_list('BlackToMagenta', ['black', 'magenta'])
black_to_cyan = LinearSegmentedColormap.from_list('BlackToCyan', ['black', 'cyan'])

original_viridis = plt.cm.viridis
original_inferno = plt.cm.inferno
original_winter = plt.cm.winter
viridis_black = create_black_bottom_colormap(original_viridis, black_ratio=0.1)
inferno_black = create_black_bottom_colormap(original_inferno, black_ratio=0.1)
winter_black = create_black_bottom_colormap(original_winter, black_ratio=0.1)


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
    rep_path_list, exp_list, fr, df_mark='cell_features_df_trenches', output_folder='output_2')

df = cell_features_all_df.copy()


'''Select position and load images'''
pos = #'xy_position_number'    # e.g. '09'

df_pos = df[(df['xy'] == pos) & (df['rep'] == rep_list[0])]

phase_path = master_folder_path+'/'+AMP_name+'/'+df['rep'].iloc[0]+'/xy'+pos+'/phase'
phase_list = [iio.imread(os.path.join(phase_path, f))
              for f in sorted(os.listdir(phase_path)) if f.endswith('.tif')]

masks_path = master_folder_path+'/'+AMP_name+'/'+df['rep'].iloc[0]+'/xy'+pos+'/masks'
masks_list = [iio.imread(os.path.join(masks_path, f))
              for f in sorted(os.listdir(masks_path)) if f.endswith('.png')]

ribo_path = master_folder_path+'/'+AMP_name+'/'+df['rep'].iloc[0]+'/xy'+pos+'/fluor1'
ribo_list = [iio.imread(os.path.join(ribo_path, f))
             for f in sorted(os.listdir(ribo_path)) if f.endswith('.tif')]

hu_path = master_folder_path+'/'+AMP_name+'/'+df['rep'].iloc[0]+'/xy'+pos+'/fluor2'
hu_list = [iio.imread(os.path.join(hu_path, f))
           for f in sorted(os.listdir(hu_path)) if f.endswith('.tif')]


'''Video settings'''
vmin, vmax = -0.01, 0.025
norm = colors.Normalize(vmin=vmin, vmax=vmax)
cmap = 'berlin'
save_video = True
export_pdf = True
flip = True    # True if cells are oriented top-to-bottom in the image

# Crop coordinates [pixels]: x1/x2 = left/right column; y1/y2 = bottom/top row margin
x1, x2, y1, y2 = 0, 0, 0, 0    # <-- set crop window for your position

AMP_label = #'AMP_label_for_figure_text'

save_path = os.path.join(
    #r"path_to_video_output_folder",
    AMP_name, rep_list[0], 'xy'+pos)
if export_pdf:
    save_path = save_path + '_pdf'


'''Generate video frames'''
for t in range(int(np.max(df_pos['time_min']))):
    labels = masks_list[int(t)]
    ribos = ribo_list[int(t)]
    hus = hu_list[int(t)]
    phase = phase_list[int(t)]

    H, W = labels.shape
    labels_crop = labels[y2:H-y1, x1:x2]
    ribos_crop = ribos[y2:H-y1, x1:x2]
    hus_crop = hus[y2:H-y1, x1:x2]
    phase_crop = phase[y2:H-y1, x1:x2]

    df_crop = df_pos[
        (df_pos['rcm_y'] >= x1) & (df_pos['rcm_y'] < x2) &
        (df_pos['rcm_x'] >= y2) & (df_pos['rcm_x'] < (H - y1))
    ].copy()
    df_crop['rcm_y_crop'] = df_crop['rcm_y'] - x1
    df_crop['rcm_x_crop'] = df_crop['rcm_x'] - y2

    # Filter out small or dim regions from the mask
    labels_filtered = labels_crop.copy()
    for region in regionprops(labels_crop, intensity_image=ribos_crop):
        if region.area < 200 or region.mean_intensity < 500:
            labels_filtered[labels_filtered == region.label] = 0
    labels_crop = labels_filtered.copy()

    # Build contour list and per-label colormap
    contours_all = []
    for label_id in np.unique(labels_crop):
        if label_id == 0:
            continue
        single_mask = (labels_crop == label_id)
        for contour in measure.find_contours(single_mask, level=0.5):
            contours_all.append((label_id, contour))

    idx, uniq = indexed_image_from_labels(labels_crop, background=0)
    cmap_ = cmap_for_labels(uniq, bg_color=(0, 0, 0, 0))
    color_for = {int(lab): cmap_(i+1) for i, lab in enumerate(uniq)}

    max_label = int(labels_crop.max())
    colors_list = [color_for.get(i, 'black') for i in range(max_label + 1)]
    cmap_synced = mpl.colors.ListedColormap(colors_list)
    cmap_synced.colors[0] = (0, 0, 0, 1)

    H_crop, W_crop = labels_crop.shape
    current_time = convert_minutes_to_hr_min(t)

    scalebar_length_um = 5
    scalebar_length_px = scalebar_length_um / px_size
    bar_x = W_crop - scalebar_length_px - 5
    bar_y = H_crop - 10

    df0 = df_crop[df_crop['time_min'] == t]
    if len(df0) > 0:
        fr_inj = df0['fr_inj'].iloc[0]
        cell_list_curr = df0['cell_id'].unique()
        df0_curr = df_crop[df_crop['time_min'] <= t]
    else:
        fr_inj = df_crop['fr_inj'].iloc[0]
        cell_list_curr = []
        df0_curr = df_crop[df_crop['time_min'] <= t]

    fig, axs = plt.subplots(1, 5, figsize=(17, 7.5), constrained_layout=True,
                            gridspec_kw={'width_ratios': [0.8, 0.8, 0.8, 0.8, 2]})

    # Panel 1: phase contrast
    ax1 = axs[0]
    ax1.imshow(phase_crop, cmap='grey')
    ax1.add_patch(patches.Rectangle((bar_x, bar_y), scalebar_length_px, 2,
                                    edgecolor='white', facecolor='white', linewidth=3))
    ax1.text(10, 10, current_time, color='white', fontsize=16, fontweight='bold', ha='left', va='top')
    if t >= fr_inj:
        ax1.text(10, 40, '+ '+AMP_label, color='yellow', fontsize=16, fontweight='bold', ha='left', va='top')
    ax1.set_xticks([]); ax1.set_yticks([])

    # Panel 2: ribosome fluorescence
    axs[1].imshow(ribos_crop, cmap=black_to_green, vmin=200, vmax=1400)
    axs[1].set_xticks([]); axs[1].set_yticks([])

    # Panel 3: nucleoid fluorescence
    axs[2].imshow(hus_crop, cmap=inferno_black, vmin=200, vmax=900)
    axs[2].set_xticks([]); axs[2].set_yticks([])

    # Panel 4: instantaneous growth rate scatter
    ax3 = axs[3]
    ax3.set_facecolor((0.9, 0.9, 0.9))
    if len(df0) > 0:
        ax3.scatter(data=df0, x='rcm_y', y='rcm_x',
                    c='norm_gr_savgol', cmap=cmap, norm=norm, s=35, alpha=1)
    ax3.set_xlim(x1, x2)
    ax3.set_ylim(y2, H-y1)
    if flip:
        ax3.invert_yaxis()
    ax3.set_aspect('equal')
    ax3.set_xticks([]); ax3.set_yticks([])

    # Panel 5: growth-rate trajectories vs. time
    ax4 = axs[4]
    cmap_obj = cm.get_cmap(cmap) if isinstance(cmap, str) else cmap
    for cid, g in df0_curr.groupby('cell_id'):
        g = g.sort_values('time_min')
        x = g['time_min'].to_numpy()
        y = (H_crop - 1 - g['rcm_x_crop'].to_numpy()) * px_size
        c = g['norm_gr_savgol'].to_numpy()
        points = np.column_stack([x, y]).reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=cmap, norm=norm, array=c[:-1],
                            linewidth=1.4, alpha=0.7)
        ax4.add_collection(lc)
        ax4.plot(x[0], y[0], marker='.', markersize=lc.get_linewidths()[0],
                 linestyle='None', color=cmap_obj(norm(c[0])), alpha=0.7, zorder=3)

    ax4.set_xlim(-5, df_crop['time_min'].max()+5)
    ax4.set_facecolor((0.9, 0.9, 0.9))
    ax4.set_ylim(0, (H-y1-y2)*px_size)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax4, location='right', shrink=1, aspect=21, pad=0.07)
    cb.ax.yaxis.set_tick_params(labelsize=16)
    cb.set_label('Normalized growth rate (1/min)', fontsize=16)
    ax4.tick_params(axis='both', labelsize=16)
    ax4.set_xlabel('Time (min)', fontsize=16)
    ax4.set_ylabel(r'Center of mass y-coordinate ($\mu$m)', fontsize=16)

    if save_video:
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        if export_pdf:
            plt.savefig(save_path+f'/frame_{t}.pdf', format='pdf')
        else:
            plt.savefig(save_path+f'/frame_{t}.png', format='png', dpi=300)

    plt.show()
