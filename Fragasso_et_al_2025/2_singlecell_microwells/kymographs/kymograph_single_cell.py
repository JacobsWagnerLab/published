# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 16:23:52 2025

@author: fragasso
"""
import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Define custom colormaps
black_to_green = LinearSegmentedColormap.from_list('BlackToGreen', ['black', 'lime'])
black_to_magenta = LinearSegmentedColormap.from_list('BlackToMagenta', ['black', 'magenta'])
black_to_cyan = LinearSegmentedColormap.from_list('BlackToCyan', ['black', 'cyan'])
'''
This function computes and saves mean intensity projections onto medial axis of pixels within px_thr (here  = 4) from the medial axis
'''
def compute_1D_mean_int_proj_single_cell(medial_proj_df, save =False, return_ = False, experiment_path='', exp_name='',
                             px_thr = 4, aug = 20, bin_size_l=0.065):
    medial_proj_df = medial_proj_df[medial_proj_df['px_dist']<px_thr*aug]
    cell_list = medial_proj_df['cell_id'].unique().tolist()
    medial_df_all= medial_proj_df.groupby(
        ['cell_id', 'frame', ]).apply(process_group).reset_index(drop=True)

    mean_proj_all_df = pd.DataFrame()
    for cell in cell_list:
        print('Computing the mean 1D intensity proj for cell ' + cell)
        medial_proj_df_cell = medial_df_all[medial_df_all['cell_id'] == cell]
        frame_list = medial_proj_df_cell['frame'].unique().tolist()
        rep = medial_proj_df_cell['rep'].iloc[0]
        for i in frame_list:
            medial_proj_frame_df = medial_proj_df_cell[medial_proj_df_cell['frame'] == i]
            length = medial_proj_frame_df['length_um'].iloc[0]
            n_bins_l = length / bin_size_l
            l_range = (max(medial_proj_frame_df['proj_abs']) + max(medial_proj_frame_df['proj_abs']) * 0.001) - (-0.0001)
            step_l = l_range / n_bins_l
            bins_l = np.arange(start=-0.00001, stop=max(medial_proj_frame_df['proj_abs']) + max(medial_proj_frame_df['proj_abs']) * 0.1, step=step_l)
            bins_l_idx = np.digitize(medial_proj_frame_df['proj_abs'], bins_l, right=True)
            medial_proj_frame_df['bins_l_idx'] = bins_l_idx
            medial_df_average = medial_proj_frame_df.groupby('bins_l_idx').apply(calculate_averages_1D).reset_index(drop=True)  
            medial_df_average['cell_id'] = cell
            medial_df_average['frame'] = i
            medial_df_average['cell_id'] = cell
            medial_df_average['exp'] = medial_proj_frame_df['exp'].iloc[0]
            medial_df_average['rep'] = medial_proj_frame_df['rep'].iloc[0]
            mean_proj_all_df = pd.concat([mean_proj_all_df,medial_df_average])
    if save:
        if save:
            save_folder = experiment_path + '/output_1D_mean_proj'          # was: experiment_path + '/'+rep+'/output_1D_mean_proj'
            if not os.path.exists(save_folder):
                os.makedirs(save_folder)
            mean_axial_1D_df_path = save_folder + '/' + 'mean_1D_axial_df_' + exp_name+'_'+rep+'_px_'+str(px_thr)+'.pkl'
            with open(mean_axial_1D_df_path, 'wb') as pickle_file:
                pickle.dump(mean_proj_all_df, pickle_file)
            print('Dataframe of 1D axial mean projections were saved in ' + save_folder)    
    if return_:
        return mean_proj_all_df
    

'''Helper functions'''

"""
Calculate Euclidean distances between consecutive coordinates
"""

def calculate_euclidean_distance(coords, pixel_size, aug, medial_axis_=[]):
    # Calculate the difference in x and y coordinates between consecutive pixels
    if len(medial_axis_)>0:   # if medial axis is provided, add start coordinate to get the correct absolute coordinates of the mean proj xx vector
        coords = np.vstack((medial_axis_[0],coords))
    dx = np.diff(coords[:, 0]) / aug
    dy = np.diff(coords[:, 1]) / aug

    # Calculate the Euclidean distance between consecutive pixels
    distances = pixel_size * np.sqrt(dx**2 + dy**2) 
    absolute_distances = np.insert(np.cumsum(distances), 0, 0)
    return distances[1:], absolute_distances[1:]

def process_group(df, px_size = 0.065841, aug = 20):
    # Calculate medial axis and projections
    m_proj_coord = np.array([np.array([s, t])
                            for (s, t) in df['m_proj_coord']])
    _, m_proj_abs = calculate_euclidean_distance(
        m_proj_coord, px_size, aug)
    m_proj_abs = np.insert(m_proj_abs, 0, 0)

    # Pixel distance and conversions
    px_dist = df['px_dist']
    px_dist_um = px_dist * px_size / 20

    # Assign new columns directly to the DataFrame
    df['px_dist_abs'] = px_dist_um
    df['proj_abs'] = m_proj_abs
    df['proj_abs_norm'] = m_proj_abs / max(m_proj_abs)
    df['px_dist_abs_norm'] = px_dist_um/max(px_dist_um)
    df['length_um'] = m_proj_abs[-1]
    df['dist_um'] = np.max(px_dist_um)
    return df

def calculate_averages_1D(df_group):
    return pd.Series({
        'px_int_phase_mean': np.mean(df_group['px_int_phase']),
        'px_int_fluor1_mean': np.mean(df_group['px_int_mScarlet']),
        'px_int_fluor2_mean': np.mean(df_group['px_int_bfp']),
        'px_int_fluor3_mean': np.mean(df_group['px_int_sytox_green']),
        'proj_abs_mean': np.mean(df_group['proj_abs']),
        'bins_l_idx': df_group['bins_l_idx'].iloc[0]
    })

def find_x1_x2_columns(arr):
    x1 = arr.shape[1]  
    x2 = 0           
    for row in range(arr.shape[0]):
        non_zero_indices = np.where(arr[row, :] > 0)[0]
        if len(non_zero_indices) > 0:
            x1 = min(x1, non_zero_indices[0])
            x2 = max(x2, non_zero_indices[-1])
    return x1, x2

'''
Imports the 1D mean-projection dataframes and builds/plots/saves a kymograph per cell,
showing raw mScarlet and BFP only.
'''
def build_export_kymographs_single_cell(experiment_path, exp_name, ch_info_list=['mScarlet', 'bfp', 'sytox_green'],
                                        ft=18, bin_size_l=0.065, fr=1, save=True, show_=False):
    output_path = experiment_path + '/output_1D_mean_proj'
    if not os.path.exists(output_path):
        print('No output_1D_mean_proj in ' + experiment_path)
        return
    df_name = [s for s in os.listdir(output_path) if 'mean_1D_axial_df' in s]
    if len(df_name) == 0:
        print('No mean_1D_axial_df in ' + output_path)
        return
    df_name = df_name[0]
    with open(output_path + '/' + df_name, 'rb') as f:
        mean_proj_all_df = pickle.load(f)
    cell_list = mean_proj_all_df['cell_id'].unique().tolist()
    for cell in cell_list:
        print(cell)
        cell_df = mean_proj_all_df[mean_proj_all_df['cell_id'] == cell]
        rep = cell_df['rep'].iloc[0]
        exp_name = cell_df['exp'].iloc[0]
        t = np.sort(cell_df['frame'].unique())
        x_int = np.arange(t[0], t[-1] + 1, step=fr) - np.min(t)
        y = max(cell_df.index) + 5

        kymo_fluor_dict = {'mScarlet': np.zeros([y, len(x_int)]),
                           'bfp': np.zeros([y, len(x_int)])}
        x_min = np.min(t)
        cell_contours = np.zeros([2, len(x_int)])
        for i in t:
            frame_df = cell_df[cell_df['frame'] == i]
            diff = y - (max(frame_df.index) + 1)
            v1 = int(diff / 2)
            v2 = diff - v1
            k = int(i - x_min)
            cell_contours[0, k] = v1
            cell_contours[1, k] = y - v2
            kymo_fluor_dict['mScarlet'][v1:y - v2, k] = frame_df['px_int_fluor1_mean']
            kymo_fluor_dict['bfp'][v1:y - v2, k] = frame_df['px_int_fluor2_mean']

        fig = plt.figure(figsize=(18, 6))
        fig.suptitle(cell, fontsize=18, y=0.98)

        def add_colorbar(im, ax, label, fontsize):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="3%", pad=0.1)
            cbar = fig.colorbar(im, cax=cax)
            cbar.ax.tick_params(labelsize=fontsize)
            cbar.set_label(label, fontsize=fontsize)

        panels = [('mScarlet', black_to_magenta, 'mScarlet'),
                  ('bfp', black_to_cyan, 'BFP')]
        for p, (key, cmap, title) in enumerate(panels):
            ax = fig.add_subplot(1, 2, p + 1)
            vmin = np.percentile(kymo_fluor_dict[key], 1)
            vmax = np.percentile(kymo_fluor_dict[key], 99)
            im = ax.imshow(kymo_fluor_dict[key], cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax))
            ax.plot(cell_contours[0, :], linewidth=1.5, linestyle='--', c='white')
            ax.plot(cell_contours[1, :], linewidth=1.5, linestyle='--', c='white')
            add_colorbar(im, ax, 'Mean Intensity', ft)
            ax.set_title(title, fontsize=20)
            ax.set_xlabel('Time (min)', fontsize=20, fontweight='bold')
            ax.set_ylabel('Length (um)', fontsize=20, fontweight='bold')

            ylim, xlim = ax.get_ylim(), ax.get_xlim()
            yt = ax.get_yticks()
            xt = ax.get_xticks()
            ax.set_yticks(yt)
            ax.set_yticklabels(np.round(yt * bin_size_l, 1), fontsize=ft)
            ax.set_xticks(xt)
            ax.set_xticklabels([int(s + x_min) for s in xt], fontsize=ft)
            ax.set_ylim(ylim)
            ax.set_xlim(xlim)

        plt.tight_layout(rect=[0, 0, 1, 0.95])

        if save:
            print('Exporting kymo for', cell)
            save_folder = os.path.join(experiment_path, 'output_single-cell_kymo')
            os.makedirs(save_folder, exist_ok=True)
            plt.savefig(os.path.join(save_folder, f'{cell}.png'), format='png', bbox_inches='tight')
        if show_:
            plt.show()
        else:
            plt.close()