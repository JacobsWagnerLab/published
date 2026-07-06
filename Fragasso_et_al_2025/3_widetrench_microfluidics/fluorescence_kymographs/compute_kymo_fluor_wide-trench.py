# -*- coding: utf-8 -*-
"""
Compute and export fluorescence kymographs for wide-trench microscopy data.

Steps:
1) Detect wide-trench positions from the phase contrast signal using Omnipose masks:
   - Background-subtract phase and fluor channels
   - Use the closed/filled Omnipose mask to restrict phase to the cell-occupied region
   - Call get_trench_xycoords() to locate individual trench x-coordinates from the
     averaged phase intensity profile; extrapolate missing trenches from the periodic pitch
2) Extract kymographs for each channel (fluor/SYTOX and phase):
   - Per-frame: extract expanded crops around each trench, average into ensemble 2D frame
     and 1D y-profile (kymograph row)
   - Per-trench: link crops across frames by x_center to build single-trench kymographs
3) Save all outputs (arrays, dataframes, metadata) to output_folder_kymo/
4) Load saved outputs, pool across positions, and plot:
   - Ensemble overlay kymograph: phase (greyscale) + fluor (green), time-aligned to AMP injection
   - Single-trench side-by-side: raw crop image overlaid with fluorescence + individual kymograph

@author: alessio fragasso
Tue Jun 20 22:09:45 2023
"""

import sys
import skimage
from scipy import io
import os
import re
import gc
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import glob
import math
import torch
import time
import cupy as cp
from scipy.signal import find_peaks, savgol_filter
from skimage import graph, morphology, measure
import imageio.v3 as iio
import imageio
from scipy.ndimage import gaussian_filter, gaussian_filter1d, binary_fill_holes, binary_closing, binary_opening, center_of_mass
from scipy import ndimage
from scipy import ndimage as ndi_np
from skimage.morphology import skeletonize, disk
from skimage.filters import threshold_otsu, threshold_local
from skimage import measure, morphology
import warnings
warnings.filterwarnings('ignore')

# Check if GPU is available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Device:", device)
print("CUDA Available:", torch.cuda.is_available())
print("Number of GPUs:", torch.cuda.device_count())

# decide which GPU to use
cp.cuda.Device(0).use()

runfile(r"path_to_analysis_functions_library.py",
        wdir=r"path_to_folder")

runfile(r"path_to_quantify_sytox_functions.py",
        wdir=r"path_to_folder")


'''Global variables'''
fr = 1    # [min] time between frame
px_size = 0.065841  # um/px   (pixel size, 100x)
offs_x = 10
offs_y = 10
save = True
bkg_sub = True
show = False
overwrite = True

disk50 = disk(50)

pos_list = create_pos_list(1, 20, double_digit=True) + create_pos_list(1, 20, double_digit=False)
ch_list = ['phase', 'fluor1', 'fluor2']
ch_fluor_list = [s for s in ch_list if 'fluor1' in s]    # only fluor1 (SYTOX/permeabilization channel)


'''Specify paths'''
master_folder_path = #r"path_to_strain_master_folder"

'''Insert list of AMP (conditions) to process'''
AMP_list = [
    ## 'AMP_condition_1',
    ## 'AMP_condition_2',
    ]

'''Insert list of lists of frames at which AMP was injected. One list per AMP condition, one value per replicate'''
fr_inj_master_list = [
    ## [fr_inj_rep_1, fr_inj_rep_2],   # for AMP_condition_1
    ## [fr_inj_rep_1, fr_inj_rep_2],   # for AMP_condition_2
    ]

'''AMP-specific trench detection parameters.
peak_p: prominence threshold for phase peak detection.
search_aft: pixel column after which to search for the trench bottom edge.
fixed_gap: expected gap width between trenches in pixels.
Tune these per condition by running with show=True on a few frames.'''
trench_params = {
    ## 'AMP_condition_1': {'peak_p': 800, 'search_aft': 1150, 'fixed_gap': 180},
    ## 'AMP_condition_2': {'peak_p': 650, 'search_aft': 1000, 'fixed_gap': 179},
}


'''Processing loop: trench detection and kymograph extraction'''
for l in range(len(AMP_list)):
    AMP_name = AMP_list[l]
    fr_inj_list = fr_inj_master_list[l]
    params = trench_params.get(AMP_name, {'peak_p': 800, 'search_aft': 1150, 'fixed_gap': 180})

    rep_list = [s for s in os.listdir(master_folder_path+'/'+AMP_name) if '.DS_Store' not in s]

    for s in range(len(rep_list)):
        rep_name = rep_list[s]
        exp_name = AMP_name + '_' + rep_name + '_' + os.path.basename(master_folder_path)
        experiment_path = master_folder_path+'/'+AMP_name+'/'+rep_name
        print(exp_name)
        fr_inj = fr_inj_list[s] if len(fr_inj_list) > 0 else 0

        for pos in pos_list:
            '''
            Get paths to all folders
            '''
            folder_path = experiment_path + '/' + pos
            cells_path = folder_path + "/cell"
            mask_path = folder_path + '/masks'

            pos_num = re.findall(r'\d+', pos)[0]
            if not os.path.exists(cells_path):
                continue
            elif len(os.listdir(cells_path)) == 0:
                continue
            print(pos)
            if os.path.exists(mask_path):
                masks_list = os.listdir(mask_path)

            paths_dict = {}
            file_list_dict = {}
            for ch in ch_list:
                paths_dict[ch] = folder_path + '/' + ch
                file_list_dict[ch] = os.listdir(paths_dict[ch])

            phase_list = file_list_dict['phase']
            phase_path = paths_dict['phase']

            if not overwrite:
                output_path = os.path.join(experiment_path, 'output_folder_kymo')
                if os.path.exists(output_path):
                    out_list = os.listdir(output_path)
                    check_name = pos + '_' + exp_name + '_sytox_kymo_arrays.npz'
                    if check_name in out_list:
                        print(check_name + ' was already computed -->> skip')
                        continue
                    else:
                        print(check_name + ' not computed already')

            for ch in ch_fluor_list:
                fluor_list = file_list_dict[ch]
                fluor_path = paths_dict[ch]

                pad = 300
                bin_size = 5

                sytox_profiles, sytox_avg_2d_frames = [], []
                phase_profiles, phase_avg_2d_frames = [], []
                frame_ids = []
                trench_rel_box = None
                trench_box_rows = []
                frame_summary_rows = []
                single_trench_records_sytox = []
                single_trench_records_phase = []

                for k in range(len(fluor_list)):
                    mask_temp = iio.imread(mask_path + '/' + masks_list[k])
                    phase_temp = iio.imread(phase_path + '/' + phase_list[k])
                    fluor_temp = iio.imread(fluor_path + '/' + fluor_list[k])

                    # Background-subtract fluor1 (SYTOX)
                    bs = int(fluor_temp.shape[0] / 16)
                    fluor_bkg_sub, _, _ = local_bkg_sub_cp(
                        fluor_temp, mask_temp, pos, exp_name, frame=str(k),
                        box_size=bs, dilations=15, sigma_=60, show=False)

                    # Background-subtract phase; restrict to Omnipose mask region
                    bs = int(fluor_temp.shape[0] / 64)
                    phase_temp_bkg_sub, mask_temp, _ = local_bkg_sub_cp(
                        phase_temp, mask_temp, pos, exp_name, frame=str(k),
                        box_size=bs, dilations=1, sigma_=25, show=False)

                    mask_dil = ndi_np.binary_closing(mask_temp > 0, structure=disk50, brute_force=True)
                    mask_dil = ndi_np.binary_fill_holes(mask_dil)
                    phase_temp_bkg_sub *= mask_dil

                    show_trench = (k % 20 == 0)
                    if show_trench:
                        print(k)
                        plt.imshow(phase_temp_bkg_sub, cmap='gray_r')
                        plt.title(AMP_name + ', frame ' + str(k))
                        plt.show()

                    title_plot = AMP_name + ', ' + pos + ', frame' + str(k)
                    trench_boxes, phase_avg_x, peaks, ensemble_y_profile = get_trench_xycoords(
                        phase_temp,
                        peak_distance=60,
                        peak_prominence=params['peak_p'],
                        trench_width_range=(70, 130),
                        gap_width_range=(150, 250),
                        edge_search_px=100,
                        min_edge_drop=100,
                        fixed_trench_width_px=125,
                        fixed_gap_px=params['fixed_gap'],
                        refined_width_tolerance_px=45,
                        extrapolate_missing=True,
                        trench_height_px=535,
                        search_after=params['search_aft'],
                        show=show_trench,
                        frame=k,
                        title=title_plot,
                    )

                    sytox_result = extract_expanded_trench_ensemble(
                        fluor_bkg_sub, trench_boxes, pad=pad,
                        kymo_use_trench_x_only=True, kymo_x_inset_px=50)

                    phase_result = extract_expanded_trench_ensemble(
                        phase_temp_bkg_sub, trench_boxes, pad=pad,
                        kymo_use_trench_x_only=True, kymo_x_inset_px=50)

                    if sytox_result is None or phase_result is None:
                        continue

                    sytox_y_profile, sytox_frame_avg_2d, sytox_expanded_crops, used_trench_boxes, trench_rel_box = sytox_result
                    phase_y_profile, phase_frame_avg_2d, phase_expanded_crops, _, _ = phase_result

                    for crop_i, (tr_id, y0, y1, x0, x1) in enumerate(used_trench_boxes):
                        sytox_crop = sytox_expanded_crops[crop_i]
                        phase_crop = phase_expanded_crops[crop_i]

                        rec_base = {
                            'frame': k, 'trench_id_raw': tr_id,
                            'x0': x0, 'x1': x1, 'y0': y0, 'y1': y1,
                            'x_center': 0.5 * (x0 + x1),
                        }
                        rec_sytox = rec_base.copy(); rec_sytox['crop'] = sytox_crop
                        rec_phase = rec_base.copy(); rec_phase['crop'] = phase_crop
                        single_trench_records_sytox.append(rec_sytox)
                        single_trench_records_phase.append(rec_phase)

                        trench_box_rows.append({
                            'exp_name': exp_name, 'AMP': AMP_name, 'rep': rep_name,
                            'pos': pos, 'frame': k, 'trench_id': tr_id,
                            'y0': y0, 'y1': y1, 'x0': x0, 'x1': x1, 'pad': pad,
                            'rel_y0': trench_rel_box['y0'], 'rel_y1': trench_rel_box['y1'],
                            'rel_x0': trench_rel_box['x0'], 'rel_x1': trench_rel_box['x1'],
                        })

                    sytox_profiles.append(sytox_y_profile)
                    sytox_avg_2d_frames.append(sytox_frame_avg_2d)
                    phase_profiles.append(phase_y_profile)
                    phase_avg_2d_frames.append(phase_frame_avg_2d)
                    frame_ids.append(k)

                    frame_summary_rows.append({
                        'exp_name': exp_name, 'AMP': AMP_name, 'rep': rep_name,
                        'pos': pos, 'frame': k,
                        'n_trenches_used': len(sytox_expanded_crops),
                        'mean_sytox': float(np.nanmean(sytox_frame_avg_2d)),
                        'max_sytox': float(np.nanmax(sytox_frame_avg_2d)),
                        'median_sytox': float(np.nanmedian(sytox_frame_avg_2d)),
                        'mean_phase': float(np.nanmean(phase_frame_avg_2d)),
                        'max_phase': float(np.nanmax(phase_frame_avg_2d)),
                        'median_phase': float(np.nanmedian(phase_frame_avg_2d)),
                    })

                single_trench_records_sytox = link_trench_records(single_trench_records_sytox, max_dx=40)
                single_trench_records_phase = link_trench_records(single_trench_records_phase, max_dx=40)

                track_outputs_sytox = build_single_trench_kymos(
                    single_trench_records_sytox, trench_rel_box=trench_rel_box,
                    min_frames=20, kymo_use_trench_x_only=True, kymo_x_inset_px=50)

                track_outputs_phase = build_single_trench_kymos(
                    single_trench_records_phase, trench_rel_box=trench_rel_box,
                    min_frames=20, kymo_use_trench_x_only=True, kymo_x_inset_px=50)

                if save and ch == 'fluor1' and len(sytox_profiles) > 0:
                    output_folder_kymo = os.path.join(experiment_path, 'output_folder_kymo')

                    save_channel_kymo_outputs(
                        output_folder_kymo=output_folder_kymo, pos=pos, exp_name=exp_name,
                        AMP_name=AMP_name, rep_name=rep_name, channel_name='sytox',
                        frame_ids=frame_ids, profiles=sytox_profiles,
                        avg_2d_frames=sytox_avg_2d_frames,
                        trench_rel_box=trench_rel_box, bin_size=bin_size)

                    save_channel_kymo_outputs(
                        output_folder_kymo=output_folder_kymo, pos=pos, exp_name=exp_name,
                        AMP_name=AMP_name, rep_name=rep_name, channel_name='phase',
                        frame_ids=frame_ids, profiles=phase_profiles,
                        avg_2d_frames=phase_avg_2d_frames,
                        trench_rel_box=trench_rel_box, bin_size=bin_size)

                    save_trench_metadata_outputs(
                        output_folder_kymo=output_folder_kymo, pos=pos,
                        exp_name=exp_name, trench_box_rows=trench_box_rows,
                        frame_summary_rows=frame_summary_rows)

                    save_single_trench_track_outputs(
                        output_folder_kymo=output_folder_kymo, pos=pos, exp_name=exp_name,
                        AMP_name=AMP_name, rep_name=rep_name, channel_name='sytox',
                        track_outputs=track_outputs_sytox, trench_rel_box=trench_rel_box)

                    save_single_trench_track_outputs(
                        output_folder_kymo=output_folder_kymo, pos=pos, exp_name=exp_name,
                        AMP_name=AMP_name, rep_name=rep_name, channel_name='phase',
                        track_outputs=track_outputs_phase, trench_rel_box=trench_rel_box)
