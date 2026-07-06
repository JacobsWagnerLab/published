# -*- coding: utf-8 -*-
"""
Quantify LL-37 fluorescence at the bottom of empty trenches.

Uses the single-trench tracked kymographs saved by compute_kymo_fluor_wide-trench.py.
For each selected track the script:
  1. Extracts the mean signal in a bottom-of-trench ROI over time
  2. Measures a lateral background ROI at the same y-position to correct for
     LED intensity fluctuations (signal / background ratio)
  3. Plots background-corrected traces normalised to their pre-injection level
  4. Shows crop images with the ROI box overlaid at selected time points

Requires: quantify_fluor_kymo_functions.py

@author: alessio fragasso
"""

import os
import numpy as np
import matplotlib.pyplot as plt

runfile(r"path_to_quantify_fluor_kymo_functions.py",
        wdir=r"path_to_functions_folder")


'''Settings'''
output_folder_kymo  = #r"path_to_output_folder_kymo"
output_folder_fig   = #r"path_to_figure_output_folder"
master_folder_path  = #r"path_to_strain_master_folder"

AMP_title   = #'AMP condition label for figure titles'
AMP_name    = #'AMP_condition_folder_name'
signal_name = #"LL-37"

channel_name  = "sytox"    # name used when saving kymo outputs
fluor_channel = "fluor2"   # raw image subfolder: "fluor1" or "fluor2"

fr_inj = #0   # frame index at which AMP was added

os.makedirs(output_folder_fig, exist_ok=True)


'''Load single-trench tracks'''
single_tracks_signal = load_single_trench_tracks(
    output_folder_kymo,
    channel_name=channel_name,
)


'''Select pos + track pairs for bottom-signal analysis'''
# List (pos, track_id) pairs of empty trenches to analyse.
pos_track_pairs = [
    ## ("xy01", 0),
    ## ("xy03", 1),
]

# Or alternatively, use track_indices directly (Python list indices):
track_indices_ratio = [
    ## 0,
    ## 1,
]


'''ROI and background parameters'''
bottom_band_px           = 40    # height of the signal ROI [px]
bottom_offset_px         = 0     # offset from detected trench bottom [px]
kymo_x_inset_px          = 50    # inset from trench x-edges for signal ROI [px]
bg_offset_from_center_px = 150   # lateral offset of background ROI from trench centre [px]
bg_width_px              = None  # width of background ROI; None = same as signal ROI
signal_bkg_sub           = False # set True only if GPU bkg-sub was applied to the saved kymo


'''Plot raw mean signal at trench bottom'''
if len(pos_track_pairs) > 0:

    bottom_signal_pdf = os.path.join(
        output_folder_fig,
        f"{AMP_name}_{signal_name}_bottom_trench_signal.pdf",
    )

    bottom_signal_df = plot_bottom_trench_signals_for_pairs(
        single_tracks_signal=single_tracks_signal,
        pos_track_pairs=pos_track_pairs,
        fr_inj=fr_inj,
        bottom_band_px=bottom_band_px,
        bottom_offset_px=bottom_offset_px,
        AMP_title=AMP_title,
        signal_name=signal_name,
        save_path=bottom_signal_pdf,
    )


'''Compute signal / background ratio (corrects LED flickering)'''
if len(track_indices_ratio) > 0:

    bottom_signal_ratio_pdf = os.path.join(
        output_folder_fig,
        f"{AMP_name}_{signal_name}_bottom_trench_signal_bg_ratio.pdf",
    )

    bottom_signal_ratio_df = plot_bottom_trench_signal_ratios_for_track_indices(
        single_tracks_signal=single_tracks_signal,
        track_indices=track_indices_ratio,
        output_folder_kymo=output_folder_kymo,
        master_folder_path=master_folder_path,
        fr_inj=fr_inj,
        fluor_channel=fluor_channel,

        bottom_band_px=bottom_band_px,
        bottom_offset_px=bottom_offset_px,
        kymo_x_inset_px=kymo_x_inset_px,
        bg_offset_from_center_px=bg_offset_from_center_px,
        bg_width_px=bg_width_px,
        signal_bkg_sub=signal_bkg_sub,

        AMP_title=AMP_title,
        signal_name=signal_name,
        save_path=bottom_signal_ratio_pdf,
    )

    '''Background-corrected trace normalised to pre-injection level'''
    corrected_pdf = os.path.join(
        output_folder_fig,
        f"{AMP_name}_{signal_name}_bottom_trench_signal_corrected.pdf",
    )

    bottom_signal_corrected_df = plot_raw_like_corrected_ratio(
        bottom_signal_ratio_df,
        AMP_title=AMP_title,
        signal_name=signal_name,
        ref_time_window=(-np.inf, 0),   # frames before AMP injection
        save_path=corrected_pdf,
    )

    '''Raw vs. background vs. ratio vs. difference (diagnostic)'''
    comparison_pdf = os.path.join(
        output_folder_fig,
        f"{AMP_name}_{signal_name}_bottom_trench_raw_bg_ratio_difference.pdf",
    )

    plot_raw_bg_ratio_difference(
        bottom_signal_ratio_df,
        AMP_title=AMP_title,
        signal_name=signal_name,
        save_path=comparison_pdf,
    )


'''Extract trench crops at selected frames and show ROI box'''
frames_to_extract = [
    fr_inj,          # t = 0
    fr_inj + 30,     # t = +30
    fr_inj + 60,     # t = +60
    fr_inj + 100,    # t = +100
    ## fr_inj + 150,
]

# which track index to show crops for
example_track_index = 0   ## change to desired index

crop_records_example = extract_tracked_trench_crops_for_frames(
    single_tracks_sytox=single_tracks_signal,
    track_index=example_track_index,
    frames_to_extract=frames_to_extract,
    output_folder_kymo=output_folder_kymo,
    master_folder_path=master_folder_path,
    fluor_channel=fluor_channel,
    pad_x=75,
    pad_top=250,
    pad_bottom=100,
    max_frame_delta=2,
    use_sorted=False,
    signal_bkg_sub=signal_bkg_sub,
)

for i, r in enumerate(crop_records_example):
    print(
        f"crop_i: {i} | target_frame: {r['target_frame']} | "
        f"used_frame: {r['frame']} | time_min: {r['frame'] - fr_inj}"
    )


'''Plot individual crops with ROI box'''
for crop_i in range(len(crop_records_example)):

    frame_used = int(crop_records_example[crop_i]["frame"])

    roi_pdf = os.path.join(
        output_folder_fig,
        f"{AMP_name}_{signal_name}_track_{example_track_index:03d}"
        f"_frame_{frame_used:04d}_bottom_roi_crop.pdf",
    )

    plot_bottom_roi_on_crop(
        crop_records=crop_records_example,
        crop_record_index=crop_i,
        AMP_title=AMP_title,
        signal_name=signal_name,
        px_size=0.065841,

        manual_signal_vmin=95,
        manual_signal_vmax=140,
        crop_phase_vmin_prc=5,
        crop_phase_vmax_prc=95,
        phase_invert=False,
        signal_alpha=1,
        zero_background="black",

        overlay_mode="colormap",
        overlay_cmap="inferno",

        bottom_band_px=bottom_band_px,
        bottom_offset_px=bottom_offset_px,
        kymo_x_inset_px=kymo_x_inset_px,

        ylims_px=(45, 900),
        ytick_step_um=10,
        scale_bar_um=5,

        roi_color="cyan",
        roi_linewidth=2,
        save_path=roi_pdf,
        show=True,
    )


'''Crop panels + corrected trace (requires corrected df computed above)'''
crop_record_indices = [0, 1, 2]   # which snapshots to show in the panel row

if len(track_indices_ratio) > 0:

    roi_traces_pdf = os.path.join(
        output_folder_fig,
        f"{AMP_name}_{signal_name}_track_{example_track_index:03d}"
        f"_bottom_ROI_crops_and_corrected_traces.pdf",
    )

    plot_bottom_roi_crops_and_corrected_traces(
        crop_records_example=crop_records_example,
        crop_record_indices=crop_record_indices,
        bottom_signal_corrected_df=bottom_signal_corrected_df,
        fr_inj=fr_inj,
        AMP_title=AMP_title,
        signal_name=signal_name,
        px_size=0.065841,

        bottom_band_px=bottom_band_px,
        bottom_offset_px=bottom_offset_px,
        kymo_x_inset_px=kymo_x_inset_px,

        low_DR=True,
        crop_phase_vmin_prc=5,
        crop_phase_vmax_prc=95,
        phase_invert=False,

        ylims_px=(45, 900),
        ytick_step_um=10,
        scale_bar_um=5,
        figsize=(5, 8),

        roi_color="cyan",
        roi_linewidth=2,
        save_path=roi_traces_pdf,
        show=True,
        xlims=(-10, 120),
    )
