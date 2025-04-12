# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 15:38:02 2025

@author: fragasso

Code for importing, cropping, and exporting wide-trenches videos (MM=mother machine).
"""

import numpy as np
import os
import re
import gc
import pandas as pd
import imageio as iio
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from matplotlib import ticker
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import mode

'''This function prompts overlay of first and last frame of the phase contrast image sequence, 
and allows user to define ROI to crop and plot a single wide trench over time'''
def get_user_roi(image1, image2, title = None):
    """Overlay first and last image with a grid and prompt user for ROI coordinates, allowing corrections."""
    while True:
        overlay = np.zeros((*image1.shape, 3), dtype=np.uint8)
        overlay[..., 1] = (image1 / image1.max() * 255).astype(np.uint8)  # Green channel
        overlay[..., 2] = (image2 / image2.max() * 255).astype(np.uint8)  # Magenta channel
        
        plt.figure(figsize=(12, 8))
        plt.imshow(overlay)
        plt.title(title if title else 'Enter x,y coordinates')
        
        # Add grid and increase tick density
        plt.grid(color='white', linestyle='-', linewidth=0.5)
        plt.xticks(np.arange(0, image1.shape[1], step=50), rotation=45)  # Increase tick density
        plt.yticks(np.arange(0, image1.shape[0], step=50), rotation=45)
        
        plt.show()
        
        x_min = int(input("Enter x_min: "))
        x_max = int(input("Enter x_max: "))
        y_min = int(input("Enter y_min: "))
        y_max = int(input("Enter y_max: "))
        
        # Replot with selected ROI
        fig = plt.figure(figsize=(12, 12))
        ax1 = fig.add_subplot(2,1,1)
        plt.imshow(image1)
        # plt.title(title if title else 'Chosen region of interest')
        plt.grid(color='white', linestyle='--', linewidth=0.5)
        plt.xticks(np.arange(0, image1.shape[1], step=50), rotation=45)
        plt.yticks(np.arange(0, image1.shape[0], step=50), rotation=45)
        plt.gca().add_patch(plt.Rectangle((x_min, y_min), x_max-x_min, y_max-y_min, edgecolor='red', facecolor='none', linewidth=2))
        
        ax1 = fig.add_subplot(2,1,2)
        plt.imshow(image2)
        # plt.title(title if title else 'Chosen region of interest')
        plt.grid(color='white', linestyle='--', linewidth=0.5)
        plt.xticks(np.arange(0, image1.shape[1], step=50), rotation=45)
        plt.yticks(np.arange(0, image1.shape[0], step=50), rotation=45)
        plt.gca().add_patch(plt.Rectangle((x_min, y_min), x_max-x_min, y_max-y_min, edgecolor='red', facecolor='none', linewidth=2))
        
        plt.show()
        
        confirm = input("Are the coordinates correct? (y/n): ").strip().lower()
        if confirm == 'y':
            break
        else:
            print("Re-enter ROI coordinates.")
    
    return y_min, y_max, x_min, x_max



'''This function fills top and bottom pixels added for padding upon drift correction from SuperSegger, with a color calculated from the top and bottom
pixel intensities.

Steps:
   1) Identify 'padding' by global mode (pad_val).
   2) Find bounding box of non-padding: [y_min, y_max].
   3) If fill_val_top is None:
        - Compute fill_val_top from top 'top_rows' lines inside bounding box (ignore leftover pad pixels).
      If fill_val_bottom is None:
        - Compute fill_val_bottom from bottom 'bottom_rows' lines inside bounding box (ignore leftover pad pixels).
   4) Overwrite rows [0 : y_min + overlap) with fill_val_top, clamped to image height if needed.
   5) Overwrite rows [y_max - overlap + 1 : end) with fill_val_bottom, ensuring we don't go above y_min.
'''
def fill_top_bottom_constant_overlap(
    image,
    fill_val_top=None,
    fill_val_bottom=None,
    overlap=3,
    top_rows=3,
    bottom_rows=3
):
   
    filled = image.astype(np.float32).copy()
    h, w = filled.shape

    # 1) Identify global 'padding' value
    pad_val = float(mode(image, axis=None).mode)

    # 2) Find bounding box of non-padding (any pixel != pad_val)
    mask_nonpad = (image != pad_val)
    rows_nonpad = np.any(mask_nonpad, axis=1)

    # If there's no real content, fill everything with top or bottom fill (or pad_val) 
    # and return
    if not np.any(rows_nonpad):
        if fill_val_top is None:
            fill_val_top = pad_val
        if fill_val_bottom is None:
            fill_val_bottom = pad_val
        filled[:, :] = fill_val_top  # or split top/bottom as needed
        return filled, fill_val_top, fill_val_bottom

    y_idx = np.where(rows_nonpad)[0]
    y_min, y_max = y_idx[0], y_idx[-1]

    # Helper function: median from bounding-box rows [start : end+1], ignoring pad pixels
    def compute_median_in_rows(start, end):
        region = filled[start:end+1, :]
        real_vals = region[region != pad_val]
        if len(real_vals) == 0:
            return pad_val
        return float(np.median(real_vals))

    # 3) Compute fill_val_top if None
    if fill_val_top is None:
        # Use up to 'top_rows' lines from [y_min : y_min+top_rows]
        top_end = min(y_max, y_min + top_rows - 1)
        fill_val_top = compute_median_in_rows(y_min, top_end)

    # Compute fill_val_bottom if None
    if fill_val_bottom is None:
        # Use up to 'bottom_rows' lines from [y_max - bottom_rows + 1 : y_max]
        bottom_start = max(y_min, y_max - bottom_rows + 1)
        fill_val_bottom = compute_median_in_rows(bottom_start, y_max)

    # 4) Overwrite the TOP region [0 : y_min + overlap)
    top_fill_end = min(y_min + overlap, h)  # don't exceed image bottom
    filled[0:top_fill_end, :] = fill_val_top
    top_pad = False
    bottom_pad = False
    if y_min>0:
        top_pad = True
    if y_max<h:
        bottom_pad = True

    # 5) Overwrite the BOTTOM region [y_max - overlap + 1 : h)
    bottom_fill_start = max(0, y_max - overlap + 1)  # don't go negative
    filled[bottom_fill_start:h, :] = fill_val_bottom

    return filled, fill_val_top, fill_val_bottom, top_pad, bottom_pad


'''
Build colormaps that gradually span from black to a given color or palette
'''
original_viridis = plt.cm.viridis
original_inferno = plt.cm.inferno
original_winter = plt.cm.winter

# Modified colormaps with black at the bottom
viridis_black = create_black_bottom_colormap(original_viridis, black_ratio=0.5)
inferno_black = create_black_bottom_colormap(original_inferno, black_ratio=0.3)
winter_black = create_black_bottom_colormap(original_winter, black_ratio=0.1)

# Define a new colormap going from black to green
black_to_green = LinearSegmentedColormap.from_list('BlackToGreen', ['black', 'lime'])
black_to_magenta = LinearSegmentedColormap.from_list('BlackToGreen', ['black', 'magenta'])
black_to_cyan = LinearSegmentedColormap.from_list('BlackToGreen', ['black', 'cyan'])

'''
Input global variables
'''
fr = 1    # [min] time between frame
px_size = 0.065841  # um/px   (pixel size in alpha scope, 100x)
save = True
show = False
concat = False
crop_ch = True

max_rpla = 1550
min_rpla = 120
max_hu = 1250
min_hu = 220
max_phase = 25000
min_phase = 2500

scalebar_length_um = 5  # Desired scale bar length in micrometers
scalebar_length_px = scalebar_length_um / px_size  # Length in pixels

pos_list = create_pos_list(1, 200, digits = 2) 
ch_list = ['phase', 'fluor1', 'fluor2']
ch_fluor_list = [s for s in ch_list if 'fluor' in s]
ch_fluor_col = ['grey', 'viridis', 'inferno']

'''Specify paths'''
master_folder_path = r"N:\Alessio_Fragasso\Microscopy_AMP\9_MM\1_CJW7753"   # from mothership2
save_video_path = r"I:\Science communication\Manuscripts\Papers\2024\AMP story\Videos\output"

AMP_label = ## insert name of the AMP/condition##

t_off = 0    ## insert time offset if videos across consecutive subfolders need to be concatenated
t_max = 26900   ## max time of the experiment (can crop it to specified time)

pos_list = ['xy01']   # insert xy position to analyse
crop_num = 1          # number the single cropped trench

ch_names = ['Phase','RplA-GFP','HupA-mCherry']    # inster channel names, they will appear on top of the image sequences


AMP_list = [
    # List of AMP/condition folders, containin the xy positions folder from SuperSegger to analyse
    ]

fr_inj_master_list = [
    # List of lists of AMP addition frames for each replicate for each AMP/condition
]


wide_tr = True          # Adjusts distances for printed labels in the image sequence
# wide_tr = False

count=0


for l in range(len(AMP_list)):
    AMP_name = AMP_list[l]
    fr_inj_list = fr_inj_master_list[l]
    rep_list = [s for s in os.listdir(master_folder_path+'/'+AMP_name) if '.DS_Store' not in s]    
    for s in range(len(rep_list)):
        rep_name = rep_list[s]
        rep_name = rep_list[1]
        exp_name = AMP_name + '_' + rep_name + '_' + os.path.basename(master_folder_path)
        experiment_path = master_folder_path+'/'+AMP_name+'/'+rep_name
        print(exp_name)
        if concat and s>0:
            t_off = len(phase_list)
        if len(fr_inj_list)==0:
            fr_inj = 0
        else:
            fr_inj = fr_inj_list[s]
        for pos in pos_list:
            '''
            Get paths to all folders
            '''
            
            folder_path = experiment_path + '/' +\
                pos  # path to position
            if not os.path.exists(folder_path):
                continue
            
            pos_num = re.findall(r'\d+', pos)[0]
            print(pos)
            paths_dict = {}
            file_list_dict = {}
            
            for ch in ch_list:   
                paths_dict[ch] = folder_path + '/'+ ch
                file_list_dict[ch] = os.listdir(paths_dict[ch]) 
                    
            '''
            Initialize dictionaries and lists
            '''
            cell_stack_dict = {}  # dictionary of cells parameters over time
            cell_ch_dict = {}  # dictionary of cells phase and mask over time
            cell_info_dict = {}   # dictionary of cell info: birth frame (starts from 1), death frame, divide or not, motherID, sisterID, daughterID
            cell_features_all_df = pd.DataFrame()   
            
            phase_list =file_list_dict['phase']
            phase_path = paths_dict['phase']
            
            fluor1_list = file_list_dict['fluor1']
            fluor1_path = paths_dict['fluor1']
            
            fluor2_list = file_list_dict['fluor2']
            fluor2_path = paths_dict['fluor2']
            
            phase_image_paths = [os.path.join(phase_path, f) for f in phase_list]
            fluor1_image_paths = [os.path.join(fluor1_path, f) for f in fluor1_list]
            fluor2_image_paths = [os.path.join(fluor2_path, f) for f in fluor2_list]
            
            phase_img_list =[iio.imread(img) for img in phase_image_paths[:t_max]]
            fluor1_img_list =[iio.imread(img) for img in fluor1_image_paths[:t_max]]
            fluor2_img_list =[iio.imread(img) for img in fluor2_image_paths[:t_max]]
            
            if crop_ch:
                y_min, y_max, x_min, x_max = get_user_roi(phase_img_list[0], phase_img_list[-1], title = exp_name + ', ' + pos)
                phase_cropped = [phase_img[y_min:y_max, x_min:x_max] for phase_img in phase_img_list]
                fluor1_cropped = [fluor1_img[y_min:y_max, x_min:x_max] for fluor1_img in fluor1_img_list]
                fluor2_cropped = [fluor2_img[y_min:y_max, x_min:x_max] for fluor2_img in fluor2_img_list]
            else:
                y_min, y_max, x_min, x_max
                phase_cropped = corrected_phase
                fluor1_cropped = corrected_fluor1
                fluor2_cropped = corrected_fluor2
            
            
            max_fluor_dict =  {'phase':0,'fluor1': 0,'fluor2': 0,'fluor3': 0}
            min_fluor_dict =  {'phase':10000,'fluor1': 10000,'fluor2': 10000}
            param_fluor_dict = {'phase':[99,5],'fluor1': [99,5],'fluor2': [99,10]}
            
            
            for i in range(len(phase_cropped)):
                if i%50==0: 
                     print(i)
                     
                     
                max_fluor_dict['phase'] = np.max([max_fluor_dict['phase'],np.percentile(phase_cropped[i], param_fluor_dict['phase'][0])])
                min_fluor_dict['phase'] = np.min([min_fluor_dict['phase'],np.percentile(phase_cropped[i], param_fluor_dict['phase'][1])])
                
                max_fluor_dict['fluor1'] = np.max([max_fluor_dict['fluor1'],np.percentile(fluor1_cropped[i], param_fluor_dict['fluor1'][0])])
                min_fluor_dict['fluor1'] = np.min([min_fluor_dict['fluor1'],np.percentile(fluor1_cropped[i], param_fluor_dict['fluor1'][1])])
                
                max_fluor_dict['fluor2'] = np.max([max_fluor_dict['fluor2'],np.percentile(fluor2_cropped[i], param_fluor_dict['fluor2'][0])])
                min_fluor_dict['fluor2'] = np.min([min_fluor_dict['fluor2'],np.percentile(fluor2_cropped[i], param_fluor_dict['fluor2'][1])])
         
            ft = 14
            # Before loop, define variables to store fill_values
            fill_val_top_phase = None
            fill_val_bottom_phase = None
            
            fill_val_top_fluor1 = None
            fill_val_bottom_fluor1 = None
            
            fill_val_top_fluor2 = None
            fill_val_bottom_fluor2 = None
            
            top_pad = False
            bottom_pad = False
            top_found = False
            bottom_found = False
            
            # while not top_pad:
            for k in range(len(phase_img_list)):
                
                # Extract images for this frame
                phase_temp = phase_cropped[k]
                fluor1_temp = fluor1_cropped[k]
                fluor2_temp = fluor2_cropped[k]
                
                _, fill_val_top_phase_temp, fill_val_bottom_phase, top_pad, bottom_pad = fill_top_bottom_constant_overlap(
                    phase_temp, 
                    fill_val_top=None, 
                    fill_val_bottom=None,
                    overlap=3,     # Overwrite 3 real rows in each direction
                    top_rows=3,    # Use top 3 bounding-box rows for fill_val_top
                    bottom_rows=3  # Use bottom 3 bounding-box rows for fill_val_bottom
                )
        
                _, fill_val_top_fluor1_temp, fill_val_bottom_fluor1, top_pad, bottom_pad = fill_top_bottom_constant_overlap(
                    fluor1_temp,
                    fill_val_top=None,
                    fill_val_bottom=None,
                    overlap=3,
                    top_rows=3,
                    bottom_rows=3
                )
        
                _, fill_val_top_fluor2_temp, fill_val_bottom_fluor2, top_pad, bottom_pad = fill_top_bottom_constant_overlap(
                    fluor2_temp,
                    fill_val_top=None,
                    fill_val_bottom=None,
                    overlap=3,
                    top_rows=3,
                    bottom_rows=3
                )
                # print(top_pad)
                
                if top_pad and not top_found:
                    top_found=True
                    fill_val_top_phase = fill_val_top_phase_temp
                    fill_val_top_fluor1 = fill_val_top_fluor1_temp
                    fill_val_top_fluor2 = fill_val_top_fluor2_temp
                    print('found top at frame '+ str(k))
                if bottom_pad and not bottom_found:
                    bottom_found=True
                    fill_val_bottom_phase = fill_val_top_phase_temp
                    fill_val_bottom_fluor1 = fill_val_top_fluor1_temp
                    fill_val_bottom_fluor2 = fill_val_top_fluor2_temp
                    print('found bottom at frame '+ str(k))
                if top_found and bottom_found:
                    print('Found both!')
                    break

            for k in range(len(phase_img_list)):
                
                # Extract images for this frame
                phase_temp = phase_cropped[k]
                fluor1_temp = fluor1_cropped[k]
                fluor2_temp = fluor2_cropped[k]
                
           
                # On SUBSEQUENT frames, reuse same fill values
                phase_filled, _, _, _, _ = fill_top_bottom_constant_overlap(
                    phase_temp,
                    fill_val_top_phase,
                    fill_val_bottom_phase,
                    overlap=3,
                    top_rows=3,
                    bottom_rows=3
                )
        
                fluor1_filled, _, _, _, _ = fill_top_bottom_constant_overlap(
                    fluor1_temp,
                    fill_val_top_fluor1,
                    fill_val_bottom_fluor1,
                    overlap=3,
                    top_rows=3,
                    bottom_rows=3
                )
        
                fluor2_filled, _, _, _, _ = fill_top_bottom_constant_overlap(
                    fluor2_temp,
                    fill_val_top_fluor2,
                    fill_val_bottom_fluor2,
                    overlap=3,
                    top_rows=3,
                    bottom_rows=3
                )
                 
                image_height = phase_temp.shape[0]
                image_width = phase_temp.shape[1]
                
                bar_x = image_width - scalebar_length_px - 5  # 10 pixels from right edge
                bar_y = image_height - 10  # 20 pixels from bottom edge
                
                current_time = convert_minutes_to_hr_min((k+t_off)*fr)
               
                if wide_tr:
                    fig = plt.figure(figsize=(10, 8))
                else: fig = plt.figure(figsize=(7, 8))
                sz_bar = '15%'
                gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1], wspace=0.01)  # Reduce space between subplots
                ax1 = fig.add_subplot(gs[0])
                vmin = min_phase
                vmax = max_phase
                im1 = ax1.imshow(phase_filled, cmap='grey', vmin=vmin, vmax=vmax)
                ax1.set_title(ch_names[0],fontsize = 14)
                
                ax1.add_patch(patches.Rectangle((bar_x, bar_y), scalebar_length_px, 2, 
                                                edgecolor='white', facecolor='white', linewidth=3))
                ax1.text(10, 10, current_time, color='black', fontsize=14, fontweight='bold', ha='left', va='top')
                if k >= fr_inj:
                    if wide_tr:
                        ax1.text(10, 40, '+ ' + AMP_label, color='red', fontsize=14, fontweight='bold', ha='left', va='top')
                    else: 
                        ax1.text(10, 30, '+ ' + AMP_label, color='red', fontsize=14, fontweight='bold', ha='left', va='top')

                divider1 = make_axes_locatable(ax1)
                cax1 = divider1.append_axes("right", size=sz_bar, pad=0.1)
                cbar1 = plt.colorbar(im1, cax=cax1)
                cbar1.ax.tick_params(labelsize=10)
                cbar1.set_label('Signal intensity (a.u.)',fontsize = 10)
                cbar1.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{x/1000:.1f}'))
                # Add the text "x10³" on top of the colorbar
                cbar1.ax.text(0.5, 1.01, 'x10$^3$', transform=cbar1.ax.transAxes,
                               ha='center', va='bottom')
                ax1.axis('off')
                
                ax2 = fig.add_subplot(gs[1])
                vmin = min_rpla
                vmax = max_rpla
                im2 = ax2.imshow(fluor1_filled, cmap=black_to_green, vmin=vmin, vmax=max_rpla, rasterized=True)
                ax2.set_title(ch_names[1],fontsize = 14)
                
                divider2 = make_axes_locatable(ax2)
                cax2 = divider2.append_axes("right", size=sz_bar, pad=0.1)
                cbar2 = plt.colorbar(im2, cax=cax2)
                cbar2.ax.tick_params(labelsize=10)
                cbar2.set_label('Signal intensity (a.u.)',fontsize = 10)
                cbar2.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{x/1000:.1f}'))
                # Add the text "x10³" on top of the colorbar
                cbar2.ax.text(0.5, 1.01, 'x10$^3$', transform=cbar2.ax.transAxes,
                               ha='center', va='bottom')
                ax2.axis('off')
                
                # ---- Subplot 3: Fluor2 ----
                ax3 = fig.add_subplot(gs[2])

                vmin = min_hu
                vmax = max_hu
                im3 = ax3.imshow(fluor2_filled, cmap='inferno', vmin=vmin, vmax=max_hu, rasterized=True)
                ax3.set_title(ch_names[2],fontsize = 14)
                
                divider3 = make_axes_locatable(ax3)
                cax3 = divider3.append_axes("right", size=sz_bar, pad=0.1)
                cbar3 = plt.colorbar(im3, cax=cax3)
                cbar3.ax.tick_params(labelsize=10)
                cbar3.set_label('Signal intensity (a.u.)',fontsize = 10)
                cbar3.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{x/1000:.1f}'))
                cbar3.ax.text(0.5, 1.01, 'x10$^3$', transform=cbar3.ax.transAxes,
                               ha='center', va='bottom')
                ax3.axis('off')
                
                # Reduce overall whitespace
                plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.05, wspace=0.10)  # Adjust width spacing
                

                
                if save:
                    save_path = save_video_path  + '/MM/' + AMP_label +'/'+ pos+'_cropped_' + str(crop_num)
                    if not os.path.exists(save_path):
                        os.makedirs(save_path)
                    if crop_ch:
                        frame_path = save_path + '/' + 'Frame_'+str(k+t_off) + '_'+exp_name
                    plt.savefig(frame_path + '.png', format='png', bbox_inches="tight")
                
                plt.show()
                
                # Clear local variables to free memory
                del phase_temp, fluor1_temp, fluor2_temp, fig
                gc.collect()
                                

            