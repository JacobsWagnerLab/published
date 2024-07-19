
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 22:09:45 2023

Code for screener cell types (e.g. swarmer vs stalked) based on cell morphology at birth

@author: A.Fragasso
"""
# import modules
from scipy.signal import lfiltic, lfilter
import skimage
from scipy import io
import os
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import math
import torch
import keyboard
import time
import cupy as cp
import json
import pickle
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import splprep, splev, splrep, interp1d, UnivariateSpline, BSpline
from scipy.signal import lfilter, savgol_filter
from scipy.stats import spearmanr
from skimage import graph, morphology, measure
import imageio.v3 as iio
from scipy.ndimage import gaussian_filter, gaussian_filter1d, binary_fill_holes, binary_closing, binary_opening
from scipy import ndimage
from PIL import Image as im
from skimage.morphology import skeletonize
from skimage.morphology import binary_closing
from skimage.filters import threshold_otsu, threshold_local
from skimage.segmentation import find_boundaries
from skimage import measure, morphology
from numpy.lib.stride_tricks import as_strided
import warnings
warnings.filterwarnings('ignore')

runfile("path_to_analysis_functions_file", wdir="path_to_analysis_folder")

'''Write paths to folders wherer the data are stored'''
experiment_path = "path_to_folder_with_data_to_screen"

'''Create list of positions to scan through'''
pos_list = create_pos_list(1, 10, double_digit = True) + create_pos_list(1, 10, double_digit = False)

save_dict_path = experiment_path + '/output'
save_dict_path_2 = experiment_path + '/output_new'
file_names_list = os.listdir(save_dict_path)
filename_0 = file_names_list[0]
idx = filename_0.find('xy')
exp_name = filename_0[:idx-1]
'''
Input global variables
'''
fr = 3     # [min] time between frame
px_size = 0.064 # px size 
save = True       # if True, save results to df
screen = True     # if True, perform screening
offs = 50         # amount of padding (px) for image display

for pos in pos_list:
    folder_path = experiment_path + '/' +\
        pos  # path to position
    cells_path = folder_path+"/cell"
    mask_path = folder_path + '/masks'
    phase_path = folder_path + '/phase'
    if not os.path.exists(cells_path):
        continue
    cells_list = os.listdir(cells_path)
    if len(cells_list)==0:
        continue
    
    masks_list = os.listdir(mask_path)
    phase_list = os.listdir(phase_path)
    '''
    Import json and pickle files
    '''
    good_cells_list_name = [s for s in file_names_list if 'good_cell' in s and pos in s and exp_name in s]
    if len(good_cells_list_name) == 0:
        continue
    with open(save_dict_path+'/'+good_cells_list_name[0], 'r') as json_file:
        good_cells_list = json.load(json_file)
    with open(save_dict_path+'/'+exp_name+'_'+pos+'_cell_ch.pkl', 'rb') as pickle_file:
        cell_ch_dict = pickle.load(pickle_file)
    with open(save_dict_path+'/'+exp_name+'_'+pos+'_cell_info.pkl', 'rb') as pickle_file:
        cell_info_dict = pickle.load(pickle_file)
    
    files = os.listdir(save_dict_path_2)
    df_name = [s for s in files if 'cell_features_screened_df' in s and pos in s]
    if len(df_name) ==0:
        continue
    
    with open(save_dict_path_2+'/'+df_name[0], 'rb') as f:
        df = pickle.load(f)
    cells_list = list(df['cell_id'].unique())
    df_new = pd.DataFrame()
    for cell in cells_list:
        birth, death, divide, motherID, sisterID, daughterID, gen, pole_1_label, pole_2_label, cell_type, p1_list, p2_list, xaxis_off_list = cell_info_dict[cell]
        cell_df = df[df['cell_id']==cell]
        frame_list = cell_df['frame'].to_list()  
        cell_ch = cell_ch_dict[cell, frame_list[0]]
        cropped_mask_first = cell_ch[0]
        w = cropped_mask_first.shape[1]/100
        h = cropped_mask_first.shape[0]/100
        fig = plt.figure(figsize=(12,9))
        fig.suptitle(cell + ' at position ' + pos + ', ' + exp_name)
        for i in frame_list:
            if i == frame_list[0]:
                phase_temp = iio.imread(phase_path+'/'+phase_list[i])
                mask_temp = iio.imread(mask_path+'/'+masks_list[i])
                cell_ch = cell_ch_dict[cell,i]
                cell_stack_temp = cell_ch[2]
                A, rcm, r_off, angle, x1, x2, y1, y2, edge_flag = get_cell_coordinates(cell_stack_temp)
                y1_n = np.max([y1-offs,0])
                y2_n = np.min([y2+offs,mask_temp.shape[0]])
                x1_n = np.max([x1-offs,0])
                x2_n = np.min([x2+offs,mask_temp.shape[1]])
                
                cell_img = phase_temp[y1_n:y2_n,x1_n:x2_n]
                mask_img = mask_temp[y1_n:y2_n,x1_n:x2_n]
                
                cell_ch = cell_ch_dict[cell, i]
                mask_idx =0
                phase_idx=1
                cell_temp = cell_ch[0]
                cropped_mask = cell_ch[mask_idx]
                cropped_phase = cell_ch[phase_idx]
     
                contours = measure.find_contours(cropped_mask, level=0.5)
                col = 'yellow'

                fig.add_subplot(2,2,2)
                
                plt.title('First frame zoom-out ' + str(i),fontsize=18) 
                plt.imshow(cell_img, cmap='gray')
                plt.scatter(rcm[0]-x1 + offs, rcm[1]-y1 + offs, marker = 'x', color = 'red')
                fig.add_subplot(2,2,1)
                plt.title('First frame ' + str(i),fontsize=18) 
                plt.imshow(cropped_phase, cmap='gray')
                if contours != False:
                    for contour in contours:
                        plt.plot(contour[:, 1], contour[:, 0],
                                 linewidth=1.5, color=col)
                plt.axis('off')  # Remove the axis
            if i == frame_list[-1]:
                phase_temp = iio.imread(phase_path+'/'+phase_list[i])
                mask_temp = iio.imread(mask_path+'/'+masks_list[i])
                cell_ch = cell_ch_dict[cell,i]
                cell_stack_temp = cell_ch[2]
                A, rcm, r_off, angle, x1, x2, y1, y2, edge_flag = get_cell_coordinates(cell_stack_temp)

                y1_n = np.max([y1-offs,0])
                y2_n = np.min([y2+offs,mask_temp.shape[0]])
                x1_n = np.max([x1-offs,0])
                x2_n = np.min([x2+offs,mask_temp.shape[1]])
                
                cell_img = phase_temp[y1_n:y2_n,x1_n:x2_n]
                mask_img = mask_temp[y1_n:y2_n,x1_n:x2_n]
                
                cropped_mask = cell_ch[mask_idx]
                cropped_phase = cell_ch[phase_idx]
                contours = measure.find_contours(cropped_mask, level=0.5)
                col = 'yellow'
                
                fig.add_subplot(2,2,4)

                plt.title('Last frame zoom-out ' + str(i),fontsize=18) 
                plt.imshow(cell_img, cmap='gray')
                plt.scatter(rcm[0]-x1 + offs, rcm[1]-y1 + offs, marker = 'x', color = 'red')
                fig.add_subplot(2,2,3)
                plt.title('Last frame ' + str(i),fontsize=18) 
                plt.imshow(cropped_phase, cmap='gray')
                if contours != False:
                    for contour in contours:
                        plt.plot(contour[:, 1], contour[:, 0],
                                 linewidth=1.5, color=col)
                plt.axis('off')  # Remove the axis
            
 
        plt.tight_layout()
        plt.show()
       
        if screen:
            print(cell+ ': Right if stalked, left if swarmer, up if not sure, down if discard')
            while True:
                if keyboard.is_pressed('right'):
                    cell_df['type'] = 'Stalked' 
                    df_new = pd.concat([df_new,cell_df])
                    print(cell+' is classified as stalked')
                    break
                elif keyboard.is_pressed('left'):
                    cell_df['type'] = 'Swarmer' 
                    df_new = pd.concat([df_new,cell_df])
                    print(cell+' is classified as swarmer')
                    break
                elif keyboard.is_pressed('up'):
                    cell_df['type'] = 'Unsure' 
                    df_new = pd.concat([df_new,cell_df])
                    print(cell+' is classified as unsure')
                    break
                elif keyboard.is_pressed('down'):
                    print(cell+' is discarded')
                    break
                elif keyboard.is_pressed('r'):
                    cell_to_reconsider = input("Insert cell to reconsider: ")
                    print('Reconsidering ' + cell_to_reconsider + '...')
                    df_new = df_new[np.isin(df_new['cell_id'],[cell_to_reconsider] ,invert=True)]
                    cell_df_new = df[df['cell_id']==cell_to_reconsider]
                    
                    type_new = input("Insert corrected classification: ")
                    cell_df_new['type'] = type_new
                    df_new = pd.concat([df_new,cell_df_new])
                    print(cell_to_reconsider +' has been re-classified as ' + type_new)
                    print(cell+ ': Right if stalked, left if swarmer, up if not sure, down if discard')

    if save:
        df_new_name = exp_name+'_'+pos+'_cell_features_screened_df'
        df_new.to_pickle(save_dict_path_2+'/'+df_new_name+'.pkl')
        print('At position ' + pos + ' the number of classified cells is: '+ str(len(df_new['cell_id_pos'].unique())))
       
 