# -*- coding: utf-8 -*-
"""
Library of functions for C.Crescentus analysis

@author: A.Fragasso
"""

import os
import pandas as pd
import numpy as np
import pickle
from pims import ND2Reader_SDK
import seaborn as sns
import matplotlib.pyplot as plt
import skimage
import math 
import cupy as cp
import cupyx.scipy.ndimage as ndi
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap
from cupyx.profiler import benchmark
from scipy.spatial.distance import pdist, squareform
from numpy.lib.stride_tricks import as_strided
from skimage.morphology import remove_small_objects, binary_erosion
from skimage.transform import resize
from skimage import filters, morphology, measure
from skimage.registration import phase_cross_correlation
from skimage.measure import label, regionprops
from skimage.feature import peak_local_max
from skimage.segmentation import watershed, expand_labels
from skimage.segmentation import relabel_sequential
from skimage.filters import threshold_otsu, threshold_local
from scipy import ndimage
from scipy.ndimage import zoom, label, center_of_mass, sum as ndi_sum
from scipy.signal import find_peaks, lfiltic, lfilter
from scipy.optimize import least_squares, curve_fit
from scipy.ndimage import gaussian_filter, gaussian_filter1d, binary_fill_holes, binary_opening, binary_dilation, binary_erosion
from scipy.interpolate import splprep, splev
from scipy.signal import lfilter, savgol_filter
import imageio.v3 as iio
import scipy.misc 
from PIL import Image as im

def format_number_as_string(number):
    formatted_number = f"cell{number:07d}"
    return formatted_number

def create_pos_list(start, end, double_digit = True):
    pos_list = []
    if double_digit:
        pos_list = [f'{i:02d}' for i in range(start, end + 1)]
    else:
        pos_list =  [str(i) for i in range(start, end + 1)]
    pos_list = ['xy'+pos for pos in pos_list]
    return pos_list

def generate_random_cmap(num_colors):
    np.random.seed(42)
    colors_list = []
    while len(colors_list) < num_colors:
        color = np.random.rand(3)
        if np.linalg.norm(color - np.ones(3)) > 0.3:
            colors_list.append(color)
    colors_list = np.vstack(([1, 1, 1], colors_list))  # Keep 0 as black
    return colors.ListedColormap(colors_list)

def curate_short_trajectories(cell_features_df, area_inc = 150, time_int = 10):
    cell_list = cell_features_df['cell_id_pos'].unique().tolist()
    cell_df_new = pd.DataFrame()
    for cell in cell_list:
        cell_df = cell_features_df[cell_features_df['cell_id_pos'] == cell]
        if  (cell_df['area'].iloc[-1] - cell_df['area'].iloc[0])>area_inc and (cell_df['time_min'].iloc[-1] - cell_df['time_min'].iloc[0])>time_int:
            cell_df_new = pd.concat([cell_df_new,cell_df])
    return cell_df_new

def curate_fast_jumps(df,gr_max = 0.2, gr_min = -0.2):
     cell_list = df['cell_id_pos'].unique().tolist()
     cell_df_new = pd.DataFrame()
     for cell in cell_list:
         cell_df = df[df['cell_id_pos'] == cell]
         gr_smooth = cell_df['gr_smooth_norm']
         if np.any(gr_smooth>gr_max) or np.any(gr_smooth<gr_min):
             continue
         else:
             cell_df_new = pd.concat([cell_df_new,cell_df])
     return cell_df_new

def get_cell_coordinates(cell_stack_temp):   # given cell stack from supersegger, extract cell coordinates
    cell_coord = cell_stack_temp['coord'][0, 0][0, 0]
    rcm = cell_coord[8][0]
    A = cell_coord[0][0]
    r_off = cell_stack_temp['r_offset'][0, 0][0]
    angle = cell_coord[0][0][0]
    edge_flag = cell_stack_temp['edgeFlag'][0, 0][0, 0]
    BB = np.zeros(4, dtype=np.uint16)
    for k in range(4):
        BB[k] = cell_stack_temp['BB'][0][0][0][k]
    x1 = BB[0]-1
    y1 = BB[1]-1
    x2 = BB[0]+BB[2]
    y2 = BB[1]+BB[3]
    return A, rcm, r_off, angle, x1, x2, y1, y2, edge_flag

def get_average_gr(df, cell_type = 'type', px_area  = 0.004096, check_sis = False):   # function for estimating cell parameters, e.g., average growth rate, for each cell
    cell_df_new = pd.DataFrame()
    cells_list = list(df['cell_id_pos'].unique())
    for cell in cells_list:
        cell_df = df[df['cell_id_pos']==cell]
        G1 = False
        DNA_rep=False
        if 'G1_duration' in cell_df.columns:
            t_end_G1 = cell_df['t_end_G1'].iloc[0]
            G1_duration = cell_df['G1_duration'].iloc[0]
            G1_duration_norm = cell_df['G1_duration']/(cell_df['time_min'].iloc[-1]-cell_df['time_min'].iloc[0])
            A_G1 = cell_df['A_G1'].iloc[0]   
            av_gr_G1 = cell_df['av_gr_G1'].iloc[0]
            area_added_G1 = cell_df['area_added_G1'].iloc[0]
            G1 = True
            if 'DNA_rep_ini' in cell_df.columns:
                DNA_rep_ini = cell_df['DNA_rep_ini'].iloc[0]
                DNA_rep = True
        type_cell = cell_df[cell_type].iloc[0]
        exp = cell_df['exp'].iloc[0]
        pos = cell_df['xy'].iloc[0]
        tc = cell_df['time_min'].iloc[-1]-cell_df['time_min'].iloc[0]
        start_t = cell_df['time_min'].iloc[0]
        av_gr = np.log(cell_df['area'].iloc[-1]/cell_df['area'].iloc[0])/tc
        av_gr_2 = (cell_df['area'].iloc[-1] - cell_df['area'].iloc[0])/tc
        av_area = np.mean(cell_df['area'])
        Ab = cell_df['area'].iloc[0]
        Ad = cell_df['area'].iloc[-1]
        sister_ID=cell_df['sisterID'].iloc[0]
        area_added = cell_df['area'].iloc[-1] - cell_df['area'].iloc[0]
        cell_df = pd.DataFrame()
        cell_df['cell_id_pos'] = pd.Series(cell)
        cell_df['xy'] = pos
        cell_df['exp'] = exp
        cell_df[cell_type] = pd.Series(type_cell)
        cell_df['average_growth_rate'] = pd.Series(av_gr)
        cell_df['average_area'] = pd.Series(av_area)
        cell_df['tc'] = pd.Series(tc)
        cell_df['exp'] = exp
        cell_df['start_t'] = start_t
        cell_df['Ab'] = Ab
        cell_df['Ab_um'] = Ab*px_area
        cell_df['Ad'] = Ad
        cell_df['Ad_um'] = Ad*px_area
        cell_df['av_gr_2'] = av_gr_2
        cell_df['sisterID'] = sister_ID
        cell_df['area_added'] = area_added
        if G1:
            cell_df['t_end_G1'] = t_end_G1
            cell_df['G1_duration'] = G1_duration
            cell_df['G1_duration_norm'] = G1_duration_norm
            cell_df['A_G1'] = A_G1
            cell_df['av_gr_G1'] = av_gr_G1
            cell_df['area_added_G1'] = area_added_G1
            if DNA_rep:
                cell_df['DNA_rep_ini'] = DNA_rep_ini
        cell_df_new = pd.concat([cell_df_new,cell_df])
    if check_sis:
        cells_list = cell_df_new['cell_id_pos'].unique().tolist()
        cell_df_new_sis = pd.DataFrame()
        for cell in cells_list:
            cell_df = cell_df_new[cell_df_new['cell_id_pos']==cell]
            exp = cell_df['exp'].iloc[0]
            pos = cell_df['xy'].iloc[0]
            sisterID = cell_df['sisterID'].iloc[0]
            sister_name = format_number_as_string(sisterID) + '_' + pos + '_' +exp
            sister_cell_list = [s for s in cells_list if sister_name[1:] in s]
            if len(sister_cell_list)>0:
                sister_cell = sister_cell_list[0]
            else:
                continue
            sister_cell_df = cell_df_new[cell_df_new['cell_id_pos']==sister_cell]
            cell_df['tc_sis'] = sister_cell_df['tc'].iloc[0]
            cell_df['Ab_sis'] = sister_cell_df['Ab'].iloc[0]
            cell_df['average_growth_rate_sis'] = sister_cell_df['average_growth_rate'].iloc[0]
            cell_df['tc_sis_ratio'] = cell_df['tc'].iloc[0]/sister_cell_df['tc'].iloc[0]
            cell_df['Ab_sis_ratio'] =  cell_df['Ab'].iloc[0]/sister_cell_df['Ab'].iloc[0]
            cell_df['av_gr_sis_ratio'] =  cell_df['average_growth_rate'].iloc[0]/sister_cell_df['average_growth_rate'].iloc[0]
            cell_df['Ad_sis_ratio'] =  cell_df['Ad'].iloc[0]/sister_cell_df['Ad'].iloc[0]
            cell_df['Ad_sis'] = sister_cell_df['Ad'].iloc[0]
            cell_df['av_gr_2_sis'] = sister_cell_df['av_gr_2'].iloc[0]
            cell_df['div_ratio'] = cell_df['Ab']/(cell_df['Ab']+sister_cell_df['Ab'])
            if 'G1_duration' in cell_df.columns:
                cell_df['tc_end_G1_sis'] = sister_cell_df['t_end_G1'].iloc[0]
                cell_df['G1_duration_sis'] = sister_cell_df['G1_duration'].iloc[0]
                cell_df['G1_duration_norm_sis'] = sister_cell_df['G1_duration_norm'].iloc[0]
                cell_df['A_G1_sis'] = sister_cell_df['A_G1'].iloc[0]
                cell_df['av_gr_G1_sis'] = sister_cell_df['av_gr_G1'].iloc[0]
                cell_df['av_gr_2_sis'] = sister_cell_df['av_gr_2'].iloc[0]
            cell_df_new_sis = pd.concat([cell_df_new_sis,cell_df])
        return cell_df_new, cell_df_new_sis
    else: return cell_df_new


def bin_cell_trajectories(cell_features_df,bin_step=0.05, skip = 2,  cell_type = 'type', check_sis = False):     # function for binning cell trajectories
    cells_list = cell_features_df['cell_id_pos'].unique().tolist()
    bins = np.arange(start = -0.0001, stop = 1, step = bin_step)
    bins[-1]=1
    cell_df_new = pd.DataFrame()
    for cell in cells_list:
        cell_df = cell_features_df[cell_features_df['cell_id_pos'] == cell]
        t_len = len(cell_df['time_min'])-6-7
        time_norm_new = np.linspace(0,1,t_len)
        time_norm_new = np.concatenate((np.zeros(7)-1, time_norm_new, np.zeros(6)+100))
        bins_idx = np.digitize(time_norm_new,bins,right=True)+1
        bins_idx[:4] = 0
        bins_idx[-3:] = 23
        cell_df['bins_idx'] = bins_idx
        df_new = pd.DataFrame()
        df_new['area_smooth'] = cell_df.groupby('bins_idx')['area_smooth'].mean()   
        df_new['gr_smooth_norm'] =cell_df.groupby('bins_idx')['gr_smooth_norm'].mean()
        df_new['gr_smooth'] =cell_df.groupby('bins_idx')['gr_smooth'].mean()
        df_new['bins_idx'] = np.unique(bins_idx)
        df_new['bins_norm'] = np.unique(bins_idx)/max(bins_idx)
        df_new['cell_id_pos'] = cell_df['cell_id_pos'].iloc[0]
        df_new['sisterID'] = cell_df['sisterID'].iloc[0]
        df_new['xy'] = cell_df['xy'].iloc[0]
        df_new['exp'] = cell_df['exp'].iloc[0]
        df_new['start_t'] = cell_df['time_min'].iloc[0]
        tc = cell_df['time_min'].iloc[-1]-cell_df['time_min'].iloc[0]
        df_new['tc'] = tc
        if 'G1_duration' in cell_df.columns:
            df_new['G1_duration'] = cell_df['G1_duration'].iloc[0]  
            df_new['G1_duration_norm'] = cell_df['G1_duration'].iloc[0]/tc
        df_new[cell_type] = cell_df[cell_type].iloc[0]
        cell_df_new = pd.concat([cell_df_new,df_new.iloc[2:-2]])
    cells_list = cell_df_new['cell_id_pos'].unique().tolist()
    cell_df_new_2 = pd.DataFrame()
    if check_sis:
        for cell in cells_list:
            cell_df = cell_df_new[cell_df_new['cell_id_pos'] == cell]
            pos = cell_df['xy'].unique()[0]
            exp_name = cell_df['exp'].unique()[0]
            sisterID = cell_df['sisterID'].iloc[0]
            sister_name = format_number_as_string(sisterID) + '_' + pos + '_' +exp_name
            sister_cell_list = [s for s in cells_list if sister_name[1:] in s]
            if len(sister_cell_list)>0:
                sister_cell = sister_cell_list[0]
            else:
                continue
            sister_cell_df = cell_df_new[cell_df_new['cell_id_pos']==sister_cell]
            cell_df['ratio_gr_norm'] = cell_df['gr_smooth_norm']/sister_cell_df['gr_smooth_norm']
            cell_df['diff_gr_norm'] = cell_df['gr_smooth_norm']-sister_cell_df['gr_smooth_norm']
            cell_df['sister_gr_norm']  = sister_cell_df['gr_smooth_norm']
            cell_df_new_2 = pd.concat([cell_df_new_2,cell_df])
    return cell_df_new_2, cell_df_new

def import_cell_df(folder_path_list, exp_list,  fr_list, window_size, df_mark = 'cell_features_df', output_folder_df = 'output', time_lapse=True):
    cell_features_all_df = pd.DataFrame()
    cell_id_list = []
    for i in range(len(folder_path_list)):
        folder_path = folder_path_list[i]
        exp = exp_list[i]
        if type(fr_list) is list:
            fr = fr_list[i]
        else: fr = fr_list
        print('Import dataframe from experiment: ' + exp)
        save_dict_path = folder_path + '/output'
        save_dict_path_2 = folder_path + '/' +output_folder_df
        if not os.path.exists(save_dict_path):
            continue
        files_2 = os.listdir(save_dict_path_2)
        cell_df_list = [s for s in files_2 if df_mark in s and exp in s]
        for i in range(len(cell_df_list)):
            df_name = cell_df_list[i]
            idx_pos = df_name.find('xy')
            if can_be_number(df_name[idx_pos+2:idx_pos+4]):
                pos = df_name[idx_pos+2:idx_pos+4]
            else: pos = df_name[idx_pos+2:idx_pos+3]
            with open(save_dict_path_2+'/'+df_name, 'rb') as f:
                df = pickle.load(f)
            print('At position: ' + pos)
            if len(df)==0:
                continue
            if 'cell_id_pos' not in df.columns:
                df['exp'] = exp
                cell_id_pos = [s + '_' + pos + '_' + exp for s in df['cell_id'].tolist()]
                df.insert(loc=2, column='cell_id_pos', value=cell_id_pos)
                df['xy'] = pos
                if time_lapse:
                    df = df[df['motherID']>0]                                      # Select only cells with detected mother and daughter (full cell cycle)
                    df = df[df['daughters']>0]
            if len(df)==0:
                print('There are N=0 cells' + ' at pos ' + str(pos) + ' for experiment ' + exp)
            else:
                if time_lapse:
                    df, _ = get_normalized_cols(df, window_size)
                    df=time_align_healthy(df) 
                    df= add_features_growth(df, fr)
                      
                cell_list = df['cell_id_pos'].unique().tolist()
                print('There are N='+ str(len(cell_list)) + ' at pos ' + str(pos) + ' for experiment ' + exp)
                cell_features_all_df = pd.concat([cell_features_all_df, df])
                cell_id_list += cell_list
    return cell_features_all_df, cell_id_list


def local_bkg_sub_cp(img, mask, pos, exp, frame='',dilations=15, box_size = 128, sigma_=60, show=False):
    binary_mask = cp.asarray(mask>0)
    img_cp = cp.asarray(img)
    img_cell_free = img_cp.copy()
    mask_dil = ndi.binary_dilation(binary_mask,iterations=dilations,brute_force=True)
    img_cell_free[cp.where(mask_dil==1)]=0  ## set to 0 values of fluorescence image within dilated masks
    img_dim = img.shape
    n_row = int(img_dim[0]/box_size)   # round to smallest integer
    n_col = int(img_dim[1]/box_size)
    dim_new = n_row*box_size-sigma_, n_col*box_size-sigma_
    mean_back_img = cp.zeros(img_dim)
    median_back_img = cp.zeros(img_dim)
    for i in range(0,n_row):
        for j in range(0,n_col):
            box_temp = img_cell_free[i*box_size:(i+1)*box_size,j*box_size:(j+1)*box_size]
            mean_back_img[i*box_size:(i+1)*box_size,j*box_size:(j+1)*box_size]=cp.median(box_temp[box_temp>0])
            median_back_img[i*box_size:(i+1)*box_size,j*box_size:(j+1)*box_size]=cp.median(box_temp[box_temp>0])
    back_img = ndi.gaussian_filter(median_back_img, sigma = sigma_)
    back_sub = img_cp - back_img
    if show:
        fig=plt.figure(figsize=(12, 6))
        fig.suptitle('Background subtraction at frame '+frame+ ', at pos ' + pos +', ' + exp)
        
        plt.subplot(3, 2, 1)
        plt.imshow(img_cp[0:dim_new[0],0:dim_new[1]].get())
        plt.title('Original Image')
        
        plt.subplot(3, 2, 2)
        plt.imshow(mask_dil[0:dim_new[0],0:dim_new[1]].get())
        plt.title('Dilated masks')
        
        plt.subplot(3, 2, 3)
        plt.imshow(img_cell_free[0:dim_new[0],0:dim_new[1]].get())
        plt.title('Cell free image')
        
        plt.subplot(3, 2, 4)
        plt.imshow(median_back_img[0:dim_new[0],0:dim_new[1]].get())
        plt.title('Local median Bkg')
        
        plt.subplot(3, 2, 5)
        plt.imshow(back_img[0:dim_new[0],0:dim_new[1]].get())
        plt.title('Local smooth Bkg')
        
        plt.subplot(3, 2, 6)
        plt.imshow(back_sub[0:dim_new[0],0:dim_new[1]].get())
        plt.title('Bkg sub img')
        plt.tight_layout()
        plt.show()
    
    return back_sub[0:dim_new[0],0:dim_new[1]].get(), back_img[0:dim_new[0],0:dim_new[1]].get(), median_back_img[0:dim_new[0],0:dim_new[1]].get()