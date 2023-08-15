# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 17:14:28 2022

@author: Alexandros Papagiannakis, Christine Jacobs-Wagner lab, Sarafan ChEM-H, Stanford University, 2022

This script includes a collection of functions that were used for the analysis of the agarose pad data.
The functions are organized in the following fields:
    1. LOAD STATISTICS AND DATAFRAMES
    2. CREATE AND ANALYZE DEMOGRAPHS
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
from skimage import filters
import scipy
from scipy.stats import pearsonr
#from pims import ND2_Reader

#--------------  1. LOAD STATISTICS AND DATAFRAMES --------------#
def assemble_cell_stats(results_path):
    
    result_files = [s for s in os.listdir(results_path) if 'mean_df' in s]
    i = 0
    for s in result_files:
        if i == 0:
            mean_df = pd.read_pickle(results_path+'/'+s, compression='zip')
            i+=1
        elif i > 0:
            pre_df = pd.read_pickle(results_path+'/'+s, compression='zip')
            mean_df = pd.concat([mean_df, pre_df])
    return mean_df

def assemble_cell_lengths(results_path):
    
    result_files = [s for s in os.listdir(results_path) if 'cell_length' in s]
    i = 0
    for s in result_files:
        if i == 0:
            length_dict = pickle.load(open(results_path+'/'+s, 'rb'))
            i+=1
        elif i > 0:
            pre_dict = pickle.load(open(results_path+'/'+s, 'rb'))
            length_dict = {**length_dict, **pre_dict}
    return length_dict

def assemble_oned_coords(results_path):
    
    result_files = [s for s in os.listdir(results_path) if 'oned_coords_dict' in s]
    i = 0
    for s in result_files:
        if i == 0:
            oned_dict = pickle.load(open(results_path+'/'+s, 'rb'))
            i+=1
        elif i > 0:
            pre_dict = pickle.load(open(results_path+'/'+s, 'rb'))
            oned_dict = {**oned_dict, **pre_dict}
    return oned_dict

def assemble_bad_cells(results_path):
    
    result_files = [s for s in os.listdir(results_path) if 'bad_cells' in s]
    i = 0
    for s in result_files:
        if i == 0:
            bad_cells_list = pickle.load(open(results_path+'/'+s, 'rb'))
            i+=1
        elif i > 0:
            pre_list = pickle.load(open(results_path+'/'+s, 'rb'))
            bad_cells_list = bad_cells_list+pre_list
    return bad_cells_list

def assemble_particle_stats(results_path):
    
    result_files = [s for s in os.listdir(results_path) if 'particles_df' in s]
    i = 0
    for s in result_files:
        if i == 0:
            par_df = pd.read_pickle(results_path+'/'+s, compression='zip')
            i+=1
        elif i > 0:
            pre_df = pd.read_pickle(results_path+'/'+s, compression='zip')
            par_df = pd.concat([par_df, pre_df])
    return par_df

def get_twod_dataframe(results_path):
    
    oned_coords_dict = assemble_oned_coords(results_path)
    
    i=0
   
    for cl in oned_coords_dict:
        if i==0:
            oned_coords_df = oned_coords_dict[cl]
            oned_coords_df['cell_id'] = cl
            i+=1
        elif i>0:
            temp_df = oned_coords_dict[cl]
            temp_df['cell_id'] = cl
            oned_coords_df = pd.concat([oned_coords_df, temp_df])
#            print(row.cell_id)

    oned_coords_df.to_pickle(results_path+'/oned_coords_df', compression='zip')
    
    return oned_coords_df



#-----------------   2. CREATE AND ANALYZE DEMOGRAPHS   ----------------#
def get_demograph_arrays(results_path, bin_n, rol_win, channel, bad_cells):
    
    cell_length_dict = assemble_cell_lengths(results_path)
    oned_coords_dict = assemble_oned_coords(results_path)
    if bad_cells == True:
        bad_cells_list = assemble_bad_cells(results_path)
    else:
        bad_cells_list = []
    # sort keys based on the values
    sorted_keys = sorted(cell_length_dict, key=cell_length_dict.get)
    print(len(list(cell_length_dict)), 'cells segmented')
    print(len(bad_cells_list), 'bad segmentations')
    mean_arrays = []
    z_arrays = []
    mean_cor_arrays = []
    length_arrays = []
    for cl in sorted_keys:
        if cl in oned_coords_dict and cl not in bad_cells_list:
#            print(cl, 'getting 1D fluorescence arrays')
            oned_coords = oned_coords_dict[cl]
            oned_coords = oned_coords[oned_coords.width.between(-3,3)]

            oned_coords['scaled_length_bin'] = pd.cut(oned_coords.scaled_length, bins=np.arange(-1,1+2/bin_n,2/bin_n), right=True, include_lowest=True, labels=False)
            mean_df = oned_coords.groupby('scaled_length_bin').mean().reindex(list(range(0,bin_n)))
            mean_df = mean_df.interpolate()
            smooth_df = mean_df.rolling(rol_win, min_periods=1, center=True).mean()
            mean_arrays.append(smooth_df[channel+'_fluor'].tolist())
            mean_cor_arrays.append((smooth_df[channel+'_fluor']/smooth_df[channel+'_fluor'].mean()).tolist())
            z_score = (smooth_df[channel+'_fluor']-smooth_df[channel+'_fluor'].mean())/smooth_df[channel+'_fluor'].std()
            z_arrays.append(z_score.tolist())
            length_arrays.append(cell_length_dict[cl])
    
    return mean_arrays, z_arrays, mean_cor_arrays, length_arrays


def plot_demograph(z_array, bin_n_2, smoothing_factor, cmap, v_range, channel, save_path):
    img = np.array(z_array)
    img_df = pd.DataFrame(img)
    img_df['fluor_bin'] = pd.cut(img_df.index, bins=bin_n_2[1], include_lowest=True, right=True)
    img = img_df.groupby('fluor_bin').mean().values
    img_df = pd.DataFrame(np.transpose(img))
    img_df['fluor_bin'] = pd.cut(img_df.index, bins=bin_n_2[0], include_lowest=True, right=True)
    img = img_df.groupby('fluor_bin').mean().values
    
    plt.figure(figsize=(5*int(bin_n_2[1]/bin_n_2[0]),5))
    plt.imshow(filters.gaussian(img,smoothing_factor), cmap=cmap, vmin=v_range[0], vmax=v_range[1])
    cbar = plt.colorbar()
    cbar.set_label(str(channel)+'statistic', rotation=90, size=14)

    # plt.axhline([bin_n_2[0]/2], color='black', linewidth=2) 
    plt.xlabel('Binned cells sorted by their length')
    plt.ylabel('Cell length bin')
    if os.path.isdir(save_path):
        plt.savefig(save_path+'/'+channel+'snapshots_demograph.eps')
    plt.show()        
           

def get_fluorescence_at_midcell(results_path):
    
    cell_length_dict = assemble_cell_lengths(results_path)
    oned_coords_dict = assemble_oned_coords(results_path)
    # sort keys based on the values
    sorted_keys = sorted(cell_length_dict, key=cell_length_dict.get)
    rpla_mid = []
    hupa_mid = []
    cell_length = []
    cell_ids = []
    rpla_cell = []
    hupa_cell = []
    rpla_std = []
    hupa_std = []
    for cl in sorted_keys:
        if cl in oned_coords_dict:
            oned_df = oned_coords_dict[cl]
            center_df = oned_df[oned_df.scaled_length.abs().between(0,0.1)]
            center_df = center_df[center_df.width.between(-3,3)]
            
            rpla_mid.append(center_df['FITC_fluor'].mean())
            hupa_mid.append(center_df['mCherry_fluor'].mean())
            
            rpla_cell.append(oned_df['FITC_fluor'].mean())
            hupa_cell.append(oned_df['mCherry_fluor'].mean())
            
            rpla_std.append(oned_df['FITC_fluor'].std())
            hupa_std.append(oned_df['mCherry_fluor'].std())
            
            cell_length.append(cell_length_dict[cl])
            cell_ids.append(cl)
            
    mid_df = pd.DataFrame()
    mid_df['cell_id'] = cell_ids
    mid_df['mid_rpla'] = rpla_mid
    mid_df['mid_hupa'] = hupa_mid
    mid_df['length_px'] = cell_length
    mid_df['cell_rpla'] = rpla_cell
    mid_df['cell_hupa'] = hupa_cell
    mid_df['std_rpla'] = rpla_std
    mid_df['std_hupa'] = hupa_std
    
    return mid_df


def get_demograph_statistics(arrays_dict, mean_df, ribo_thres, nuc_thres, width_range=[45,56],  
                             channels=['RplA-mEOS2 z-score', 'DAPI z-score']):
    """
    This function returns the relative timing of polysome accumulation, nucleoids segregation,
    as well as the maximum level of polysome accumulation at mid-cell, from demographs.
    """
    max_ribo_dict = {}
    nucleoid_seg_timing = {}
    ribo_acum_timing = {}
    max_ribo_der = {}
    min_nuc_der = {}
    
    cell_area_dict = mean_df.groupby('condition').mean().cell_area_px.to_dict()
    if mean_df['strain'].values[0]==7323:
        exception_cond=mean_df.cell_area_px.mean()
    else:
        exception_cond=0
    
#     col_n = 6
#     if len(list(arrays_dict[channels[0]].keys()))%col_n==0:
#         row_n = int(len(list(arrays_dict[channels[0]].keys()))/col_n)
#     else:
#         row_n = int(len(list(arrays_dict[channels[0]].keys()))/col_n)+1
    
#     fig = plt.figure(figsize=(col_n*5,row_n*3))
    
#     gs = matplotlib.gridspec.GridSpec(row_n, col_n, width_ratios = col_n*[1], height_ratios = row_n*[1], hspace=0.4, wspace=0.6)
#     gss = []
#     for r in range(row_n):
#         for c in range(col_n):
#             gss.append(gs[r,c])
    
#     condition_axes = {}
    i = 0
    for cond in arrays_dict[channels[0]]:
        
#         condition_axes[cond+'_l']=plt.subplot(gss[i])
#         condition_axes[cond+'_r'] = condition_axes[cond+'_l'].twinx()
#         condition_axes[cond+'_r'].tick_params(axis='y', direction='in', size=10, labelsize=12, labelcolor='tomato', color='tomato')
#         condition_axes[cond+'_l'].tick_params(axis='y', direction='in', size=10, labelsize=12, labelcolor='royalblue', color='royalblue')
#         condition_axes[cond+'_l'].tick_params(axis='x', direction='in', size=10, labelsize=12, labelcolor='black', color='black')
# #         condition_axes[cond+'_r'].set_ylim(-1,1)
# #         condition_axes[cond+'_l'].set_ylim(0,1.5)
        
#         condition_axes[cond+'_l'].set_title(cond, fontname='Arial', fontweight='bold', fontsize=20)
        
#         print(cond)
        
        cell_n = np.array(arrays_dict[channels[0]][cond]).shape[0]
#         print(cell_n,'cells')
        five_per = int(0.05*cell_n)
#         five_per = 0
        nint_per = int(0.95*cell_n)
#         nint_per=cell_n
        
        
        data_ribo_mid = np.median(np.array(arrays_dict[channels[0]][cond])[five_per:nint_per,width_range[0]:width_range[1]], axis=1)
        data_nuc_mid = np.median(np.array(arrays_dict[channels[1]][cond])[five_per:nint_per,width_range[0]:width_range[1]], axis=1)
#         data_ribo_mid = np.percentile(np.array(arrays_dict[channels[0]][cond])[:,width_range[0]:width_range[1]], 75, axis=1)
#         data_nuc_mid = np.percentile(np.array(arrays_dict[channels[1]][cond])[:,width_range[0]:width_range[1]], 25, axis=1)

#        data_length = len(list(data_ribo_mid))
#        print(data_length)
#        data_ribo_mid = np.percentile(np.array(arrays_dict[channels[0]][cond])[:,width_range[0]:width_range[1]], 75, axis=1)
#        data_nuc_mid = np.percentile(np.array(arrays_dict[channels[1]][cond])[:,width_range[0]:width_range[1]], 25, axis=1)

        data_df = pd.DataFrame()
        data_df['data_nuc'] = data_nuc_mid#-((data_nuc_left+data_nuc_right)/2)
        data_df['data_ribo'] = data_ribo_mid#-((data_ribo_left+data_ribo_right)/2)
        data_df['cell_index'] = np.arange(0,data_df.shape[0])/data_df.shape[0]
        data_df['cell_bin'] = pd.cut(data_df.cell_index,100)
        mean_df = data_df.groupby('cell_bin').mean()

        fit_ribo = scipy.interpolate.UnivariateSpline(mean_df.cell_index, mean_df.data_ribo, k=5, s=6) # 25
        fit_nuc = scipy.interpolate.UnivariateSpline(mean_df.cell_index, mean_df.data_nuc, k=4, s=2) # 50

        fit_ribo_der = fit_ribo.derivative()
        fit_x = np.arange(0,1.0001,0.0001)
        ribo_der = fit_ribo_der(fit_x)
#         ribo_x =  solve_for_y(fit_x, ribo_der, np.max(ribo_der))[0]
        ribo_sol = fit_ribo(fit_x)
        nuc_sol = fit_nuc(fit_x)
        
        fit_df = pd.DataFrame()
        fit_df['x'] = fit_x
        fit_df['r'] = ribo_sol
        fit_df['n'] = nuc_sol
        
        early_fit = fit_df[fit_df.x.between(0,0.8)]
        min_ribo = early_fit.r.min()
        early_fit['ribo_dif'] = (early_fit.r-min_ribo).abs()
        min_index = early_fit[early_fit.ribo_dif==early_fit.ribo_dif.min()].x.values[0]
        early_fit = fit_df[fit_df.x>=min_index]
        max_ribo = early_fit.r.max()
        early_fit['ribo_dif'] = (early_fit.r-max_ribo).abs()
        max_index = early_fit[early_fit.ribo_dif==early_fit.ribo_dif.min()].x.values[0]
        early_fit = fit_df[fit_df.x.between(min_index, max_index)]
        
#         ribo_half = (early_fit.r.min()+early_fit.r.max())/2
        ribo_half = ribo_thres
#         print('ribo half z:', ribo_half)
        early_fit['ribo_dif'] = (early_fit.r-ribo_half).abs()
        ribo_x = early_fit[early_fit.ribo_dif==early_fit.ribo_dif.min()].x.values[0]
#         print('cell cycle at half ribo:', ribo_x)
#         nuc_half = (fit_df.n.min()+fit_df.n.max())/2
        nuc_half=nuc_thres
#         print('nuc half z:', nuc_half)
        fit_df['nuc_dif'] = (fit_df.n-nuc_half).abs()
        nuc_x = fit_df[fit_df.nuc_dif==fit_df.nuc_dif.min()].x.values[0]
#         print('cell cycle at half ribo:', nuc_x)
        
        fit_nuc_der = fit_nuc.derivative()
        nuc_der = fit_nuc_der(fit_x)

#         condition_axes[cond+'_l'].plot(mean_df.cell_index, mean_df.data_ribo, color='powderblue')
#         condition_axes[cond+'_r'].plot(mean_df.cell_index, mean_df.data_nuc, color='lightsalmon')
#         condition_axes[cond+'_l'].plot(fit_x, ribo_sol, linewidth=1, color='royalblue')
#         condition_axes[cond+'_r'].plot(fit_x, nuc_sol, linewidth=1, color='tomato')
#         condition_axes[cond+'_r'].plot(nuc_x, nuc_half, 'o', markersize=10, color='tomato')
#         condition_axes[cond+'_l'].plot(ribo_x, ribo_half, 'o', markersize=10, color='royalblue')
#         condition_axes[cond+'_l'].axhline(np.max(ribo_sol), linewidth=0.5, color='royalblue')
        

#         max_ribo_dict[cond] = np.max(ribo_sol)
        max_ribo_dict[cond] = max_ribo
        nucleoid_seg_timing[cond] = nuc_x
        ribo_acum_timing[cond] = ribo_x
        max_ribo_der[cond] = np.max(ribo_der)
        min_nuc_der[cond] = np.min(nuc_der)
        i+=1
#     plt.show()
    
    results_df = pd.DataFrame()
    results_df['condition'] = list(max_ribo_dict.keys())
    results_df['max_ribo'] = results_df.condition.map(max_ribo_dict)
    results_df['nuc_time'] = results_df.condition.map(nucleoid_seg_timing)
    results_df['ribo_time'] = results_df.condition.map(ribo_acum_timing)
    results_df['ribo_der'] = results_df.condition.map(max_ribo_der)
    results_df['nuc_der'] = results_df.condition.map(min_nuc_der)
    results_df['cell_area_px'] = results_df.condition.map(cell_area_dict)
    if exception_cond!=0:
        results_df['cell_area_px']=exception_cond
    
    return results_df


def get_depletion(arrays_dict, mean_df):
    """
    This function returns the relative polysome accumulation at mid-cell and 
    the relative nucleoid depletion from the same region from the average demograph projection.
    """
    def get_linescan_difference(linescan, species):
        middle = list(linescan[48:53])
#         edges = list(linescan[25:31])+list(linescan[70:76])
        edgel = list(linescan[20:51])
        edger = list(linescan[50:81])
        
        if species == 'ribo':
            return np.max(middle) - np.mean([np.min(edgel),np.min(edger)])
        elif species == 'nuc':
            return np.min(middle) - np.mean([np.max(edgel),np.max(edger)])
    
    def get_linescan_middle(linescan):
        middle = list(linescan[48:53])
        return np.mean(middle)
    

    ribo_list = []
    hu_list = []
    pearson_list = []
    significance_list = []
    area_list = []
    
    ribo_key = list(arrays_dict.keys())
    if 'DAPI z-score' in ribo_key:
        ribo_key.remove('DAPI z-score')
        nuc_key = 'DAPI z-score'
    elif 'HupA-mCherry z-score' in ribo_key:
        ribo_key.remove('HupA-mCherry z-score')
        nuc_key = 'HupA-mCherry z-score'
    ribo_key = ribo_key[0]
    
    
    for cond in arrays_dict[nuc_key]:
        print(mean_df.strain.values[0], cond)
        area = mean_df[mean_df.condition==cond].cell_area_px.mean()*0.066**2
        ribo_linescan = np.mean(arrays_dict[ribo_key][cond], axis=0)
        hu_linescan = np.mean(arrays_dict[nuc_key][cond], axis=0)
#         plt.plot(ribo_linescan)
#         plt.plot(hu_linescan)
#         plt.show()
        pears = pearsonr(ribo_linescan[25:76], hu_linescan[25:76])
#         if pears[1]<0.05:
        pearson_list.append(pears[0])
        significance_list.append(pears[1])
#         else:
#             pearson_list.append(np.nan)
        ribo_list.append(get_linescan_difference(ribo_linescan, 'ribo'))
        hu_list.append(get_linescan_difference(hu_linescan, 'nuc'))
        area_list.append(area)
    return ribo_list, hu_list, pearson_list, significance_list, area_list
