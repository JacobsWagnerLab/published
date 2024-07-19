
"""
Code for importing data pre-processed by Supersegger-Omnipose and export:
- cell features dataframes
- dictionaries with lineage information 
- dictionaries with cropped images 

@author: A.Fragasso
"""
# import modules

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
from scipy.interpolate import splprep, splev, splrep, interp1d, UnivariateSpline, BSpline
from scipy.signal import lfiltic, lfilter, savgol_filter
from scipy.stats import spearmanr
from skimage import graph, morphology, measure
import imageio.v3 as iio
from scipy.ndimage import gaussian_filter, gaussian_filter1d, binary_fill_holes, binary_closing, binary_opening, center_of_mass
from scipy import ndimage
from PIL import Image as im
from skimage.morphology import skeletonize
from skimage.morphology import binary_closing
from skimage.filters import threshold_otsu, threshold_local
from skimage.segmentation import find_boundaries
from skimage import measure, morphology
from numpy.lib.stride_tricks import as_strided

runfile("path_to_analysis_functions_file", wdir="path_to_analysis_folder")

'''
Input global variables
'''
fr = 3     # [min] time between frames (write the specific one for each experiment)
px_size = 0.064 # px size
save = True
plot = True

'''Write paths to folders wherer the data are stored'''
folder_path_list = [
    'path_to_experiment_1',
    'path_to_experiment_2',
    'path_to_experiment_3'
    ]


'''Assign names to experiments'''
exp_list = [
    'name_experiment_1',
    'name_experiment_2',
    'name_experiment_3'
        ]

'''Create list of positions to scan through'''
pos_list = create_pos_list(1, 10, double_digit = True) + create_pos_list(1, 10, double_digit = False)

count=0     # counter for number of cells
for s in range(len(folder_path_list)):
    
    experiment_path = folder_path_list[s]
    exp_name = exp_list[s]
    print(exp_name)
    for pos in pos_list:
        
        '''
        Get paths to all folders
        '''
        folder_path = experiment_path + '/' +\
            pos  # path to position
        
        cells_path = folder_path+"/cell"
        if not os.path.exists(cells_path):
            continue
        elif len(os.listdir(cells_path))==0:
            continue
        
        print(pos)
        mask_path = folder_path + '/masks'
        phase_path = folder_path + '/phase'
        cells_list = os.listdir(cells_path)
        if os.path.exists(mask_path):
            masks_list = os.listdir(mask_path)
        phase_list = os.listdir(phase_path)
        '''
        Initialize dictionaries and lists
        '''
        cell_stack_dict = {}  # dictionary of cells parameters over time
        cell_ch_dict = {}  # dictionary of cells phase and mask over time
        cell_info_dict = {}   # dictionary of cell info: birth frame (starts from 1), death frame, divide or not, motherID, sisterID, daughterID
        good_cells_list = []
        cell_features_all_df = pd.DataFrame()
        for cell in cells_list:
            cell_id = cell[:11]
            good_cells_list.append(cell_id)
            C = io.loadmat(cells_path+'/'+cell)                         # dictionary with all meshes and fluor stats for one cell over time from Supersegger
            ### get cell info and store in cell_dict
            birth = C['birth'].flatten()[0]   
            death = C['death'].flatten()[0]
            divide = C['divide'].flatten()[0]
            motherID = C['motherID'].flatten()[0]
            sisterID = C['sisterID'].flatten()[0]
            daughterID = C['daughterID'].flatten()
            cell_stack = dict(enumerate(C['CellA'].flatten(), 1))
            cell_stack_dict[cell_id] = cell_stack
            if motherID == 0:
                gen = 0
                pole_1_label = -1
                pole_2_label = -1
                cell_type = 'first_gen_unknown'
                pole_inherited = np.nan                
            else:
                mother_cell = format_number_as_string(motherID)
                mother_name = [s[:11] for s in cells_list if mother_cell[1:] in s][0]
                birth_m, _, _, _, _, _, gen_m, pole_1_label_m, pole_2_label_m, _ , _, _, _ = cell_info_dict[mother_name]
                gen = gen_m + 1
               
                cell_stack_m = cell_stack_dict[mother_name]
                lm = list(cell_stack_m.keys())[-1]   # select last frame cell_stack from mother
                cell_stack_temp_m = cell_stack_m[lm]
                cell_coord_m = cell_stack_temp_m['coord'][0,0][0,0]
                xaxis_m = cell_coord_m[4]
                p1_m = [xaxis_m[0][0],xaxis_m[0][1]]
                p2_m = [xaxis_m[1][0],xaxis_m[1][1]]
            
            
                fd = list(cell_stack.keys())[0]   # select last frame cell_stack from mother
                cell_stack_temp_d = cell_stack[fd]
                cell_coord_d = cell_stack_temp_d['coord'][0,0][0,0]
                xaxis_d = cell_coord_d[4]
                p1_d = [xaxis_d[0][0],xaxis_d[0][1]]
                p2_d = [xaxis_d[1][0],xaxis_d[1][1]]
                
                d1 = np.sqrt((p1_d[0]-p1_m[0])**2 + (p1_d[1]-p1_m[1])**2) 
                d2 = np.sqrt((p2_d[0]-p2_m[0])**2 + (p2_d[1]-p2_m[1])**2) 
                
                d3 = np.sqrt((p1_d[0]-p2_m[0])**2 + (p1_d[1]-p2_m[1])**2) 
                d4 = np.sqrt((p2_d[0]-p1_m[0])**2 + (p2_d[1]-p1_m[1])**2) 
                
                if d1 + d2 <= d3 + d4:   ## if poles are passed on directly (no flipping of the major axis coordinates)
                    if d1 <= d2:        ## if p1_d is closer to p1_m than p2_d to p2_m, then consider p1_m as the pole inherited by p1_d (old pole), p2_m will be assigned to the other daughter, and p2_d will be the new pole 
                        pole_1_label = 0 
                        pole_2_label = 1
                        pole_inherited = pole_1_label_m
                        p_coord_inherited = p1_m
                    elif d2 < d1:       ## otherwise, consider p2_m as the inherited old pole and pass it on to p2_d, while p1_d will be hte new pole
                        pole_1_label = 1 
                        pole_2_label = 0
                        pole_inherited = pole_2_label_m
                        p_coord_inherited = p2_m
                else:  ## if poles coordinates are flipped
                    if d3 <= d4:  # p2_m will be passed on to p1_d as old pole
                        pole_1_label = 0
                        pole_2_label = 1
                        pole_inherited = pole_2_label_m
                        p_coord_inherited = p2_m
                    elif d4 < d3:
                        pole_1_label = 1 
                        pole_2_label = 0
                        pole_inherited = pole_1_label_m
                        p_coord_inherited = p1_m
                
                if pole_inherited == 0: 
                    cell_type = 'type_A'
                elif pole_inherited == 1:
                    cell_type = 'type_B'
                elif pole_inherited == -1:
                    cell_type = 'unknown'
     
            area_list=[]
            frame_list=[]
            edge_flag_list=[]
            p1_list = []
            p2_list = []
            p_inherited_list = []
            
            xaxis_off_list=[]
            daughterID_list = []
            for i in cell_stack:
                A, rcm, r_off, angle, x1, x2, y1, y2, edge_flag = get_cell_coordinates(cell_stack[i])
                cell_stack_temp = cell_stack[i]
                cropped_mask = cell_stack_temp['mask'][0][0]
                cropped_phase = cell_stack_temp['phase'][0][0]
                cell_coord = cell_stack_temp['coord'][0,0][0,0]
                xaxis = cell_coord[4]
                p1 = [xaxis[0][0],xaxis[0][1]]
                p2= [xaxis[1][0],xaxis[1][1]]
                daughterID_list.append(daughterID)
                p1_list.append(p1)
                p2_list.append(p2)
                if gen>0:   
                    p_inherited_list.append(p_coord_inherited)
                else : 
                    p_inherited_list.append(np.nan)
                cell_coord = cell_stack_temp['coord'][0,0][0,0]
                xaxis, yaxis = cell_coord[4], cell_coord[5]
                BB = np.zeros(4, dtype=np.uint16)
                for k in range(4):
                    BB[k] = cell_stack_temp['BB'][0][0][0][k]
                x1 = BB[0]
                y1 = BB[1]
                x2 = BB[0]+BB[2]
                y2 = BB[1]+BB[3]
                xaxis_off = np.array([[xaxis[0][0]-x1, xaxis[0][1]-y1], [xaxis[1][0]-x1, xaxis[1][1]-y1]])
                frame = i-2+birth
                cell_ch_dict[(cell_id, frame)] = cropped_mask, cropped_phase, cell_stack[i]
                area = np.size(np.where(cropped_mask>0)[0])
                edge_flag_list.append(edge_flag)
                frame_list.append(frame)
                area_list.append(area)
                xaxis_off_list.append(xaxis_off)
            
            cell_info_dict[cell_id] = birth, death, divide, motherID, sisterID, daughterID, gen, pole_1_label, pole_2_label, cell_type, p1_list, p2_list, xaxis_off_list
            cell_features_df = pd.DataFrame()
            cell_features_df['area'] = area_list
            cell_features_df['frame'] = frame_list
            cell_features_df['time_min'] = [fr*s for s in frame_list]
            cell_features_df['cell_id'] = cell_id
            cell_features_df['birth'] = birth
            cell_features_df['death'] = death
            cell_features_df['motherID'] = motherID
            cell_features_df['sisterID'] = sisterID
            cell_features_df['daughterID'] = daughterID_list
            cell_features_df['daughters'] = len(daughterID)
            cell_features_df['edge_flag'] = edge_flag_list
            cell_features_df['gen'] = gen
            cell_features_df['pole_1_label'] = pole_1_label
            cell_features_df['pole_2_label'] = pole_2_label
            cell_features_df['type_2'] = cell_type
            cell_features_df['p_inherited'] = pole_inherited
            cell_features_df['p1_coord'] = p1_list
            cell_features_df['p2_coord'] = p2_list
            cell_features_df['p_coord_inherited'] = p_inherited_list
            
            cell_features_all_df = pd.concat([cell_features_all_df,cell_features_df])
            
        if save and len(good_cells_list)>0:
            cell_ch_dict_name = exp_name + '_' + pos + '_cell_ch'
            cell_info_dict_name = exp_name + '_' + pos + '_cell_info'
            good_cells_name = exp_name + '_' + pos + '_good_cell_list'
            df_basename = exp_name+'_'+pos+'_cell_features_df'        
            save_dict_path = experiment_path + '/output'
            if not os.path.exists(save_dict_path):
                os.makedirs(save_dict_path)
        
            dict_cell_path = os.path.join(save_dict_path, cell_ch_dict_name+'.pkl')
            cell_info_path = os.path.join(save_dict_path, cell_info_dict_name+'.pkl')
            good_cells_path = os.path.join(save_dict_path, good_cells_name+'.json')
        
            with open(dict_cell_path, 'wb') as pickle_file:
                pickle.dump(cell_ch_dict, pickle_file)
            with open(cell_info_path, 'wb') as pickle_file:
                pickle.dump(cell_info_dict, pickle_file)
            with open(good_cells_path, 'w') as json_file:
                json.dump(good_cells_list, json_file)
            
            cell_features_all_df.to_pickle(save_dict_path+'/'+df_basename+'.pkl')
            
            count += len(good_cells_list)
            print('At position '+pos+' there are ' +
                  str(len(good_cells_list))+' cells')


print('In total there are '+
      str(count)+' cells')            
                