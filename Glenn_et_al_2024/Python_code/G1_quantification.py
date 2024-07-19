# -*- coding: utf-8 -*-
"""
Code for quantification of G1 in MipZ-eYFP expressing cells

@author: A.Fragasso
"""
import warnings
warnings.filterwarnings('ignore')
import skimage
from scipy import io
import os
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import warnings
warnings.filterwarnings('ignore')

runfile("path_to_analysis_functions_file", wdir="path_to_analysis_folder")

'''Import dataframes'''

runfile("path_to_analysis_functions_file", wdir="path_to_analysis_folder")

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

exp_name = 'experiment_name'
window_size = 5         # half of the sliding for calculating instantaneous GR
fr = 1.5       # time interval [min]
px_size = 0.064  # um/px  
cell_features_all_df, cell_id_list = import_cell_df(folder_path_list, exp_list, fr, window_size,  df_mark='cell_features_screened', output_folder_df = 'output_new')

'''Function to curate or re-assign MipZ particle 1 or particle 2 based on total intensities or areas'''
def curate_bad_particles(df, cell_to_remove = []):
    cells_list = df['cell_id_pos'].unique().tolist()
    df_new = pd.DataFrame()
    for cell in cells_list:
        if cell in cell_to_remove:
            continue
        cell_df = df[df['cell_id_pos'] == cell]
        cell_df_new = cell_df[np.isin(cell_df['num_p'],[2])].reset_index()
        idx_1 = cell_df_new[cell_df_new['area_particle_2']<15]['index']
        idx_2 = cell_df_new[cell_df_new['tot_fluor_particle_2']<600]['index']
        idx = pd.concat([idx_1,idx_2]).unique()
        cell_df['area_particle_2'].iloc[idx]=0
        cell_df['tot_fluor_particle_2'].iloc[idx]=0
        cell_df['mean_fluor_particle_2'].iloc[idx]=0
        cell_df['max_fluor_particle_2'].iloc[idx]=0
        cell_df['centroid_particle_2'].iloc[idx]=np.nan
        cell_df['dist_p'].iloc[idx]=0
        cell_df['num_p'].iloc[idx]=(cell_df['area_particle_1'].iloc[idx]>0).astype(int) + (cell_df['area_particle_2'].iloc[idx]>0).astype(int)
        idx_1 = cell_df_new[cell_df_new['area_particle_1']<10]['index']
        idx_2 = cell_df_new[cell_df_new['tot_fluor_particle_1']<500]['index']
        idx = pd.concat([idx_1,idx_2]).unique()
        cell_df['area_particle_1'].iloc[idx]=cell_df['area_particle_2'].iloc[idx]
        cell_df['area_particle_2'].iloc[idx]=0
        cell_df['tot_fluor_particle_1'].iloc[idx]=cell_df['tot_fluor_particle_2'].iloc[idx]
        cell_df['tot_fluor_particle_2'].iloc[idx]=0
        cell_df['mean_fluor_particle_1'].iloc[idx]=cell_df['mean_fluor_particle_2'].iloc[idx]
        cell_df['mean_fluor_particle_2'].iloc[idx]=0
        cell_df['max_fluor_particle_1'].iloc[idx]=cell_df['max_fluor_particle_2'].iloc[idx]
        cell_df['max_fluor_particle_2'].iloc[idx]=0
        cell_df['centroid_particle_1'].iloc[idx]=cell_df['centroid_particle_2'].iloc[idx]
        cell_df['centroid_particle_2'].iloc[idx] = np.nan
        cell_df['dist_p'].iloc[idx]=0
        cell_df['num_p'].iloc[idx]=(cell_df['area_particle_1'].iloc[idx]>0).astype(int) + (cell_df['area_particle_2'].iloc[idx]>0).astype(int) 
        df_new = pd.concat([df_new,cell_df])
    return df_new


'''Function to estimate G1 phase parameters: G1 duration, last frame at G1, area at end of G1, average growth rate during G1, area added during G1'''
def get_G1(df, px_area):
    cells_list = df['cell_id_pos'].unique().tolist()
    df_new = pd.DataFrame()
    def find_four_2_out_of_five(array, window=5):
        for i in range(len(array) - window):  
            if array[i] == 2:
                count_two = sum(array[i:i+5]==2)
                if count_two>window-2:
                    return i  
        return -1  
    for cell in cells_list:
        cell_df = df[df['cell_id_pos'] == cell]
        cell_df_new = cell_df[cell_df['num_p'].notna()].reset_index()
        num_p = np.array(cell_df_new['num_p']).T
        idx = find_four_2_out_of_five(num_p)
        if idx==-1:
            cell_df['t_end_G1'] = np.nan
            cell_df['G1_duration'] = np.nan
        else:
            t_end_G1 = cell_df_new['time_min'].iloc[np.max([idx-1,0])]
            t_G1 = t_end_G1 - cell_df['time_min'].iloc[0]
            cell_df['t_end_G1'] = t_end_G1
            cell_df['G1_duration'] = t_G1
            cell_df['A_G1'] = cell_df_new['area'].iloc[np.max([idx-1,0])]*px_area       
            cell_df['av_gr_G1'] = np.log(cell_df_new['area'].iloc[np.max([idx-1,0])]/cell_df_new['area'].iloc[0])/t_G1
            cell_df['area_added_G1'] = (cell_df_new['area'].iloc[np.max([idx-1,0])] - cell_df_new['area'].iloc[0])*px_area
            if cell_df['G1_duration'].iloc[0] == (cell_df['time_min'].iloc[-1]- cell_df['time_min'].iloc[0]):
                print(cell+' has G1 same as tc and will be excluded, for this idx='+str(idx))
        df_new = pd.concat([df_new,cell_df])
    return df_new


df = cell_features_all_df.copy()
df = df[np.isin(df['type'],['Swarmer', 'Stalked'])]
df = curate_bad_particles(df)
df = get_G1(df, px_area)
df=df[df['G1_duration'].notna()]
