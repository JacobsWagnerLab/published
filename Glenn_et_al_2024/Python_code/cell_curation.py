# -*- coding: utf-8 -*-
"""
Curation script. Plots cell trajectories and curate too short or too quikly varying trajectories 

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
cell_features_all_df, cell_id_list = import_cell_df(folder_path_list, exp_list, fr, window_size,  df_mark='cell_features', output_folder_df = 'output')

df = cell_features_all_df.copy()
df = curate_short_trajectories(df, area_inc = 150, time_int = 30)   # function to curate short trajectories
df = curate_fast_jumps(df,gr_max = 0.015, gr_min = -0.0025)


''''
Plot cell area, growth
'''
fig = plt.figure(figsize=(12, 10))
fig.suptitle(exp_name + ', N=' + str(len(df['cell_id_pos'].unique().tolist())))
ax1 = fig.add_subplot(2,2, 1)
x_var = 'time_min'
y_var = 'area_smooth'
df_frames = cell_start_end(df,'cell_id_pos',x_var)   # data frame with start and end frame for each cell
df_areas = cell_start_end(df,'cell_id_pos',y_var)
df_ends = pd.concat([df_frames,df_areas],axis=1)
sns.lineplot(df, x = x_var, y = y_var, estimator=None, lw=0.5, units="cell_id_pos", alpha=0.3,hue='type')
plt.xlabel('Time (min)', fontsize=14, fontweight = 'bold')
plt.ylabel('Area (px)', fontsize=14, fontweight = 'bold')

ax2 = fig.add_subplot(2,2, 2)
x_var = 'time_min_aligned'
y_var = 'gr_smooth_norm'
sns.lineplot(df, x = x_var, y = y_var, estimator=None, lw=0.2, units="cell_id_pos", color = 'black',alpha=0.2)
plt.xlabel('Time (min)', fontsize=14, fontweight = 'bold')
plt.ylabel('Normalized growth Rate (1/min)', fontsize=14, fontweight = 'bold')
plt.tight_layout()
plt.show()


