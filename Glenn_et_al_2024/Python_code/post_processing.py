# -*- coding: utf-8 -*-
"""
Code for quantifying cell parameters and plotting

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

df = cell_features_all_df.copy()
df = df[np.isin(df['type'],['Swarmer', 'Stalked'])]
df_binned, df_binned_sis = bin_cell_trajectories(df, bin_step=0.05, skip=2)
av_gr_df = get_average_gr(df)
av_st = av_gr_df[av_gr_df['type']=='Stalked']
av_sw = av_gr_df[av_gr_df['type']=='Swarmer']

''''
Example of plotting cell area, average and instantaneous growth rate (1/min)
'''
fig = plt.figure(figsize=(12, 10))
fig.suptitle(exp_name + ', N=' + str(len(df['cell_id_pos'].unique().tolist())))
ax1 = fig.add_subplot(2,2, 1)
x_var = 'time_min'
y_var = 'area_smooth'
sns.lineplot(df, x = x_var, y = y_var, estimator=None, lw=0.5, units="cell_id_pos", alpha=0.3,hue='type')
plt.xlabel('Time (min)', fontsize=14, fontweight = 'bold')
plt.ylabel('Area (px)', fontsize=14, fontweight = 'bold')

ax2 = fig.add_subplot(2,2, 2)
y_var = 'area_smooth'
x_var = 'bins_norm'
sns.lineplot(df_binned, x = x_var, y = y_var, hue = 'type')
plt.xlabel('Normalized cell cycle', fontsize=14, fontweight = 'bold')
plt.ylabel('Area (px)', fontsize=14, fontweight = 'bold')

ax3= fig.add_subplot(2,2,3)
plt.title('N_stalked='+str(len(av_st))+', N_swarmer='+str(len(av_sw)))
sns.histplot(av_gr_df, x = 'average_growth_rate', hue = 'type',edgecolor="white", stat='probability',kde=True, common_norm=False)
plt.xlabel('Average growth rate (1/min)', fontsize=14, fontweight = 'bold')
plt.ylabel('Probability', fontsize=14, fontweight = 'bold')
plt.xlim(0.002,0.010)

ax4= fig.add_subplot(2,2,4)
y_var = 'gr_smooth_norm'
x_var = 'bins_norm'
sns.lineplot(df_binned, x = x_var, y = y_var, hue = 'type')
plt.xlabel('Normalized cell cycle', fontsize=14, fontweight = 'bold')
plt.ylabel('Normalized growth Rate (1/min)', fontsize=14, fontweight = 'bold')

plt.tight_layout()
plt.show()

