# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 16:52:42 2025

@author: fragasso
"""
import pandas as pd
import numpy as np

'''This function is used to create new time axes based on evaluating peak values ('max' or 'min') of the markers (e.g. 'area' or 'SCF') over the current time_axis, using a defined evaluation time window
Can be used iteratively, e.g. first aligning by area over time_min, then align by SCF over t_align_area.
'''

def time_align(df, markers = ['area','SCF'], peak_values = ['max','max'], time_axis = 'time_min', window = [0,100]):
    def shift_time_axis(df_group, marker, window, time_axis, peak_value):
        df_group_w = df_group[(df_group[time_axis]>window[0]) & (df_group[time_axis]<window[1])]
        if peak_value == 'min':
            off_time = df_group.loc[df_group_w[marker].idxmin(), 'time_min']
        elif peak_value == 'max': 
            off_time = df_group.loc[df_group_w[marker].idxmax(), 'time_min']
        df_group['t_align_'+marker] = df_group['time_min'] - off_time
        return df_group
    for i in range(len(markers)):
        marker = markers[i]
        peak_value = peak_values[i]
        df = df.groupby('cell_id').apply(shift_time_axis, marker = marker, window = window, time_axis = time_axis, peak_value=peak_value).reset_index(drop=True)
    return df