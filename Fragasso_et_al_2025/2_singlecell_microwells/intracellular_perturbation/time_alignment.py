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

def time_align(df, markers=['area', 'SCF'], peak_values=['max', 'max'], time_axis='time_min_aligned', window=[0, 100]):
    df = df.copy()
    for marker, peak_value in zip(markers, peak_values):
        def off_time_for_group(g):
            gw = g[(g[time_axis] > window[0]) & (g[time_axis] < window[1])]
            if len(gw) == 0:
                return np.nan
            idx = gw[marker].idxmin() if peak_value == 'min' else gw[marker].idxmax()
            return g.loc[idx, 'time_min']
        off_times = df.groupby('cell_id').apply(off_time_for_group, include_groups=False)
        df['t_align_' + marker] = df['time_min'] - df['cell_id'].map(off_times)
    return df