# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 14:49:20 2025

@author: fragasso
"""



'''Function to calculate normalized growth rate 'gr_smooth_norm' used for curation steps. For this study, px_size = 0.065841  um/px.
It is calculated by smoothing the area signal using a 12-point moving average, and then calculating the first derivative of the smoothened signal.
'''
def add_features_growth(df, fr, px_area = px_size**2, window_size = 5):
    def dynamic_smooth(y, window_size=12):
        """
        Smooths the input array y using a dynamic moving average filter that adjusts
        the window size at the edges of the array.
        
        Parameters:
        - y: The input array to be smoothed.
        - window_size: The size of the moving window. Default is 3.
        
        Returns:
        - smoothed_y: The smoothed array.
        """
        n = len(y)
        smoothed_y = np.zeros(n)
        
        for i in range(n):
            # Adjust the window size for edges
            half_window = window_size // 2
            start_index = max(0, i - half_window)
            end_index = min(n, i + half_window + 1)
            
            # Compute the average for the current window
            window_elements = y[start_index:end_index]
            smoothed_y[i] = np.mean(window_elements)
        
        return smoothed_y
    df['area_smooth'] = df.groupby('cell_id')['area'].transform(lambda x: dynamic_smooth(x.values))
    cell_df_new = pd.DataFrame()
    for cell in list(df['cell_id'].unique()):
        cell_df = df[df['cell_id']==cell]
        cell_df['gr_smooth'] = np.insert(np.diff(cell_df['area_smooth']),0,0)/fr
        cell_df['gr'] = np.insert(np.diff(cell_df['area_um']),0,0)/fr
        cell_df['gr_smooth_norm'] = cell_df['gr_smooth']/cell_df['area_smooth']
        cell_df_new = pd.concat([cell_df_new,cell_df])
    df = cell_df_new.copy()
    return df

'''Curation functions'''

'''This functions eliminates trajectories shorter than min_fr, here set to 30 frames.'''
def curate_short_trajectories(df, min_fr=30):
    grouped = df.groupby('cell_id')
    filtered_df = grouped.filter(lambda x: (x['frame'].iloc[-1] - x['frame'].iloc[0]) > min_fr)
    return filtered_df


'''This functions is used to filter out part of the trajectory where growth rate is higher or lower than gr_max or gr_min. Values are adjusted based 
on specific dataset'''
def curate_fast_jumps(df, gr_max=0.2, gr_min=-0.2, window=[], x_var='time_min'):
    def filter_trajectory(group, gr_min, gr_max, x_var, window=[]):
        if window:
            window_mask = (group[x_var] >= window[0]) & (group[x_var] <= window[1])
            window_data = group[window_mask]
        else:
            window_data = group
        out_of_range = (window_data['gr_smooth_norm'] > gr_max) | (window_data['gr_smooth_norm'] < gr_min)
        if out_of_range.any():
            cutoff_index = out_of_range.idxmax()
            group = group.loc[:cutoff_index-1]
        return group
    grouped = df.groupby('cell_id')
    filtered_df = grouped.apply(filter_trajectory, gr_min=gr_min, gr_max=gr_max, x_var=x_var, window=window)
    filtered_df = filtered_df.drop(columns='cell_id', errors='ignore').reset_index(drop=False).drop(columns='level_1')
    return filtered_df


'''This functions is used to cells with too large area (>max_area) or too small (<min_area)'''
def curate_area(df, max_area = 2000, min_area = 300):
    grouped = df.groupby('cell_id')
    filtered_df = grouped.filter(lambda x: np.max(x['area'])<max_area and np.min(x['area'])>min_area)
    return filtered_df


'''This functions is used to curate cells with negative normalized growth rate prior to the frame of  AMP addition (fr_inj)'''
def curate_early_shrinkage(df):
    df_new = pd.DataFrame()
    for rep in df['rep'].unique().tolist():
        df_rep = df[df['rep']==rep]
        fr_inj = df_rep['fr_inj'].iloc[0]
        fr = (df_rep['time_min'].iloc[1]-df_rep['time_min'].iloc[0])/(df_rep['frame'].iloc[1]-df_rep['frame'].iloc[0])
        df_early = df_rep[(df_rep['time_min']/fr <= fr_inj - 6) & (df_rep['time_min']/fr > 7)]
        grouped = df_early.groupby('cell_id')
        bad_cells_df = grouped.filter(lambda x: np.min(x['gr_smooth_norm'])<0)
        df_temp = df_rep[~np.isin(df_rep['cell_id'],bad_cells_df['cell_id'])]
        df_new = pd.concat([df_new,df_temp])
    return df_new

'''This functions is used to curate cells with normalized growth rate <gr_min within the first window (=10) points of the cell trajectory'''
def curate_early_shrinkage_2(df, window=10, gr_min = 0.005):     ### eliminates cells not growing well for at least the first window time [min] of their cell cycle
    df_early = df[(df['time_min_aligned']<window) & (df['time_min_aligned']>0)]
    grouped_early = df_early.groupby('cell_id')
    filtered_early = grouped_early.filter(lambda x: np.min(x['gr_smooth_norm'])>gr_min)
    df = df[(df['cell_id'].isin(filtered_early['cell_id']))]
    return df

'''This functions is used to curate cells with normalized growth rate <gr_min or > gr_max within the last the window (=10) points of the cell trajectory'''
def curate_not_killed(df, window=10, gr_max=0.005, gr_min=-0.005):
    grouped = df.groupby('cell_id')
    filtered = grouped.filter(lambda x: (np.max(x['gr_smooth_norm'].iloc[-window:]) < gr_max) and  (np.min(x['gr_smooth_norm'].iloc[-window:]) > gr_min))
    return filtered