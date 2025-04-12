# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 16:09:22 2025

@author: fragasso
"""

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter, find_peaks


'''Combined derivatives (normalized growth rates) with weights for detecting growth rate stop'''
def savgol_normalized_derivative(time, area, window_length, polyorder=2):
    n = len(time)
    # Make sure window_length is valid (odd, >=3, <= n, etc.)
    # You can clamp if you want, e.g.:
    if window_length > n:
        window_length = n if n % 2 == 1 else n - 1
    if window_length < 3:
        window_length = 3
    if window_length % 2 == 0:
        window_length -= 1
    # Time step (assuming uniform sampling)
    dt = time[1] - time[0]
    # SavGol derivative of area
    dA_dt = savgol_filter(area, window_length=window_length, 
                          polyorder=polyorder, deriv=1, delta=dt)
    # SavGol smoothed area (0th derivative)
    A_smooth = savgol_filter(area, window_length=window_length, 
                             polyorder=polyorder, deriv=0)
    # Normalized derivative = (1/A) * (dA/dt)
    # Avoid division by zero if A_smooth has near-zero values
    norm_deriv = dA_dt / A_smooth
    return norm_deriv

def multi_scale_savgol_derivative(
    time, 
    area, 
    window_sizes,
    polyorder=2, 
    weight_mode='inverse',   # can be 'inverse', 'uniform', a custom list, etc.
    normalized=True
):
   
    window_sizes = sorted(window_sizes)

    dt = time[1] - time[0]  # assumes uniform spacing
    n = len(time)
    
    derivs_dict = {}
    for w in window_sizes:
        w_clamped = _clamp_sg_window(w, n)
        
        A_smooth = savgol_filter(area, w_clamped, polyorder, deriv=0)
        dA_dt = savgol_filter(area, w_clamped, polyorder, deriv=1, delta=dt)
        
        if normalized:
            deriv = dA_dt / A_smooth
        else:
            deriv = dA_dt
        
        derivs_dict[w_clamped] = deriv
    stacked = np.vstack([derivs_dict[w] for w in window_sizes])  
    weights = _build_weights(window_sizes, weight_mode)
    combined_deriv = np.average(stacked, axis=0, weights=weights)
    return combined_deriv, derivs_dict, weights

def _clamp_sg_window(w, n):
    """Ensure w is odd, >=3, and <= n."""
    if w > n:
        w = n
    if w < 3:
        w = 3
    if w % 2 == 0:
        w = max(3, w - 1)
    return w

def _build_weights(window_sizes, weight_mode):
    """Return a normalized array of weights for the given window_sizes."""
    n_scales = len(window_sizes)
    if isinstance(weight_mode, (list, np.ndarray)):
        w_arr = np.array(weight_mode, dtype=float)
        if len(w_arr) != n_scales:
            raise ValueError("Length of custom weights does not match number of window_sizes.")
    elif weight_mode == 'inverse':
        w_arr = 1.0 / np.array(window_sizes, dtype=float)
    elif weight_mode == 'uniform':
        w_arr = np.ones(n_scales, dtype=float)
    else:
        print(f"Unknown weight_mode='{weight_mode}', using uniform weights.")
        w_arr = np.ones(n_scales, dtype=float)
    w_sum = w_arr.sum()
    if w_sum == 0:
        w_arr = np.ones(n_scales, dtype=float) / n_scales
    else:
        w_arr /= w_sum
    return w_arr

def detect_shrink_stably(time, derivative, 
                         threshold_strict=0.0, 
                         threshold_stable=0.002, 
                         stable_window=30):
    n = len(time)
    if n < stable_window:
        return None  # not enough data to confirm stable crossing
    
    for i in range(1, n):
        # print(i)
        if derivative[i-1] >= threshold_strict and derivative[i] < threshold_strict:
            end_index = np.min([n, i + stable_window])
            if end_index > n:
                # not enough points left to confirm stability
                return None
            
            # If all derivative[i..(end_index-1)] < threshold_stable, we confirm
            if np.all(derivative[i:end_index] < threshold_stable):
                return time[i]
    
    return None


def process_cell_group(
    group, 
    window_sizes, 
    threshold_strict=0.0, 
    threshold_stable = 0.002,
    threshold_strict_2=0.01, 
    threshold_stable_2 = 0.012,
    stable_window = 30,
    polyorder=2,
    weight_mode='inverse', 
    normalized=True
):
   
    group = group.sort_values('time_min').copy()
    
    time_arr = group['time_min'].values
    area_arr = group['area'].values

    # Weighted multi-scale derivative
    combined_deriv, derivs_dict, weights = multi_scale_savgol_derivative(
        time_arr, 
        area_arr,
        window_sizes=window_sizes,
        polyorder=polyorder,
        weight_mode=weight_mode,
        normalized=normalized
    )
    
    event_time = detect_shrink_stably(time_arr, combined_deriv,
                                           threshold_strict=threshold_strict,
                                           threshold_stable=threshold_stable,
                                           stable_window=stable_window)
    
    event_time_2 = detect_shrink_stably(time_arr, combined_deriv,
                                           threshold_strict=threshold_strict_2,
                                           threshold_stable=threshold_stable_2,
                                           stable_window=stable_window)
    
    if event_time:
        offset_event_time = event_time - time_arr[0]
    else:
        offset_event_time = None
    if event_time_2:
        offset_event_time_2 = event_time_2 - time_arr[0]
    else:
        offset_event_time_2 = None
    if len(combined_deriv) == len(group):
        group['combined_deriv'] = combined_deriv
    else:
        group['combined_deriv'] = np.nan
    
    group['event_time'] = event_time
    group['offset_event_time'] = offset_event_time
    group['event_time_2'] = event_time_2
    group['offset_event_time_2'] = offset_event_time_2
    return group

def apply_weighted_derivative_detection(df, window_sizes, threshold_strict=0.0, threshold_stable = 0.002,threshold_strict_2=0.01, threshold_stable_2 = 0.012, stable_window = 30, polyorder=2, weight_mode='inverse', normalized=True):
    """
    Wrapper that applies process_cell_group via groupby('cell_id').
    """
    df_out = df.groupby('cell_id').apply(
        process_cell_group,
        window_sizes=window_sizes,
        threshold_strict=threshold_strict, 
        threshold_stable = threshold_stable,
        threshold_strict_2=threshold_strict_2, 
        threshold_stable_2 = threshold_stable_2,
        stable_window = stable_window,
        polyorder=polyorder,
        weight_mode=weight_mode,
        normalized=normalized
    ).reset_index(drop=True)
    return df_out


