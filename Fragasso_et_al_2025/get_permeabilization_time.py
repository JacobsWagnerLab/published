# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 15:55:57 2025

@author: fragasso

Function to get the timings of IM and OM permeabilization based on flurescent protein (mScarlet-I, mTagBFP2) signals
"""
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

 """
 Smooths the input array y using a dynamic moving average filter that adjusts
 the window size at the edges of the array.
 """

def dynamic_smooth(y, window_size=12):
   
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


'''This function was used to compute the timings of IM and OM permeabilization, given dataframe containing mean fluorescences of the three channels 
(mScarlet-I, mTagBFP2, and SYTOX Green), their derivatives, as well as peri_cyto_diff(mScarlet-I), the difference between mean periplasmic and cytoplasmic
mScarlet-I signal.

Returns the tp_df dataframe that contains all cell informations from df, plus timings of IM, OM permeabilization, growth arrest and start of inhibition.
'''
def get_permeabilization_time(df, ch_info_df, ch_info_list = ['mScarlet','BFP','Sytox'],   w_Sy = 5, fr = 1, peri_fluor = 'fluor1', thr_IM = -25, t_max = 241):
    tp_df = pd.DataFrame()
    for cell in df['cell_id'].unique().tolist():
        cell_df = df[df['cell_id']==cell].reset_index(drop=True)
        cell_df = cell_df[cell_df['time_min']<=t_max]
        
        tp_temp = pd.DataFrame()
        tp_temp['cell_id'] = pd.Series(cell)
        tp_temp['AMP'] = cell_df['AMP'].iloc[0]
        tp_temp['rep'] = cell_df['rep'].iloc[0]
        tp_temp['fr_inj'] = cell_df['fr_inj'].iloc[0]
        ch_list = ch_info_df[(ch_info_df['AMP']==tp_temp['AMP'].iloc[0]) & (ch_info_df['rep']==tp_temp['rep'].iloc[0])]['mean_fluors'].iloc[0]      
        peri_fluor = ch_info_df[(ch_info_df['AMP']==tp_temp['AMP'].iloc[0]) & (ch_info_df['rep']==tp_temp['rep'].iloc[0])]['peri_fluor'].iloc[0]
        cell_df['mean_'+ch_info_list[0]] = cell_df[ch_list[0]]
        cell_df['mean_'+ch_info_list[1]] = cell_df[ch_list[1]]
        cell_df['mean_'+ch_info_list[2]] = cell_df[ch_list[2]]
        cell_df['mean_'+ch_info_list[2]+'_smooth'] = dynamic_smooth(cell_df['mean_'+ch_info_list[2]], window_size=3)
        cell_df['gr_'+ch_info_list[0]] = cell_df['gr_'+ch_list[0]+'_norm']
        cell_df['gr_'+ch_info_list[1]] = cell_df['gr_'+ch_list[1]+'_norm']
        cell_df['gr_'+ch_info_list[2]] = cell_df['gr_'+ch_list[2]+'_norm']
        
    
        tp_temp['t_area_max'] = pd.Series(cell_df['time_min'].iloc[np.argmax(cell_df['area'])]) 
        tp_temp['fr_area_max'] = pd.Series(cell_df['frame'].iloc[np.argmax(cell_df['area'])]) 
        tp_temp['area_max'] = pd.Series(np.max(cell_df['area']))
        tp_temp['t_gr_stop'] = pd.Series(cell_df['event_time'].iloc[0])
        tp_temp['t_gr_stop_offset'] = cell_df['event_time'].iloc[0]- tp_temp['fr_inj'].iloc[0]*fr
        tp_temp['t_gr_stop_early'] = pd.Series(cell_df['event_time_2'].iloc[0])
        if cell_df['time_min'].iloc[-1] - tp_temp['t_gr_stop'].iloc[0] <40:    ## there needs to be at least 40 minutes after the growth has stopped
            continue
        tp_temp['t_gr_stop_early_offset'] = tp_temp['t_gr_stop_early'].iloc[0]- tp_temp['fr_inj'].iloc[0]*fr
        tp_temp['t_aligned_gr_stop'] = pd.Series(cell_df['offset_event_time'].iloc[0])
        tp_temp['t_aligned_gr_stop_early'] = pd.Series(cell_df['offset_event_time_2'].iloc[0])
        
        
        # Search for early IM permeabilization using peri_cyto_diff function
        peri_cyto_diff = cell_df[peri_fluor+'_peri_cyto_diff']
        gr_peri_cyto_diff = cell_df['gr_peri_cyto_diff_'+peri_fluor+'_norm']
        
        inverted_signal = -gr_peri_cyto_diff
        peaks, properties = find_peaks(inverted_signal, height=0.15)
        negative_peak_indices = peaks
        negative_peak_amplitudes = gr_peri_cyto_diff.iloc[negative_peak_indices]
        valid_peaks = []
        for peak_index in peaks:
            search_region = peri_cyto_diff.iloc[peak_index -1  : np.min([peak_index + 20,len(peri_cyto_diff)])]
            search_stable_region = peri_cyto_diff.iloc[peak_index+5: np.min([peak_index + 30,len(peri_cyto_diff)])]
            if (search_region < thr_IM).any() and not (search_stable_region>thr_IM+10).any():
                valid_peaks = [peak_index]
                if len(valid_peaks)>0:
                    break
        if len(valid_peaks)>0:
            IM_idx = valid_peaks[0] -cell_df.index.tolist()[0]
            tp_temp['tp_IM'] = pd.Series(cell_df['time_min'].iloc[IM_idx])
            tp_temp['tp_aligned_IM'] = pd.Series(cell_df['time_min_aligned'].iloc[IM_idx])
            tp_temp['fp_IM'] = pd.Series(cell_df['frame'].iloc[IM_idx])
            tp_temp['pc_diff_IM'] = pd.Series(cell_df[peri_fluor+'_peri_cyto_diff'].iloc[IM_idx])
            
        # Go through all the three channels to find timings of permeabilization. For SYTOX Green (two thresholds were used for mild or fast uptake)
        # though not used in the final analysis
        for i in range(len(ch_info_list)):
            mean_fluor = cell_df['mean_'+ch_info_list[i]]
            if ch_info_list[i] == 'Sytox':                  ## probe permeabilization to Sytox Green
                fluor_thr1 =50
                fluor_thr2 = 1000
    
                # Detect the first event where the signal exceeds fluor_thr1
                if np.any(mean_fluor >= fluor_thr1):
                    idx1 = np.where(mean_fluor >= fluor_thr1)[0][0]
                    sytox1 = mean_fluor.iloc[idx1]
                    tp_temp['tp_Sytox1'] = pd.Series(cell_df['time_min'].iloc[idx1])
                    tp_temp['tp_aligned_Sytox1'] = pd.Series(cell_df['time_min_aligned'].iloc[idx1])
                    tp_temp['Sytox1'] = pd.Series(cell_df['mean_'+ch_info_list[i]].iloc[idx1])
                else:
                    idx1 = -1
                    sytox1 = 0
            
                # Detect the second event where the signal exceeds fluor_thr2
                gr_fluor_norm = cell_df['gr_'+ch_info_list[i]]
                peaks, properties = find_peaks(gr_fluor_norm, height=0.04)
                peak_indices = peaks
                peak_amplitudes = gr_fluor_norm.iloc[peak_indices]
                valid_peaks = []
                for peak_index in peaks:
                    search_region = mean_fluor.iloc[peak_index -1  : np.min([peak_index + 10,len(mean_fluor)])]
                    if (search_region > fluor_thr2).any():
                        valid_peaks = search_region[search_region > fluor_thr2].index.tolist()
                        if len(valid_peaks)>0:
                            break
                if len(valid_peaks)>0:
                    idx2 = valid_peaks[0] -cell_df.index.tolist()[0]
                    sytox2 = mean_fluor.iloc[idx2]
                    tp_temp['tp_Sytox2'] = pd.Series(cell_df['time_min'].iloc[idx2])
                    tp_temp['tp_aligned_Sytox2'] = pd.Series(cell_df['time_min_aligned'].iloc[idx2])
                    tp_temp['Sytox2'] = pd.Series(cell_df['mean_'+ch_info_list[i]].iloc[idx2])
                else:
                    idx2 = -1
                    sytox2 = 0

            else:
                mean_fluor_norm = cell_df['mean_'+ch_info_list[i]]/(np.median(cell_df['mean_'+ch_info_list[i]].iloc[:5]))
                gr_fluor_norm =  cell_df['gr_'+ch_info_list[i]]
                
                # Search for OM or IM permeabilization using mScarlet or BFP
                inverted_signal = -gr_fluor_norm
                peaks, properties = find_peaks(inverted_signal, height=0.05)
                negative_peak_indices = peaks
                negative_peak_amplitudes = gr_fluor_norm.iloc[negative_peak_indices]
                valid_peaks = []
                for peak_index in peaks:
                    search_region = mean_fluor_norm.iloc[peak_index -1  : np.min([peak_index + 20,len(mean_fluor_norm)])]
                    if (search_region < 0.5).any():
                        valid_peaks = search_region[search_region < 0.5].index.tolist()
                        idx = valid_peaks[0] -cell_df.index.tolist()[0]
                        if len(valid_peaks)>0 and cell_df['time_min_aligned'].iloc[idx]>10:
                            break
                if len(valid_peaks)>0 and cell_df['time_min_aligned'].iloc[idx]>5:
                    idx = valid_peaks[0] -cell_df.index.tolist()[0]
                    if 'pc_diff_IM' in tp_temp.columns and ch_info_list[i] == 'BFP':    # IM perm was detected before OM and I'm detecting BFP loss
                        tp_temp['tp_OM'] = cell_df['time_min'].iloc[idx]                # then I'm detecting OM perm
                        tp_temp['tp_aligned_OM'] = pd.Series(cell_df['time_min_aligned'].iloc[idx])
                        tp_temp['fp_OM'] = pd.Series(cell_df['frame'].iloc[idx])
                    elif not 'pc_diff_IM' in tp_temp.columns:                           # IM perm was not detected 
                        if ch_info_list[i] == 'BFP':                                     # If I'm detecting BFP loss, then this IM perm
                            tp_temp['tp_IM'] = cell_df['time_min'].iloc[idx]
                            tp_temp['tp_aligned_IM'] = pd.Series(cell_df['time_min_aligned'].iloc[idx])
                            tp_temp['fp_IM'] = pd.Series(cell_df['frame'].iloc[idx])
                        elif ch_info_list[i] == 'mScarlet':                              # I'm detecting mScarlet, and don't know about BFP necessarily, then this is OM
                            tp_temp['tp_OM'] = cell_df['time_min'].iloc[idx]
                            tp_temp['tp_aligned_OM'] = pd.Series(cell_df['time_min_aligned'].iloc[idx])
                            tp_temp['fp_OM'] = pd.Series(cell_df['frame'].iloc[idx])
                else: 
                    idx = -1
            
            tp_temp['fr_inj'] = pd.Series(cell_df['fr_inj'].iloc[0])
            if idx>=0 and ch_info_list[i] != 'Sytox':
                tp_temp['tp_'+ch_info_list[i]] = pd.Series(cell_df['time_min'].iloc[idx])
                tp_temp[ch_info_list[i]] = pd.Series(cell_df['mean_'+ch_info_list[i]].iloc[idx])
                tp_temp[ch_info_list[i]+'_norm'] = pd.Series(cell_df['mean_'+ch_info_list[i]].iloc[idx]/np.max(cell_df['mean_'+ch_info_list[i]]))
                tp_temp['idx_'+ch_info_list[i]] = pd.Series(idx)
                tp_temp['fp_'+ch_info_list[i]] = pd.Series(cell_df['frame'].iloc[idx])
                tp_temp['tp_aligned_'+ch_info_list[i]] = pd.Series(cell_df['time_min_aligned'].iloc[idx])
                tp_temp['fr_delay_'+ch_info_list[i]]=  cell_df['frame'].iloc[idx] - cell_df['fr_inj'].iloc[0]
                tp_temp['tp_delay_'+ch_info_list[i]]=  cell_df['time_min'].iloc[idx] - (cell_df['fr_inj'].iloc[0]-2)*fr
                
        # Check false positives, by ensuring that if IM perm was detected by mTagBFP2 loss, then also mScarlet-I was lost   
        if 'tp_IM' in tp_temp.columns: 
            if not 'pc_diff_IM' in tp_temp.columns and not 'tp_OM' in tp_temp.columns:   # if tp_IM was detected not through peri_cyto_diff, but tp_OM was not detected, then that was not a IM perm. event
                tp_temp.drop('tp_IM', axis=1, inplace=True)  
                tp_temp.drop('tp_aligned_IM', axis=1, inplace=True)  
                tp_temp.drop('fp_IM', axis=1, inplace=True)  
    
            elif tp_temp['tp_IM'].iloc[0] < tp_temp['t_gr_stop'].iloc[0]-15 or tp_temp['tp_aligned_IM'].iloc[0]<5:
                tp_temp.drop('tp_IM', axis=1, inplace=True)  
                tp_temp.drop('tp_aligned_IM', axis=1, inplace=True)  
                tp_temp.drop('fp_IM', axis=1, inplace=True)  
            tp_temp['tp_IM_offset'] = tp_temp['tp_IM'].iloc[0]- tp_temp['fr_inj'].iloc[0]*fr
        if 'tp_OM' in tp_temp.columns:
            tp_temp['tp_OM_offset'] = tp_temp['tp_OM'].iloc[0]- tp_temp['fr_inj'].iloc[0]*fr
        tp_df = pd.concat([tp_df,tp_temp])
    return tp_df.reset_index(drop=True)
