# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 13:56:53 2025

@author: Alexandros Papagianakis, Christine Jacobs-Wagner lab
"""

import numpy as np
import warnings
warnings.filterwarnings('ignore')
import random


def get_alpha2_dictionary(long_track_df, scale_, frame_lag_):
    """
    This function uses the alpha2_estimation function to return a dictionary
    which includes the particle IDs as keys linked to the non-gaussian parameter values.
    
    The track_df (x,y,t positions of trakced particles) and the MSD dataframe of the long trajectories (>=200frames)
    are the inputs.
    """
    
    def alpha2_estimation(dataframe_curated, scale=scale_, frame_lag=frame_lag_):
        """
        non-Gaussian parameter estimation
        https://journals.aps.org/pre/pdf/10.1103/PhysRevE.100.033211
        """
        dataframe_curated = dataframe_curated.sort_values(['particle_complex', 'frame'])
        dataframe_curated['particle'] = dataframe_curated.particle_complex.astype('category').cat.codes

        x_use = 'x'
        y_use = 'y'
            
        dx = (dataframe_curated[x_use].shift(-frame_lag) - dataframe_curated[x_use])*scale
        dy = (dataframe_curated[y_use].shift(-frame_lag) - dataframe_curated[y_use])*scale
        df = dataframe_curated['frame'].shift(-frame_lag) - dataframe_curated['frame']
        dparticle = dataframe_curated['particle'].shift(-frame_lag) - dataframe_curated['particle']
    
        dataframe_curated['dx'] = dx
        dataframe_curated['dy'] = dy
        dataframe_curated['df'] = df
        dataframe_curated['dparticle'] = dparticle
        
        dataframe_curated = dataframe_curated[(dataframe_curated['df']==frame_lag) & (dataframe_curated['dparticle']==0)]
        displacements_sqr = (dataframe_curated['dx']**2 + dataframe_curated['dy']**2).values
        displacements = np.sqrt(displacements_sqr)
        al2 = np.mean(displacements**4)/(2*(np.mean(displacements**2)**2))-1
        
        return al2
    
    alpha2_df = long_track_df.groupby(['particle_complex']).apply(alpha2_estimation)
    alpha2_dict = alpha2_df.to_dict()
    
    return alpha2_dict



def alpha_dictionary_estimation(particle_msd_df, window_):
    
    def alpha_estimation(particle_msd_df, window=window_):
        """
        This function estimates the alpha for an MSD/lag-time window specified as a global variable
        """
        return np.polyfit(np.log10(particle_msd_df.lag_time_sec[window[0]:window[1]]), 
                          np.log10(particle_msd_df.msd[window[0]:window[1]]), 1)[0]
    
    alpha_df = particle_msd_df.groupby(['particle_complex']).apply(alpha_estimation)
    alpha_dict = alpha_df.to_dict()
    
    return alpha_dict



def bootstrap_gaussian_stats(gaussian_df, statistic, sample_size, iterations, random_seed):
    
    sample_dict = {}
    
    stat_list = gaussian_df[statistic].to_list()
    random.seed(random_seed)
    for i in range(iterations):
        sample_dict[i] = random.choices(stat_list, k=sample_size)
    
    return sample_dict



def get_alpha_ratio(sample_list, alpha_threshold):
    
    return np.nonzero(np.array(sample_list)>=alpha_threshold)[0].shape[0]/np.array(sample_list).shape[0]



def collect_alpha_ratios(boots_dict, alpha_threshold):
    
    alpha_ratios = []
    
    for i in boots_dict:
        sample_list = boots_dict[i]
        alpha_ratios.append(get_alpha_ratio(sample_list, alpha_threshold))
        
    return alpha_ratios



def analyze_percentage_of_moving_particles(gaussian_df, random_seed, alpha_threshold):  
    boots_dict = bootstrap_gaussian_stats(gaussian_df, 'alpha', 50, 1000, random_seed)
    return collect_alpha_ratios(boots_dict, alpha_threshold)
        


        