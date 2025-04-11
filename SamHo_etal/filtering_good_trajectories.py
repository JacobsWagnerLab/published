# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 14:05:47 2025

@author: Alexandros Papagianakis, Christine Jacobs-Wagner lab
"""



from sklearn.preprocessing import StandardScaler
from sklearn import mixture
import gaussian_statistics as gaus
import data_plotting as datap
import numpy as np
import pandas as pd
import seaborn as sns

def get_long_trajectories(track_df, min_length):
    good_trajectories = []
    length_dict = track_df.groupby('particle_complex').gaussian_amplitude.size().to_dict()
    for traj in length_dict:
        if length_dict[traj]>=min_length:
            good_trajectories.append(traj)
    return good_trajectories



def gaussian_classification(dataframe, feat=['log_msd', 'log_alpha2']):
    """
    The MSD intercept dataframe containing the particle Brownian statistics is the input dataframe.
    
    feat: a list of the features to be used for the classification. The last feature is the particle ID 
        which will be used to link the classes to the particles. The first two variables will be used for the 2D
        gaussian classification
    
    datafrane: pandas dataframe which includes the featured data.
    
    Returns:
        Dictionary which links the particle ID (key) to the class (value)
    """
    data=dataframe.copy()
    features = feat
    data_featured = data.loc[:, features]
    data_featured = data_featured.dropna()
    featured_particles = data_featured['particle_complex'].tolist()
    features = features[0:2]
    data_featured = data_featured.loc[:, features]
    
    x = data_featured.values
    scaled_x = StandardScaler().fit_transform(x)
    # Gaussian mixture classification
    gmm = mixture.GaussianMixture(n_components=2,
                                  covariance_type='tied', max_iter=500, n_init=10, verbose=0)
    gmm.fit(scaled_x)
    gmm.bic(scaled_x)
    prediction = gmm.predict(scaled_x)
    data_featured['mixed_gaussian_class'] = prediction

    data_featured['featured_particles'] = featured_particles
    # gaussian_dict = dict(zip(data_featured.featured_particles, data_featured.mixed_gaussian_class))
    gaussian_dict = pd.Series(data_featured.mixed_gaussian_class.values, 
                              index=data_featured.featured_particles).to_dict()
    
    return gaussian_dict



def filter_long_gaussian_trajectories(track_df, msd_df, title, min_length=150, classification='gaussian'):
    
    long_traj_list =  get_long_trajectories(track_df, min_length)
    long_track_df = track_df[track_df.particle_complex.isin(long_traj_list)]
    alpha_two_dict = gaus.get_alpha2_dictionary(long_track_df, scale_=0.066, 
                                                frame_lag_=1)
    long_msd_df = msd_df[msd_df.particle_complex.isin(long_track_df.particle_complex)]
    msd_intercept_df = long_msd_df[long_msd_df.lag_time_sec == long_msd_df.lag_time_sec.min()]
    msd_intercept_df['alpha2'] = msd_intercept_df.particle_complex.map(alpha_two_dict)
    msd_intercept_df['log_msd'] = np.log10(msd_intercept_df.msd)
    msd_intercept_df['log_alpha2'] = np.log10(msd_intercept_df.alpha2)
    
    if classification == 'gaussian':
        class_dict = gaussian_classification(msd_intercept_df, feat=['log_msd', 'log_alpha2', 'particle_complex'])
        msd_intercept_df['mixed_gaussian_class'] = msd_intercept_df.particle_complex.map(class_dict)
        
        class_zero_particles = msd_intercept_df[msd_intercept_df.mixed_gaussian_class==0].alpha2.mean()
        class_one_particles = msd_intercept_df[msd_intercept_df.mixed_gaussian_class==1].alpha2.mean()
        
        class_verb_dict = {}
        if class_zero_particles < class_one_particles:
            class_verb_dict[0] = 'included'
            class_verb_dict[1] = 'excluded'
        elif class_zero_particles > class_one_particles:
            class_verb_dict[1] = 'included'
            class_verb_dict[0] = 'excluded'
        
        msd_intercept_df['mixed_gaussian_class_verb'] = msd_intercept_df.mixed_gaussian_class.map(class_verb_dict)
    
    elif type(classification) == int or type(classification) == float:
        bad_trajectories = msd_intercept_df[msd_intercept_df.alpha2>classification].particle_complex.to_list()
        msd_intercept_df['mixed_gaussian_class_verb'] = 'included'
        msd_intercept_df.loc[msd_intercept_df.particle_complex.isin(bad_trajectories), 'mixed_gaussian_class_verb'] = 'excluded'
        
    else:
        raise ValueError("The parameter 'classification' should be set to 'gaussian' for an unsupervised 2D Gaussian classification or a number if a hard threshold is applied on the non-Gaussian parameter")
    
    good_df = msd_intercept_df[msd_intercept_df.mixed_gaussian_class_verb=='included']
    bad_df = msd_intercept_df[msd_intercept_df.mixed_gaussian_class_verb=='excluded']

    print(f"{good_df.shape[0]} trajectories with Gaussian behavior and {bad_df.shape[0]} with non-Gaussian behavior, longer than {min_length-1} frames")

    datap.plot_gaussian_classification(msd_intercept_df, title)
    
    alpha_dict = gaus.alpha_dictionary_estimation(msd_df, window_=(0,8))
    good_df['alpha'] = good_df.particle_complex.map(alpha_dict)
    
    return good_df
    
    
