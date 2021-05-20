# -*- coding: utf-8 -*-
"""
Created on Tue May 18 13:47:17 2021

@author: Alexandros Papagiannakis
         Chrisitne Jacobs-Wagner lab, Stanford University
"""

import pandas as pd
import Sanchez_etal_figure_functions as paper_fig # from the JacobsWagnerLab Github:
# https://github.com/JacobsWagnerLab/published/raw/master/Sanchez_2020/Sanchez_etal_figure_functions.py


# load the microtubule position/angle Pandas DataFrames from the JacobsWagnerLab Github
control_path = 'https://github.com/JacobsWagnerLab/published/blob/master/Sanchez_2020/control_df?raw=true'
mutant_10Bgut_path = 'https://github.com/JacobsWagnerLab/published/blob/master/Sanchez_2020/VAB_10Bgut_minus_df?raw=true'
mutant_62gut_path = 'https://github.com/JacobsWagnerLab/published/blob/master/Sanchez_2020/WDR_62gut_minus_df?raw=true'
control_df = pd.read_csv(control_path, compression='zip')
mutant_10Bgut_df = pd.read_csv(mutant_10Bgut_path, compression='zip')
mutant_62gut_df = pd.read_csv(mutant_62gut_path, compression='zip')
# concatenate the data
comet_df = pd.concat([control_df, mutant_10Bgut_df, mutant_62gut_df])


# Download the central lines of the embryos from the JacobsWagnerLab Github
# https://github.com/JacobsWagnerLab/published/blob/master/Sanchez_2020/Central_line_coordinates.zip
# The line_coords_path corresponds to the path of the unzipped downloaded folder
# example
line_coords_path = '/Downloads/Central_line_coordinates'
# uncomment line 78 to visualize all the masked microtubule positions


#--------------- PRE=PROCESSING ------------------#
# drop columns that are not used for the analysis
comet_df = comet_df.drop(columns=['central_line_distance_slope', 'angle_corrected_360'])
# select only the comets that belong into the embryo of interest
embryo_df = comet_df[comet_df.embryo_position==1]

def angle_wrap(angle):
    """
    This function is used to wrap all angles between 0 and 360 degrees
    
    Input:
        angle: float
    Returns:
        The corrected angle
    """
    if angle >=180:
        return angle-180
    elif angle<180:
        return angle
    
def angle_correction(angle):
    """
    This function is used to wrap all angles between 0 and 180 degrees
    
    Input:
        angle: float
    Returns:
        The corrected angle
    """
    if angle < 0:
        return angle+360
    elif angle>360:
        return angle-360
    else:
        return angle
    
# for wrapping the angles around 360 degrees   
embryo_df['circle_angle'] = embryo_df.angle_corrected.apply(angle_correction)    
# for wrapping the angles around 180 degrees
embryo_df['circle_angle_wrap'] = embryo_df.circle_angle.apply(angle_wrap)
# get only the microtubules in the gut
gut_df = embryo_df[embryo_df.gut_position==True]
# get all the embryo/gut masks and the central lines for inspection
# line_coords_path corresponds to the path to the directory DOWNLOADED from the JacobsWagnerLab GitHub repository:
# https://github.com/JacobsWagnerLab/published/blob/master/Sanchez_2020/Central_line_coordinates.zip
#paper_fig.plot_embryo_masks(embryo_df, line_coords_path) # WARNING, download the central line coordinates to visualize the microtubule positions in all the segmented embryos


# ------------------------ FIGURE 4E --------------------------#
paper_fig.plot_angle_distance_wrap(gut_df, mutant='control', min_length=6, plot_mean=False, save_path=None)
paper_fig.plot_angle_distance_wrap(gut_df, mutant='10Bgut', min_length=6, plot_mean=False, save_path=None)
paper_fig.plot_angle_distance_wrap(gut_df, mutant='62gut', min_length=6, plot_mean=False, save_path=None)


# ------------------------ FIGURE S6H --------------------------#
paper_fig.plot_angle_distance_wrap(gut_df, mutant='control', min_length=6, plot_mean=True, save_path=None)
paper_fig.plot_angle_distance_wrap(gut_df, mutant='10Bgut', min_length=6, plot_mean=True, save_path=None)
paper_fig.plot_angle_distance_wrap(gut_df, mutant='62gut', min_length=6, plot_mean=True, save_path=None)


# ------------------------ FIGURE S6G --------------------------#
slopes, distance_list, time_list = paper_fig.plot_trajectories(distance_range=(0,2), min_length=4, gut_df=gut_df, save_path=None)
slopes, distance_list, time_list = paper_fig.plot_trajectories(distance_range=(2,5), min_length=4, gut_df=gut_df, save_path=None)

