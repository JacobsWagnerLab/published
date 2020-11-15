#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 23:08:48 2020

@author: Alexandros Papagiannakis, Christine-Jacobs Wagner lab, Stanford University, 2020
"""

from skimage import measure
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cmath



# The analysis pipeline
# Run the class file
# If the class scipt and this script are in the same folder, the class script can be imported here

image_path = 'insert the path with the stream acquisition images'
snapshots_path = 'Insert the path with the embryo mid-line or gut markers"
save_path = 'The path of the folder where the dataframes will be saved'
experiment = 'The experiment ID'

# initializing the class
microtubule = microtub_tracking(image_path, snapshots_path, save_path, experiment)

# If you want to check the segmentation parameters
#microtubule.check_segmentation_parameters(0, 'TEST 488')

# Some parameters to try
#{'log': (2.0, 40, 98.5),
# 'adaptive_threshold': (1.0, 7, -45.0),
# 'erosion': 0,
# 'max_spot_area': 40,
# 'min_spot_area': 5,
# 'aspect_ratio': 1.0,
# 'log_zoom': (2.0, 30, 90.0)}


# get the mask of the embryo gut
# Parameters to run:
#   The channel of the embryo gut marker (e.g. 'TEST 405') as specified in the .nd2 file
#   The parameters for the adaptive filter (2, 71, -20)
#   The iterations of the binary closing process (9)
#   The rounds of dilation (0)
#   The minimum number of pixels of the gut mask (7000)
gut_mask = microtubule.get_gut_mask('TEST 405', (2,71,-20), 9, 0, 7000)

# get the mask of the embryo
# Parameters to run:
#   The channel of the embryo marker (e.g. 'TEST 488') as specified in the .nd2 file
#   The parameters for the adaptive filter (2, 501, -1)
#   The iterations of the binary closing process (8)
#   The rounds of dilation (0)
#   The minimum number of pixels of the gut mask (5000)
embryo_mask = microtubule.get_embryo_mask('TEST 488', (2,501,-1), 8, 0, 5000)

# get the central of the embryo gut
# Need to use an appropriate matplotlib backend and not inline plotting.
# The Qt5 backend is recommended > matplotlib.use('Qt5Agg') to plot the midline image in a separate window and select coordinates.
# Parameters to run:
#   The channel of the mideline marker (e.g. 'TEST 488') as specified in the .nd2 file
    # if 'MAX' is selected, the max projection of the stram acquisition images is used
#   The order of the fitted univariate smoothing spline
#   The ammount of smoothing
#   The LUT of the image used to draw the central line
#   If the embryo is positioned vertically, set rotate to True
central_line = microtubule.get_central_line('MAX', order=1, smoothing=1000, lut = (0,2500), rotate=False)

# run the microtubule segmentation
microtubule.run_segmentation(channel='TEST 488', log_parameters=(2, 40, 98.5), 
                             adaptive_threshold_parameters=(1, 7, -45), 
                             erosion_iterations=0, min_spot_size=5, max_spot_size=40, spot_aspect_ratio=1, 
                             log_parameters_zoom=(2,30,90))

# If you want to check the segmentation
# frames: the range of frames to check
# channel: the channel on which the segmented spots are shown
# save: the path where the segmentation images are saved in .jpeg format
#microtubule.show_segmentation(channel='TEST 488', frames=(10,11), save='none')

# tracking the segmented microtubules
# Parameters:
    # The search radius in pixels (4)
    # The memory of linkeage (2)
microtubule.microtubule_tracking(4,2)

# If you want to check the tracking
# Parameter:
#   The minimum length of the trajecotires to be displayed
microtubule.show_microtubule_trajectories(5)


# get the angle and the speed of the microtubules
# WARNING: need to run the get_embryo_mask and get_gut_mask and get_central_line functions before running this one
# Parameter:
#   Time-lag: the time lag for smoothing the comet trajectories (window of smoothing average)
angle_df = microtubule.get_microtubule_speed_and_angles(time_lag=4)

# get the microtubule dataframes
#microtubule_df = pd.read_pickle(microtubule.save_path+'/'+microtubule.experiment+'_segmented_microtubules_df', compression='zip')
#microtubule_df = pd.read_pickle(microtubule.save_path+'/'+microtubule.experiment+'_tracked_microtubules_df', compression='zip')
#microtubule_df = pd.read_pickle(microtubule.save_path+'/'+microtubule.experiment+'_angle_speed_df', compression='zip')

# Get the microtubule dataframe and all the microtubules in the embryo
microtubule_df = microtubule_df[microtubule_df['embryo_position']==True]
# In the gut
gut_df = microtubule_df[microtubule_df ['gut_position']==True]
# Outside of the gut
non_gut_df = microtubule_df[microtubule_df ['gut_position']==False]


# get a plot with the lineages
# jet or rainbow are some good colormaps
def plot_trajectory_lines(save_path, experiment, length, scalebar, scalebar_position, lut, colormap_input=plt.cm.rainbow):
    """
    This function plots the trajectory lines in the C. elegans embryo.
    Each trajectory line is color-coded according to the frames which it consists of in the stream acquisition.
    The colormap ranges from 0 sec to the maximum duration of the stream acquisition.
    An interval of 100msec is assumed. 
    The controur of the embruo gut is also plotted.
    
    The figure is saved in the specified folder (save_path) in pdf format (vector graph)
    
    Input:
        save_path: string - the path where the analysis dataframes are saved and where the generated figure 
        experiment: string - the id of the experiment (must match the experiment id of the analysis)
        length: integer - the minimum length of the trajectories
        scalebar: float or integer - the size of the scalebar in μm
        scalebar_position: tuple of integers - the x, y pixel coordinates of the left edge of the scalebar
        lut: tuple of integers - the min and the max pixel values of the LUTs used to set the contrast and brightness of the image.
        colormap_input: matplotlib.pytplot.cm object - the colormap to be used
    """
    microtubule_df = pd.read_pickle(microtubule.save_path+'/'+microtubule.experiment+'_angle_speed_df', compression='zip')
    
#    with open(save_path+'/'+experiment+'_central_line_fitted_coordinates', 'rb') as handle:
#        central_line_coordinates = pickle.load(handle)
    
    with open(save_path+'/'+experiment+'_gut_mask', 'rb') as handle:
        gut_mask = pickle.load(handle)
    
    with open(save_path+'/'+experiment+'_embryo_mask', 'rb') as handle:
        embryo_mask = pickle.load(handle)
    
    gut_contours = measure.find_contours(gut_mask, 0.5)[0]
    embryo_contours = measure.find_contours(embryo_mask, 0.5)[0]
    
    gut_y, gut_x = map(list, zip(*gut_contours))
    embryo_y, embryo_x = map(list, zip(*embryo_contours))
    
    
    image = microtubule.max_projection()
    
#    plt.figure(figsize=(24,20))
    plt.imshow(microtubule.image_arrays['TEST 488'][0], cmap='gray')
    plt.imshow(image, cmap='gray', vmin=lut[0], vmax=lut[1])
    for traj in microtubule_df['microtubule_id'].unique():
        traj_df = microtubule_df[microtubule_df['microtubule_id']==traj]
        print(traj)
        if traj_df.shape[0] >= length:
            
            
            frames = np.array(traj_df['frame'])
            x = np.array(traj_df['x'])
            y = np.array(traj_df['y'])
            
            for i in range(len(frames)-1):
                color_map = colormap_input(i/(len(frames)-1))
                plt.plot(x[i:i+2], y[i:i+2], color=color_map, linewidth=1, solid_capstyle='round')
            
            
    scaled_map = plt.cm.ScalarMappable(cmap=colormap_input, norm=plt.Normalize(vmin=microtubule_df['time_sec'].min(), vmax=microtubule_df['time_sec'].max()))
    plt.colorbar(scaled_map)
    plt.plot(gut_x, gut_y, '--', color='white', linewidth=0.5)
#    plt.plot(central_line_coordinates[0], central_line_coordinates[1], color='black', linewidth=1.5)
#    plt.xlim(0,512)
#    plt.ylim(0,512)
    # 0.11 μm per pixel or 10 μm per 90.9 pixels.
    
    plt.plot((scalebar_position[0],scalebar_position[0]+(10/microtubule.scale)), (scalebar_position[1], scalebar_position[1]), color='white', linewidth='3')
    plt.savefig(save_path+'/'+experiment+'_trajectory_lines.pdf')
    plt.show()
    
# Run this line to plot and save the trajectory lines
plot_trajectory_lines(save_path, experiment, length=8, scalebar=10, scalebar_position = (400,450), lut = (0,5000), colormap_input=plt.cm.rainbow)





