# -*- coding: utf-8 -*-
"""
Created on Tue May 18 12:59:41 2021

@author: Alexandros Papagiannakis
         Chrisitne Jacobs-Wagner lab, Stanford University
"""


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pickle

#FUNCTION TO PLOT ALL THE EMBRYO/GUT MASKS AND MIDLINES
def plot_embryo_masks(embryo_df):
    """
    This function is used to plot all the masked microtubule positions and the fitted midlines.
    
    embryo_df: pandas DataFrame including all the microtubule positions inside and outside of the gut for all
               control and mutant embryos.
    """
    central_line_dict = {}
    for exp in embryo_df.experiment.unique():
        if 'control' in exp:
            path = r"N:\Collaborations\Jessica Feldman\ArianaSanchez_microtubule_tracking_in_embryo\Final_analysis\Ariana_Sanchez\Results"
        elif '62gut' in exp:
            path = r"N:\Collaborations\Jessica Feldman\ArianaSanchez_microtubule_tracking_in_embryo\Final_analysis\Ariana_Sanchez\Results_62gut"
        elif '10Bgut' in exp:
            path = r"N:\Collaborations\Jessica Feldman\ArianaSanchez_microtubule_tracking_in_embryo\Final_analysis\Ariana_Sanchez\Results_10Bgut"
        with open(path+'/'+exp+'_central_line_fitted_coordinates', 'rb') as handle:
            fitted_line_coordinates = pickle.load(handle)
        central_line_dict[exp] = fitted_line_coordinates
    # test if the microtubules in and outside of the gut are properly mapped
    for exp in embryo_df.experiment.unique():
        exp_df = embryo_df[embryo_df.experiment==exp]
        en_gut_df = exp_df[exp_df.gut_position==1]
        en_out_gut_df = exp_df[exp_df.gut_position==0]
        
        plt.figure(figsize=(8,8))
        plt.plot(en_gut_df.x, en_gut_df.y, 'o', label='inside the gut')
        plt.plot(en_out_gut_df.x, en_out_gut_df.y, 'o', label='outside the gut')
        plt.plot(central_line_dict[exp][0], central_line_dict[exp][1], color='black', linewidth=2)
        plt.legend()
        plt.xlim(0,512)
        plt.ylim(0,512)
        plt.xlabel('x coordinates (px)')
        plt.ylabel('y coordinates (px)')
        plt.title(exp)
        plt.show()
    

# FUNCTION TO PLOT FIGURE 4E AND FIGURE S6
def plot_angle_distance_wrap(gut_df, mutant, min_length, plot_mean, save_path):
    """
    This function can be used to ploy the 2D histogram of the angles versus the distance from the central line of the embryo
    for the instantaneous comet displacements (Figure 4E).
    It can also be used to plot the average trajectory angle versus its intial distance from the central line (Figure S6H).
    Only the trajectories that start being tracked closer than 5um from the midline are considered.
    
    Input:
        gut_df: panadas DataFrame including the instantaneous comet angles from 0 to 180 degrees and their distance from the central line
                Only the trajectories in the gut are considwered.
        mutant: string - the ID of the control or the mutant embryo ('control', '10Bgut', '62gut')
                'control' for the control embryos
                '10Bgut' for the VAB-10Bgut(-)
                '62gut' or the WDR-62gut(-)'
        min_length: integer - the minimum length of the comet trajectories considered. The length corresponds to the number of frames
        plot_mean: bool - True to plot the mean angles of the trajectoris on top of the heatmap, as in Figure S6H
                          False to plot the heatmap of the instantaneous displacements as in figure 4E
        save_path: string - the path to the directory where the figures will be saved
                   Choose None to not save the figure 
        
    Returns:
        plots and saves the specofied 2D histograms
    """
    control_df = gut_df[gut_df.exp_type.str.contains(mutant)]
    long_control_dict =  control_df.groupby('unique_id').size().to_dict()
    control_df['trajectory_length'] = control_df.unique_id.map(long_control_dict)
    long_control_df = control_df[control_df.trajectory_length>=min_length]

    def angle_estimation(dataframe):
        """
        This function can be used to estimate the mean slope of a comet trajectory and convert a slope into angle degrees
        The mean slope corresponds to the slope of a linear regression fitted to the trajectory x,y positions.
        
        Input:
            dataframe: pandas DataFrame - a dataframe including the x and y coordinates of the trajectory positions in two separate columns "x" and "y"
        
        Returns:
            angle: the mean angle of the comet trajectory in degrees
        """
        slope =  np.polyfit(dataframe['x'], dataframe['y'], 1)[0]
        
        angle = np.rad2deg(np.arctan(slope))
        
        return angle
    
    def angle_correction(angle):
        """
        Wraps the angle between 0 and 180 degrees.
        
        Input:
            angle: float 
        Returns:
            The corrected angle
        """
        if angle < 0:
            return angle+180
        elif angle>180:
            return angle-180
        else:
            return angle
    
    def initial_distance(dataframe):
        """
        Returns the initial 1D potition of a comet trajectory.
        The 1D position corresponds to the minimum distance from the central line.
        """
        dataframe = dataframe[dataframe.frame==dataframe.frame.min()]
        
        return dataframe.central_distance_um.values[0]
    
    def initial_central_angle(dataframe):
        """
        Returns the angle of the central line at the initial position of the comet trajectory.
        This angle can be used to correct the mean angle of the trajectory.
        """
        dataframe = dataframe[dataframe.frame==dataframe.frame.min()]
        
        return dataframe.central_angle.values[0]
        

    mean_angle = long_control_df.groupby(['unique_id']).apply(angle_estimation)
    mean_distance = long_control_df.groupby(['unique_id']).mean().central_distance_um
    mean_central_angle = long_control_df.groupby(['unique_id']).mean().central_angle
    initial_distance = long_control_df.groupby(['unique_id']).apply(initial_distance)
    initial_angle = long_control_df.groupby(['unique_id']).apply(initial_central_angle)
    
    mean_angle_dict = pd.Series(mean_angle.values,index=mean_angle.index).to_dict()
    mean_distance_dict = pd.Series(mean_distance.values,index=mean_distance.index).to_dict()
    central_angle_dict = pd.Series(mean_central_angle.values,index=mean_central_angle.index).to_dict()
    initial_distance_dict = pd.Series(initial_distance.values,index=initial_distance.index).to_dict()
    initial_central_angle_dict = pd.Series(initial_angle.values,index=initial_angle.index).to_dict()
    
    mean_df = pd.DataFrame()
    mean_df['unique_id'] = list(mean_angle_dict.keys())
    mean_df['mean_comet_angle'] = mean_df['unique_id'].map(mean_angle_dict)
    mean_df['mean_central_angle'] = mean_df['unique_id'].map(central_angle_dict)
    mean_df['initial_central_angle'] = mean_df['unique_id'].map(initial_central_angle_dict)
    mean_df['mean_distance_um'] = mean_df['unique_id'].map(mean_distance_dict)
    mean_df['mean_angle_corrected'] = mean_df['mean_comet_angle'] - mean_df['initial_central_angle']
    mean_df['mean_angle_corrected_180'] = mean_df['mean_angle_corrected'].apply(angle_correction)
    mean_df['initial_distance_um'] = mean_df['unique_id'].map(initial_distance_dict)
    
    mean_df = mean_df[mean_df.initial_distance_um<=5] # select the range of trajectories from the midline given their initial position
    long_control_df = long_control_df[long_control_df.unique_id.isin(mean_df['unique_id'].tolist())]

    print('displacements:', long_control_df.shape[0])
    print('long trajectories:', mean_df.shape[0])
    
    if plot_mean == True:
        plt.figure(figsize=(5,4))
        plt.plot(mean_df.initial_distance_um, 
                 mean_df.mean_angle_corrected_180, 'o', markersize=4, color='red')
        plt.xlabel('Distance from annotated midline (um)')
        plt.ylabel('mean trajectory angle (degrees)')
        plt.title(mutant)
    elif plot_mean == False:
        plt.figure(figsize=(8,5))
        sns.kdeplot(
                data=long_control_df, x="central_distance_um", y="circle_angle_wrap",
                fill=True, thresh=0, bw_adjust=1, levels=30, cmap="mako", cbar=True
                )
        plt.title(mutant)
    
    plt.xlim(-0.5,6)
    plt.ylim(0,180)
    
    if plot_mean==True:
        if save_path != None:
            plt.savefig(save_path + '/' + mutant+'_'+'angle_distance_with_mean.pdf', transparent=True)
    elif plot_mean==False:
        if save_path != None:
            plt.savefig(save_path + '/' + mutant+'_'+'angle_distance_without_mean.pdf', transparent=True)
    plt.show()


# FUNCTION TO PLOT FIGURE S6G
def plot_trajectories(distance_range, min_length, gut_df, save_path):
    """
    This function plots the 1D trajectories at different segmentes from the central line.
    The 1D distance from the central line of each individual trajectory is plotted over time.
    The time and 1D distance for each trajectory is set to zero for the first frame.
    Thus all individual trajectories start from the origin.
    
    Input:
        distance_range: tuple of 2 floats or integers indicating the min and max bound of the segment 
                        where the comet trajectories begin. The min and max values correspond to the distance
                        from the central line of the embryo.
                        To visualize trajectories that start at 1 to 5um from the central line, the tuple (1,5)
                        should be used.
        min_length: integer - the minimum length of the plotted trajectories in frames
        gut_df: pandas DataFrame - that includes the 1D position of the comet relative to the central line and time
                In our analysis the 1D position is specified as "central_distance_um" and the "frame" column which includes the time frames
                of the experiment is used to get the time of the comet trajectory. The frame column is multiplied by the frame interval in seconds (0.1sec)
        save_path: string - the path to the directory where the figures will be saved
                   Choose None to not save the figure 
                   
    Returns:
        Plots the 1D trajectories over time
        slopes: a dictionary of lists:
            key: string - the ID of the experiment (e.g. "control" or "mutant")
                list - a list of the slopes of each indovidual trajectory. Each slope corresponds to the slope of a linear regression fitted to the 1D trajectory.
        distance_list: a dictionary of lists:
            key: string - the ID of the experiment (e.g. "control" or "mutant")
                list - a list of the 1D distance each comet travels. The 1D distance corresponds to the 1D position at the end of the comet trajectory minus the 1D position at the start
        time: a dictionary of lists:
            key: string - the ID of the experiment (e.g. "control" or "mutant")
                list - a list of the time (longevity) of each comet trajectory. The time corresponds to the frame at the end of the trajectory minus the frame at the start multiplied by the frame interval in seconds        
    """
    slopes = {}
    distance_list = {}
    time_list = {}
    for exp_tp in gut_df.exp_type.unique():
        exp_df = gut_df[gut_df.exp_type==exp_tp]
    
        slopes[exp_tp] = []
        distance_list[exp_tp] = []
        time_list[exp_tp] = []
        plt.figure(figsize=(7,5))
        for mt in exp_df.unique_id.unique():
            mt_df = exp_df[exp_df.unique_id==mt]
            if mt_df.shape[0] >= min_length:
                if  mt_df.central_distance_um.values[0]>=distance_range[0] and mt_df.central_distance_um.values[0]<distance_range[1]:
                    plt.plot((mt_df.frame-mt_df.frame.values[0])*0.1, mt_df.central_distance_um-mt_df.central_distance_um.values[0], color='gray', alpha=0.2)
                    slopes[exp_tp].append(np.polyfit((mt_df.frame-mt_df.frame.values[0])*0.1, mt_df.central_distance_um-mt_df.central_distance_um.values[0], 1)[0])
                    distance_list[exp_tp].append(mt_df.central_distance_um.values[-1] - mt_df.central_distance_um.values[0])
                    time_list[exp_tp].append((mt_df.frame.values[-1]-mt_df.frame.values[0])*0.1)
                
        plt.xlabel('time (sec)')
        plt.ylabel('distance travelled from the central line (μm)')
        plt.ylim(-1,4)
        plt.xlim(0,8)
        title_string = exp_tp+', distance range:'+str(distance_range[0])+' to '+str(distance_range[1])+'um'
        plt.title(title_string)
        if save_path != None:
            plt.savefig(save_path+'/'+'kymograph_'+exp_tp+'_from'+str(distance_range[0])+'_to_'+str(distance_range[1])+'_um_.pdf')
        plt.show()
    
    return slopes, distance_list, time_list




