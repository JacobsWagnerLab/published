# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 01:04:01 2021

@author: Alexandros Papagiannakis, Christine Jacobs-Wagner lab, Sarafan ChEM-H, Stanford University, 2023

This script includes a collection of functions that were used for the analysis of the microfluidics data.
The functions are organized in the following fields:
    1. LOADING IMAGES, LABELS, MEDIAL AXES
    2. BACKGROUND CORRECTION FOR MICRODLUIDICS
    3. NUCLEOID SEGMENTATION AND POSITION 
    4. COMPLETE TRAJECTORIES AND CELL DIVISION CYCLE STATS 
    5. MEDIAL AXIS CALCULATION
    6. CELL POLARITY ESTIMATION BASED ON THE MEDIAL AXIS AND THE CELL DIVISION CYCLE 
    7. RELATIVE PIXEL COORDINATES 
    8. CELL AND DETECTED OBJECT (e.g. NUCLEOID) FLUROESCENCE STATISTICS
    9. NUCLEOID TRACKING
    10. CALCULATION OF CELL DIVISION CYCLE AND INSTANTANEOUS RATES (including growth rate estimations)
    11. DEFINITION OF CELL AGES
"""

import pandas as pd
import pickle
import matplotlib.pyplot as plt
import numpy as np
import os
from skimage import filters, io, measure, morphology  # transform
from scipy import ndimage, interpolate
from PIL import Image
import multiprocessing as mp
import warnings




# ------------------   1. LOADING IMAGES, LABELS, MEDIAL AXES  --------------------- #
def get_long_trajectories(cell_dataframe, length_range):
    """
    This function is used to collect the cell trajectories IDs that are larger and smaller 
    than the specified trajectory length range in frames
    """
    
    size_df = cell_dataframe.groupby('cell_trajectory_id').size()
    long_cc_list = size_df[(size_df.values>=length_range[0])&(size_df.values<=length_range[1])].index.tolist()
    
    return long_cc_list


def get_cell_labels(results_path, curated=True):
    """
    This function returns all cell labels.
    """
    cell_labels = {}

    dir_list = os.listdir(results_path)
    if curated == True:
        label_dir = [s for s in dir_list if 'curated_watershed_dict' in s]
    elif curated == False:
        label_dir = [s for s in dir_list if 'watershed_labels_dict' in s]
    for s in label_dir:
        experiment = s[0:s.find('pos')-1]
        try:
            position = int(s[s.find('pos')+3:s.find('pos')+5])
        except ValueError:
            position = int(s[s.find('pos')+3:s.find('pos')+4])
        print(experiment, position)
        
        if experiment not in cell_labels:
            cell_labels[experiment] = {}
        
        cell_labels[experiment][position]=pickle.load(open(results_path+'/'+s, 'rb'))
    
    return cell_labels


def get_specific_cell_labels(results_path, experiment, position, curated=True):
    """
    This function returns the cell labels of a specific microfluidic channel position.
    """
    dir_list = os.listdir(results_path)
    if curated == True:
        label_dir = [s for s in dir_list if 'curated_watershed_dict' in s]
    elif curated == False:
        label_dir = [s for s in dir_list if 'watershed_labels_dict' in s]
    exp_dir = [s for s in label_dir if experiment in s]
    pos_dir = [s for s in exp_dir if 'pos'+str(position)+'_' in s]

    if len(pos_dir)==0:
        print('There is no cell label definition for this experiment or position')
        cell_labels=np.nan
    elif len(pos_dir)==1:
        print('Loading cell labels:',pos_dir[0])
        cell_labels = {}
        cell_labels[experiment] = {}
        cell_labels[experiment][position] = pickle.load(open(results_path+'/'+pos_dir[0], 'rb'))
    elif len(pos_dir)>1:
        print('Multiple cell label definitions found', pos_dir)
        cell_labels=np.nan
    
    return cell_labels


def get_fluorescence_image(images_path, xy_position, channel, position, frame, edge_width):
    """
    This function returns one specific background corrected fluorescence image, 
    for a specific channel, microfluidic channel position and frame.
    
    If background correction is required set edge_width to an integer above zero 
    """
    experiment_folders = os.listdir(images_path)
    cropped_path = images_path+'/'+[s for s in experiment_folders if xy_position+'_cropped' in s][0]+'/-'+str(position)+'/c'+str(channel)+'Exp'
    cropped_images = os.listdir(cropped_path)
    time_stamp_length = cropped_images[0].find('c'+str(channel)) - (cropped_images[0].find(xy_position)+len(xy_position)+1)
    frame_stamp_length = len(str(frame+1))
    frame_stamp = (time_stamp_length-frame_stamp_length)*'0'+str(frame+1)
    read_image = [s for s in cropped_images if xy_position+'t'+frame_stamp in s]
    if len(read_image)==1:
        fluor_image = io.imread(cropped_path+'/'+read_image[0])
        if edge_width > 0:
            fluor_image = background_correction(fluor_image, edge_width)
        return fluor_image
    elif len(read_image)==0:
        return np.nan
    

    
def get_fluorescence_images(images_path, channel, frame_range, edge_width):
    """
    This function returns all background corrected fluorescence images along a range of time-points and a specific channel. 
    
    If background correction is required set edge_width to an integer above zero 
    """
    fluorescence_images_dict = {}
    experiment_folders = os.listdir(images_path)
    cropped_images = [s for s in experiment_folders if 'cropped' in s]
    for exp_img in cropped_images:
        xy_position = exp_img[exp_img.find('xy'):exp_img.find('xy')+3]
        fluorescence_images_dict[xy_position] = {}
        cropped_images_2 = os.listdir(images_path+'/'+exp_img)
        cropped_images_2 = [s for s in cropped_images_2 if 'c' not in s]
        for pos in cropped_images_2:
#            print(pos)
            if pos != '.DS_store' and pos[0]=='-':
                position = int(pos[1:])
                channel_string = 'c'+str(channel)+'Exp'
                print(xy_position, position, channel_string)
                fluorescence_images_dict[xy_position][position]={}
                for frame in range(frame_range[0], frame_range[1]):
                    fluorescence_images_dict[xy_position][position][frame]=get_fluorescence_image(images_path, xy_position, channel, position, frame, edge_width)
    return fluorescence_images_dict
            

def get_channel_fluorescence_images(images_path, xy_position, channel, position, frame_range, edge_width):
    """
    This function returns background corrected fluorescence images for a specific channel, xy position,
    microfluidic channel and  range of time-points. 
    
    If background correction is required set edge_width to an integer above zero 
    """
    fluorescence_images_dict = {}
    for frame in range(frame_range[0], frame_range[1]):
        fluorescence_images_dict[frame]=get_fluorescence_image(images_path, xy_position, channel, position, frame, edge_width)
    return fluorescence_images_dict
    
    
# ------------------   2. BACKGROUND CORRECTION FOR MICRODLUIDICS  --------------------- #
def background_correction(fluor_image, edge_width):
    """
    Used in the get_fluorescence images function
    """
    return fluor_image - np.mean(np.concatenate((fluor_image[:, 0:edge_width],fluor_image[:, -edge_width:]),axis=1), axis=1)[:,None]



# ------------------   3. NUCLEOID SEGMENTATION AND POSITION  --------------------- #
def log_adaptive_filter(image, parameters): 
    """
    This fucntion constructs an LoG filer (Laplace of Gaussian) as well as an adaptive filter to segment particles
    
    Parameters
    ---------
    image: numpy array - the image to be filtered and thresholded (usually this is the background subtracted image)
    particle_detection_parameters: list - the parameters for the LoG/adaptive filter
        [0] Smoothing factor for the Gaussian filter (https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.gaussian)
        [1] Laplace threshold (https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.laplace)
        [2] Hard threshold on the LoG image as a percentile of the brightest pixels
        [3] Gaussian smoothing factor before applying the adaptive threshold (https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.gaussian)
        [4] Block_size for the adaptive threshold (https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.threshold_local)
        [5] Offset for the adaptive threshold (https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.threshold_local)
        [6] Mask erosion applied after the adaptive threshold (it can be set to zero). Otherwise it is a positive integer (erosion rounds).
    
    Returns
    -------
    [0] The thresholded image after the LoG filter
    [1] The image after the adaptive threshold
    [2] The filter after multiplying the thresholded smoothed image with the adaptively thresolded image
    """
    
#    image[image<0]=0
    # LoG filter with hard threshold
    image_smoothed = ndimage.gaussian_filter(image, sigma = parameters[0])
    image_laplace = filters.laplace(image_smoothed, parameters[1])
    image_pixel_intensities = image_laplace.ravel()
    sorted_pixel_intensities = np.sort(image_pixel_intensities)
    pixel_intensity_threshold = sorted_pixel_intensities[int((parameters[2]/100)*len(sorted_pixel_intensities))]
    log_image = image_laplace > pixel_intensity_threshold
    # Adaptive smoothing
    image_smoothed_2 = ndimage.gaussian_filter(image, sigma=parameters[3])
    adaptive_threshold = filters.threshold_local(image_smoothed_2, block_size=parameters[4], offset=parameters[5])
    adaptively_masked_image = image_smoothed_2 > adaptive_threshold
    masked_image = adaptively_masked_image * log_image
    if parameters[6] > 0:
        adaptively_masked_image = ndimage.morphology.binary_erosion(masked_image, iterations = parameters[6])
        final_image = adaptively_masked_image
    elif parameters[6] == 0:
        final_image = masked_image
    
    return log_image, adaptively_masked_image, final_image


def segment_nucleoid(nucleoid_image, log_adaptive_params, min_mask_size, hole_fill):
    """
    This function uses an LoG/Adaptive filter to define the nucleoid boundaries and compute the nucleoid mask.
    
    Parameters
    ----------
    log_adaptive_params - list of 7 parameters:
        [0] Smoothing factor for the Gaussian filter (https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.gaussian)
        [1] Laplace threshold (https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.laplace)
        [2] Hard threshold on the LoG image as a percentile of the brightest pixels
        [3] Gaussian smoothing factor before applying the adaptive threshold (https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.gaussian)
        [4] Block_size for the adaptive threshold (https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.threshold_local)
        [5] Offset for the adaptive threshold (https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.threshold_local)
        [6] Mask erosion applied after the adaptive threshold (it can be set to zero). Otherwise it is a positive integer (erosion rounds).
    
    min_mask_size: integer (pixels) - the minimum mask sizes considered (set to 10)
    hole_fil: integer (pixels) - the max size of the holes to be filed in the mask (set to 3)
    
    Returns
    -------
    The background corrected nucleoid image
    The nucleoid masks
    """
    #log_adaptive_params=[2,1000,80,2,9,-2,0])
    hu_mask = log_adaptive_filter(nucleoid_image, parameters=log_adaptive_params)[2]
    hu_mask = morphology.remove_small_objects(hu_mask, min_size=min_mask_size)
    hu_mask = morphology.remove_small_holes(hu_mask, area_threshold=hole_fill)
    
    return hu_mask


def show_nucleoid_movie(nucleoid_images_dict, log_adaptive_params, min_area, hole_fill, save_path):
    """
    This function can be used to inspect the nucleoid segmentation.
    """
    for img_index in nucleoid_images_dict:
        if type(nucleoid_images_dict[img_index])!=float:
            fig, (ax1, ax2)  = plt.subplots(1,2, figsize=(5,13))
            ax1.imshow(nucleoid_images_dict[img_index])
            ax1.text(10,10,'0.066um/px', fontsize=12, color='white')
            ax1.axes.xaxis.set_ticks([])
            ax2.imshow(segment_nucleoid(nucleoid_images_dict[img_index], log_adaptive_params, min_area, hole_fill))
            ax2.axes.yaxis.set_ticks([])
            ax2.axes.xaxis.set_ticks([])
            ax2.text(10,10,'Time: '+str(img_index)+'min', fontsize=12, color='white')
            if os.path.isdir(save_path):
                plt.savefig(save_path+'/'+str(img_index)+'.jpeg')
            plt.show()
        

def segment_all_nucleoids(images_path, channel, cell_dataframe,
                          log_adaptive_params, min_mask_size, hole_fill, edge_width,
                          save_path):
    """
    This function can be used to segment all nucleoids.
    See example in the Jupyter notebook: Microfluidics segmentation and tracking
    """
    print('segmenting the nucleoids')
    hu_masks_dict = {}
    for exp in cell_dataframe.experiment.unique():
        xy_pos = exp[exp.find('xy'):] 
        print(exp)
        xy_pos_df = cell_dataframe[cell_dataframe.experiment.str.contains(xy_pos)]
        hu_masks_dict[xy_pos]= {}
        for pos in xy_pos_df.position_int.unique():
            pos_df = xy_pos_df[xy_pos_df.position_int==pos]
            hu_masks_dict[xy_pos][pos] = {}
            print(xy_pos, pos)
            for fr in range(pos_df.frame.min(), pos_df.frame.max()+1):
                try:
                    nucleoid_image = get_fluorescence_image(images_path, xy_pos, channel, pos, fr, edge_width)
                except IndexError:
                    xy_pos_2 = 'xy'+str(int(xy_pos[2:]))
                    nucleoid_image = get_fluorescence_image(images_path, xy_pos_2, channel, pos, fr, edge_width)
                if type(nucleoid_image)!=float:
                    hu_masks_dict[xy_pos][pos][fr] = segment_nucleoid(nucleoid_image, log_adaptive_params, min_mask_size, hole_fill)
                else:
                    hu_masks_dict[xy_pos][pos][fr] = np.nan
    
    pickle.dump(hu_masks_dict, open(save_path+'/nucleoid_masks_dict', 'wb'))
    
    return hu_masks_dict


# ------------------   4. COMPLETE TRAJECTORIES AND CELL DIVISION CYCLE STATS  --------------------- #
def get_mother_cells(cell_dataframe, results_path, area_percentage=(1.65, 2.35)):
    """
    This function returns the mother cell trajectory IDs of all lineages.
    """
    cell_labels = get_cell_labels(results_path, curated=True)
    
    mother_cell_dict = {}
    for traj in cell_dataframe.cell_trajectory_id.unique():
        
        traj_df = cell_dataframe[cell_dataframe.cell_trajectory_id==traj]
        traj_df = traj_df.sort_values('frame')
        experiment = traj_df.experiment.values[0]
        position = traj_df.position_int.values[0]
        min_frame = traj_df.frame.min()
        min_area = traj_df.area_px.values[0]
        min_x = traj_df.x.values[0]
        min_y = traj_df.y.values[0]
        #print(position, min_frame, min_area)
        
        if min_frame == 0:
            mother_cell_dict[traj] = np.nan
        elif min_frame > 0:
            previous_labels = cell_labels[experiment][position][min_frame-1]
            mother_traj_label = previous_labels[int(min_y)][int(min_x)]
            if mother_traj_label > 0:
                mother_area_px = previous_labels[previous_labels==mother_traj_label].shape[0]
                if mother_area_px > area_percentage[0]*min_area and mother_area_px < area_percentage[1]*min_area:
                    try:
                        mother_cell_dict[traj] = cell_dataframe[(cell_dataframe.position_int==position)&(cell_dataframe.experiment==experiment)&
                                        (cell_dataframe.frame==(min_frame-1))&
                                        (cell_dataframe.label==mother_traj_label)].cell_id.values[0]
                    except IndexError:
                        mother_cell_dict[traj] = np.nan

    cell_dataframe['mother_cell'] = cell_dataframe.cell_trajectory_id.map(mother_cell_dict)
    
    return cell_dataframe, mother_cell_dict


def get_mothers_with_mothers(cell_dataframe):
    """
    This function is used to get complete cell division cycle trajectories:
        mother cells that have mothers and divide to daughters.
    """
    list_of_mothers = cell_dataframe.mother_cell.dropna().unique().tolist()
    cell_dataframe = cell_dataframe[~cell_dataframe.mother_cell.isnull()]
    
    mothers_with_mothers = []

    for traj in cell_dataframe.cell_trajectory_id.unique():
        
        traj_df = cell_dataframe[cell_dataframe.cell_trajectory_id==traj]
        traj_df = traj_df.sort_values('frame')
        if traj_df.cell_id.values[-1] in list_of_mothers:
            #print(traj)
            mothers_with_mothers.append(traj)
        
    return mothers_with_mothers


def get_sister_cell_pairs(cell_dataframe):
    """
    Returns a dictionary with the mother cell as key and a list of the sister cells as value.
    The list may contain one or two sister cells but no more.
    
    This function was not used eventually.
    """
    
    # the mother cells need to be recognized first
    # run get_mother_cells function
    return cell_dataframe.groupby('mother_cell').cell_trajectory_id.unique().to_dict()


def get_cell_cycle_phases(cell_dataframe):
    """
    Returns the cell division cycle phase for each trajectory,
    with 0 being birth and 1 division.
    """
    cell_dataframe['min_trajectory_frame'] = cell_dataframe.cell_trajectory_id.map(cell_dataframe.groupby('cell_trajectory_id').frame.min().to_dict())
    cell_dataframe['max_trajectory_frame'] = cell_dataframe.cell_trajectory_id.map(cell_dataframe.groupby('cell_trajectory_id').frame.max().to_dict())
    cell_dataframe['zeroed_frame'] = cell_dataframe.frame - cell_dataframe.min_trajectory_frame
    cell_dataframe['max_zeroed_frame'] = cell_dataframe.cell_trajectory_id.map(cell_dataframe.groupby('cell_trajectory_id').zeroed_frame.max().to_dict())
    cell_dataframe['cc_phase'] = cell_dataframe.zeroed_frame / cell_dataframe.max_zeroed_frame
    
    return cell_dataframe




# ------------------   5. MEDIAL AXIS CALCULATION  --------------------- #
def crop_cell_mask(cell_mask, crop_frame=3):
    """
    Crops the cell mask
    """
    max_y = np.nonzero(cell_mask)[0].max() + crop_frame + 1
    min_y = np.nonzero(cell_mask)[0].min() - crop_frame
    max_x = np.nonzero(cell_mask)[1].max() + crop_frame + 1
    min_x = np.nonzero(cell_mask)[1].min() - crop_frame
    
    cropped_cell_mask = cell_mask[min_y:max_y, min_x:max_x]
    
    return cropped_cell_mask, (min_y, max_y, min_x, max_x)


def get_medial_axis(cell_mask, crop_frame, 
                    resize_factor, order, 
                    smoothing, end_cap, show):
    """
    This function construct the medial axis of a signle cell, 
    as well as the relative coordinates of the cell from one pole to the other.
    """    

    def resize_rotate_image(cropped_cell_mask, resize_factor):
        
        resized_cell_mask = np.array(Image.fromarray\
                                     (cropped_cell_mask).resize((cropped_cell_mask.shape[1]*resize_factor, \
                                     cropped_cell_mask.shape[0]*resize_factor), resample=Image.NEAREST))
        
        rotated_cell_mask = np.transpose(resized_cell_mask)
        skel, dist = morphology.medial_axis(rotated_cell_mask, return_distance=True)

        return dist, rotated_cell_mask
    
    
    def fit_medial_axis(rotated_distance_mask, end_cap, resize_factor, order, smoothing):
        
        x_coords = []
        y_coords = []
        
        i = 0
        for x in np.nditer(rotated_distance_mask, flags=['external_loop'], order='F'):
            if rotated_distance_mask[np.argmax(x)][i] > end_cap * resize_factor:
                y_coords.append(np.argmax(x))
                x_coords.append(i)
            i+=1
        
        lin_fit = interpolate.UnivariateSpline(x_coords, y_coords, k=order, s=smoothing)
        
        return x_coords, list(lin_fit(x_coords)), i

    
    def extend_cell_caps(x_coords, y_coords, max_mask_coord):
        
        lcap_fit = np.polyfit(x_coords[0:10], y_coords[0:10], 1)
        rcap_fit = np.polyfit(x_coords[-10:], y_coords[-10:], 1)
        
        l_cap = np.polyval(lcap_fit, np.arange(0,min(x_coords)))
        r_cap = np.polyval(rcap_fit, np.arange(max(x_coords)+1,max_mask_coord))
        
        x_coords_cap = list(np.arange(0,min(x_coords))) + x_coords + list(np.arange(max(x_coords)+1,max_mask_coord))
        y_coords_cap = list(l_cap) + y_coords + list(r_cap)
        
        return np.array(x_coords_cap), np.array(y_coords_cap)
    
    
    def rotate_coords(point):
        return point[1], point[0]
    
    
    def trim_central_line(x_coords, y_coords, cropped_cell_mask):
        
        non_zero_y = np.nonzero(cropped_cell_mask)[0]
        cell_y_coords = np.intersect1d(non_zero_y, y_coords.astype(int))
        cell_coords_mask = np.isin(y_coords.astype(int), cell_y_coords)
        x_coords = x_coords[cell_coords_mask]
        y_coords = y_coords[cell_coords_mask]
        
        return x_coords[:-9], y_coords[:-9]
    
    
    def center_medial_axis(x_coords, y_coords, crop_pad):
        
        x_coords_adj = x_coords + crop_pad[2]
        y_coords_adj = y_coords + crop_pad[0]
        
        medial_axis_df = pd.DataFrame()
        medial_axis_df['x'] = x_coords_adj
        medial_axis_df['y'] = y_coords_adj
        # get the arch length of the medial axis
        delta_x_sqr = (x_coords_adj[1:] - x_coords_adj[0:-1])**2
        delta_y_sqr = (y_coords_adj[1:] - y_coords_adj[0:-1])**2
        disp_array = np.sqrt(delta_x_sqr + delta_y_sqr)
        disp_list = [0]
        for disp in disp_array:
            disp_list.append(disp_list[-1]+disp)
        medial_axis_df['arch_length_centered'] = disp_list - np.max(disp_list)/2
        medial_axis_df['arch_length_scaled'] = medial_axis_df['arch_length_centered'] / medial_axis_df['arch_length_centered'].max()
        
        return medial_axis_df
    
    
    def estimate_medial_axis(cell_mask, crop_frame, 
                             resize_factor, order, 
                             smoothing, end_cap, show):
        
        cropped_cell_mask, crop_pad = crop_cell_mask(cell_mask, crop_frame)
        rotated_distance_mask, rotated_cell_mask = resize_rotate_image(cropped_cell_mask, resize_factor)
        x_coords, y_coords, max_mask_coord = fit_medial_axis(rotated_distance_mask, end_cap, resize_factor, order, smoothing)
        x_coords, y_coords = extend_cell_caps(x_coords, y_coords, max_mask_coord)
#        x_coords, y_coords = rotate_coords((x_coords/resize_factor, y_coords/resize_factor), -rotation_angle)
#        x_coords = x_coords + cropped_cell_mask.shape[1]
        x_coords, y_coords = rotate_coords((x_coords/resize_factor, y_coords/resize_factor))
        x_coords, y_coords = trim_central_line(x_coords, y_coords, cropped_cell_mask)
        medial_axis_df = center_medial_axis(x_coords, y_coords, crop_pad)
        
        if show == True:
            plt.imshow(cropped_cell_mask)
            plt.plot(medial_axis_df.x-crop_pad[2], medial_axis_df.y-crop_pad[0])
            plt.show()
        
        return medial_axis_df
    
    
    return estimate_medial_axis(cell_mask, crop_frame, resize_factor, order, smoothing, end_cap, show)


def all_medial_axis(results_path, index,
                    crop_frame=3, 
                    resize_factor=10, order=2, 
                    smoothing=3000, end_cap=3, show=False):
    """
    This function was used to draw all medial axes for all cells.
    See example in the Jupyter notebook: Microfluidics segmentation and tracking
    """    
    
    cell_dataframe = pd.read_pickle(results_path+'/cell_dataframe', compression='zip')
    
    exp_pos_list = []
    for exp in cell_dataframe.experiment.unique():
        for pos in cell_dataframe[cell_dataframe.experiment==exp].position_int.unique():
            exp_pos_list.append((exp, pos))
    
    position_dataframe = cell_dataframe[(cell_dataframe.experiment==exp_pos_list[index][0])&
                                        (cell_dataframe.position_int==exp_pos_list[index][1])]

    print('getting the cell labels...')
    cell_labels_dict = get_specific_cell_labels(results_path, exp_pos_list[index][0], exp_pos_list[index][1], curated=True)
    
    prefix = exp_pos_list[index][0] + '_pos'+str(exp_pos_list[index][1])
    
    if os.path.exists(results_path +'/'+prefix+'_medial_axis_dict'):
        print('Medial axis definition found. Loading...')
        medial_axis_dict = pickle.load(open(results_path +'/'+prefix+'_medial_axis_dict', 'rb'))
        bad_cells = pickle.load(open(results_path +'/'+prefix+'_bad_cell_axis', 'rb'))
    else:
        print('Drawing medial axes for experiment:',exp_pos_list[index][0],' in position:',str(exp_pos_list[index][1]))
        medial_axis_dict = {}
        bad_cells = []
    
        for index, row in position_dataframe.iterrows():
            cell_id = row.cell_id
            label = row.label
            frame = row.frame
            experiment = row.experiment
            position = row.position_int
            cell_mask = cell_labels_dict[experiment][position][frame]==label
            print(cell_id)
            
            try:
                medial_axis_dict[cell_id] = get_medial_axis(cell_mask, crop_frame, 
                                                            resize_factor, order, 
                                                            smoothing, end_cap, show)
            except:
                medial_axis_dict[cell_id] = np.nan
                bad_cells.append(cell_id)
                
                
        pickle.dump(medial_axis_dict, open(results_path +'/'+prefix+'_medial_axis_dict', 'wb'))
        pickle.dump(bad_cells, open(results_path +'/'+prefix+'_bad_cell_axis', "wb"))
            
    return medial_axis_dict, bad_cells


def check_medial_axis(results_path):
    """
    This funciton cna be used for the visual inspection of the medial axes. 
    """
    cell_dataframe = pd.read_pickle(results_path+'/cell_dataframe', compression='zip')
    medial_axis_dict = load_medial_axes(results_path)
#    bad_cells = pickle.load(open(results_path+'/bad_cell_axis', 'rb'))
    cell_labels_dict = get_cell_labels(results_path, curated=True)
            
    for index, row in cell_dataframe.iterrows():
        cropped_mask, crop_pad = crop_cell_mask(cell_labels_dict[row.experiment][row.position_int][row.frame]==row.label, crop_frame=3)
        plt.figure(figsize=(2,10))
        plt.imshow(cropped_mask)
#        plt.imshow(ce1ll_labels_dict[row.experiment][row.position_int][row.frame]==row.label)
        if row.cell_id in medial_axis_dict:
            plt.plot(medial_axis_dict[row.cell_id].x-crop_pad[2], medial_axis_dict[row.cell_id].y-crop_pad[0])
#        plt.plot(medial_axis_dict[row.cell_id].x, medial_axis_dict[row.cell_id].y)
        plt.show()
        input()


def load_medial_axes(results_path):
    """
    Returns a dictionary with all medial axes.
    The keys are the cell IDs and the values correspond to Pandas dataframes
    with the absolute and the relative coordinates of each medial axis. 
    """
    medial_axis_dict = {}
    dir_list = os.listdir(results_path)
    label_dir = [s for s in dir_list if 'medial_axis_dict' in s]
    for s in label_dir:
        print(s)
        pre_dict = pickle.load(open(results_path+'/'+s, 'rb'))
        medial_axis_dict = {**medial_axis_dict, **pre_dict}

    return medial_axis_dict


def load_specific_medial_axes(results_path, experiment, position):
    """
    Returns medial axes for a specific xy position and microfluidic channel position.
    """
    dir_list = os.listdir(results_path)
    label_dir = [s for s in dir_list if 'medial_axis_dict' in s]
    exp_dir = [s for s in label_dir if experiment in s]
    pos_dir = [s for s in exp_dir if 'pos'+str(position)+'_' in s]

    if len(pos_dir)==0:
        print('There is no medial axis definition for this experiment or position')
        medial_axis_dict=np.nan
    elif len(pos_dir)==1:
        print('Loading medial axis:',pos_dir[0])
        medial_axis_dict =pickle.load(open(results_path+'/'+pos_dir[0], 'rb'))
    elif len(pos_dir)>1:
        print('Multiple medial axis definitions found', pos_dir)
        medial_axis_dict=np.nan
    
    return medial_axis_dict


def get_cell_length(cell_dataframe, medial_axis_dict):
    """
    Returns a dictionary of the cell lengths (values - the arch length of the medial axis), 
    to each corresponding cell ID (keys).
    """
#    cell_dataframe = pd.read_pickle(results_path+'/cell_dataframe', compression='zip')
#    medial_axis_dict = pickle.load(open(results_path+'/medial_axis_dict', 'rb'))
    cell_length_dict = {}
    for cell in medial_axis_dict:
#        print(cell)
        if type(medial_axis_dict[cell])!=float:
            cell_length_dict[cell] = medial_axis_dict[cell].arch_length_centered.max()*2
    cell_dataframe['cell_length'] = cell_dataframe.cell_id.map(cell_length_dict)
#    cell_dataframe.to_pickle(results_path+'/cell_dataframe', compression='zip')
    return cell_dataframe




# ------------------   6. CELL POLARITY ESTIMATION BASED ON THE MEDIAL AXIS AND THE CELL DIVISION CYCLE  --------------------- #
def get_polarity(mother_cell_centroid, daughter_medial_axis_df):
    """
    Returns the polarity of each cell division cycle trajectory (1 or -1).
    """
    daughter_medial_axis_df['mother_distance'] = np.sqrt((daughter_medial_axis_df.x-mother_cell_centroid[0])**2 + (daughter_medial_axis_df.y-mother_cell_centroid[1])**2)
    
    if daughter_medial_axis_df.mother_distance.values[0] > daughter_medial_axis_df.mother_distance.values[-1]:
        if daughter_medial_axis_df.arch_length_scaled.values[0] > 0:
            polarity = -1
        elif daughter_medial_axis_df.arch_length_scaled.values[0] < 0:
            polarity = 1
    elif daughter_medial_axis_df.mother_distance.values[0] < daughter_medial_axis_df.mother_distance.values[-1]:
        if daughter_medial_axis_df.arch_length_scaled.values[0] > 0:
            polarity = 1
        elif daughter_medial_axis_df.arch_length_scaled.values[0] < 0:
            polarity = -1
    
    return polarity


def apply_polarity(cell_dataframe, medial_axis_dict):
    """
    Applies the get_polarity function.
    See example in the Jupyter notebook: Microfluidics segmentation and tracking
    """
    if 'mother_cell' not in cell_dataframe:
        print('get the mother cells first...')
        return()
    if 'zeroed_frame' not in cell_dataframe:
        print('get the cell cycle phases first...')
        return()
        
    polarity_dict = {}
    
    zero_df = cell_dataframe[~cell_dataframe.mother_cell.isnull()]
    zero_df = zero_df[zero_df.cc_phase==0]
    for index, row in zero_df.iterrows():
        daughter_medial_axis_df = medial_axis_dict[row.cell_id]
        if type(daughter_medial_axis_df)==pd.core.frame.DataFrame:
            #print('getting polarity for:', row.cell_trajectory_id)
            mother_cell_centroid = cell_dataframe.loc[cell_dataframe.cell_id==row.mother_cell, ['x', 'y']].values[0]
            polarity_dict[row.cell_trajectory_id]=get_polarity(mother_cell_centroid, daughter_medial_axis_df)
        elif type(daughter_medial_axis_df)==float:
            print('This daughter cell does not have a medial axis:', row.cell_trajectory_id)
    
    cell_dataframe['polarity'] = cell_dataframe.cell_trajectory_id.map(polarity_dict)
    cell_dataframe = cell_dataframe.sort_values(['cell_trajectory_id', 'frame'])
    
    return cell_dataframe, polarity_dict


def check_polarity(cell_dataframe, results_path, 
                   save_path):
    cell_labels_dict = get_cell_labels(results_path, curated=True)
    """
    This function can be used to visualy inspect the calculated polarities. 
    """
    for exp in cell_labels_dict:
        for pos in cell_labels_dict[exp]:
            for fr in cell_labels_dict[exp][pos]:
                plt.figure(figsize=(3,10))
                plt.text(x=5, y=10, s=str(fr)+' frame', color='white')
                plt.imshow(cell_labels_dict[exp][pos][fr])
                frame_df = cell_dataframe[(cell_dataframe.experiment==exp)&
                                          (cell_dataframe.position_int==pos)&
                                          (cell_dataframe.frame==fr)]
                for index, row in frame_df.iterrows():
                    plt.plot(row.x, row.y, 'o', color='black')
                    plt.text(x=row.x+5, y=row.y, s=str(row.polarity), color='white')
                if os.path.isdir(save_path):
                    plt.savefig(save_path+'/'+row.experiment+'_pos'+str(row.position_int)+'_fr'+str(row.frame)+'.jpeg')
                plt.show()
                

# ------------------   7. RELATIVE PIXEL COORDINATES  --------------------- #
def get_oned_coordinates(cell_mask, medial_axis_df, half_window): 
    """
    This function returns the projection of each fluorescent pixel on the medial axis 
    and its position across the cell width.
    """    
    cell_mask_df = pd.DataFrame()
    cell_mask_df['x'] = np.nonzero(cell_mask)[1]
    cell_mask_df['y'] = np.nonzero(cell_mask)[0]
#    cell_mask_df['z'] = fluor_image[np.nonzero(cell_mask)]

    def get_pixel_projection(pixel_x, pixel_y, medial_axis_df, half_window):
        
        medial_axis_df['pixel_distance'] = np.sqrt((medial_axis_df.x-pixel_x)**2+(medial_axis_df.y-pixel_y)**2)
        min_df = medial_axis_df[medial_axis_df.pixel_distance == medial_axis_df.pixel_distance.min()]
        min_arch_centered_length = min_df.arch_length_centered.values[0]
        min_arch_scaled_length =  min_df.arch_length_scaled.values[0]
        min_distance_abs = min_df.pixel_distance.values[0]
        min_index = min_df.index.values[0]
        medial_axis_coords = (min_df.x.values[0], min_df.y.values[0])
        
        def get_relative_distance(min_distance_abs, medial_axis_df, min_index, medial_axis_coords, pixel_x, pixel_y, half_window):
    
            if min_index>=half_window and min_index<medial_axis_df.index.max()-half_window:
                index_range = (min_index-half_window, min_index+half_window)
            elif min_index<half_window and min_index<medial_axis_df.index.max()-half_window:
                index_range = (0, min_index+half_window)
            elif min_index>=half_window and min_index>=medial_axis_df.index.max()-half_window:
                index_range = (min_index-half_window, medial_axis_df.index.max())
            
            delta_x = (medial_axis_df.iloc[index_range[1]].x -  medial_axis_df.iloc[index_range[0]].x)
            delta_y = (medial_axis_df.iloc[index_range[1]].y -  medial_axis_df.iloc[index_range[0]].y)
            medial_axis_vector = [delta_x, delta_y]
            
            delta_x = pixel_x - medial_axis_coords[0]
            delta_y = pixel_y - medial_axis_coords[1]
            pixel_vector = [delta_x, delta_y]
            
            cross_product = np.cross(medial_axis_vector, pixel_vector)
            if cross_product != 0:
                min_distance = np.sign(cross_product)*min_distance_abs
#                return min_distance
            elif cross_product == 0:
#                half_window+=1
#                get_relative_distance(min_distance_abs, medial_axis_df, min_index, medial_axis_coords, pixel_x, pixel_y, half_window)
                min_distance = 0
            return min_distance
        
        min_distance = get_relative_distance(min_distance_abs, medial_axis_df, min_index, medial_axis_coords, pixel_x, pixel_y, half_window)
    
        return min_arch_centered_length, min_arch_scaled_length, min_distance
        
            
    cell_mask_df['oned_coords'] = cell_mask_df.apply(lambda x: get_pixel_projection(x.x, x.y, medial_axis_df, half_window=5), axis=1)
    cell_mask_df[['arch_length', 'scaled_length', 'width']] = pd.DataFrame(cell_mask_df.oned_coords.to_list(), index=cell_mask_df.index)
    cell_mask_df = cell_mask_df.drop(['oned_coords'], axis=1)
    
    return cell_mask_df
    

def apply_oned_coordinates(position_dataframe, cell_labels_dict, medial_axis_dict, images_path,
                           results_path, fluor_channels, edge_width, half_window, every_nth):
    """
    Applies the get_oned_coordinates functios.
    """
    exp = position_dataframe.experiment.values[0]
    pos = position_dataframe.position_int.values[0]
    xy_pos = exp[exp.find('xy'):]
    
    
    polar_df = position_dataframe[~position_dataframe.polarity.isnull()]
    polar_df = polar_df[polar_df.frame%every_nth==0]
    oned_coords_dict = {}
#    pickle.dump(oned_coords_dict, open(results_path+'/'+exp +'_pos'+str(pos)+'_oned_coords_dict', 'wb'))
    
    for index, row in polar_df.iterrows():
        try:
            print('mapping the pixel 1D coordinates for cell:', row.cell_id)
            cell_mask = cell_labels_dict[row.experiment][row.position_int][row.frame]==row.label
            medial_axis_df = medial_axis_dict[row.cell_id]
            medial_axis_df['arch_length_centered'] = medial_axis_df.arch_length_centered*row.polarity
            medial_axis_df['arch_length_scaled'] = medial_axis_df.arch_length_scaled*row.polarity
            cell_mask_df = get_oned_coordinates(cell_mask, medial_axis_df, half_window)
            cell_mask_df['width'] = cell_mask_df.width*row.polarity
            for ch in fluor_channels:
                try:
                    fluor_image = get_fluorescence_image(images_path, xy_pos, ch, pos, row.frame, edge_width)
                except IndexError:
                    xy_pos_2 = 'xy'+str(int(xy_pos[2:]))
                    fluor_image = get_fluorescence_image(images_path, xy_pos_2, ch, pos, row.frame, edge_width)
                cell_mask_df['fluor_pixels_'+str(ch)] = fluor_image[np.nonzero(cell_mask)]
            oned_coords_dict[row.cell_id] = cell_mask_df
            
        except AttributeError:
            print('This cell does not have a medial axis')
            
    pickle.dump(oned_coords_dict, open(results_path+'/'+exp +'_pos'+str(pos)+'_oned_coords_dict', 'wb'))
    
    return oned_coords_dict
        

def parallel_oned_estimation(results_path, images_path, 
                             fluorescent_channels, edge_width, half_window, every_nth):
    """
    Parallel application of the apply_oned_coordinates function
    """
    print('getting the cell dataframe...')
    cell_dataframe = pd.read_pickle(results_path+'/cell_dataframe', compression='zip')

    print('getting the medial axes...')
    medial_axis_dict = load_medial_axes(results_path)
#    bad_cells = pickle.load(open(results_path+'/bad_cell_axis', 'rb'))
    print('getting the cell labels...')
    cell_labels_dict = get_cell_labels(results_path, curated=True)
    
    exp_pos_list = []
    for exp in cell_dataframe.experiment.unique():
        for pos in cell_dataframe[cell_dataframe.experiment==exp].position_int.unique():
            exp_pos_list.append((exp, pos))
    
    def parallel_oned_estimation(index):
        position_dataframe = cell_dataframe[(cell_dataframe.experiment==exp_pos_list[index][0])&
                                            (cell_dataframe.position_int==exp_pos_list[index][1])]
        apply_oned_coordinates(position_dataframe, cell_labels_dict, medial_axis_dict, 
                               images_path, results_path, fluorescent_channels, edge_width, half_window, every_nth)

    if __name__ == '__main__':   
        pool = mp.Pool(mp.cpu_count())
        # Create a multiprocessing Pool
        pool.map(parallel_oned_estimation, range(0,len(exp_pos_list)))  # process data_inputs iterable with pool
    pool.close() 


#    @njit(parallel=True)
#    def parallel_estimation():
#        
#        for index in prange(len(exp_pos_list)):
#            print(index)
#            print(exp_pos_list[index])
#            position_dataframe = cell_dataframe[(cell_dataframe.experiment==exp_pos_list[index][0])&
#                                                (cell_dataframe.position_int==exp_pos_list[index][1])]
#            apply_oned_coordinates(position_dataframe, cell_labels_dict, medial_axis_dict, 
#                                   fluorescence_images_dict, results_path, half_window=5)
#
#    parallel_estimation()



def run_oned_estimation_single(results_path, images_path, half_window, edge_width, fluorescent_channels, index, every_nth):
    """
    Non-parallel application of the apply_oned_coordinates function
    """
    print('getting the cell dataframe...')
    cell_dataframe = pd.read_pickle(results_path+'/cell_dataframe', compression='zip')
    
    exp_pos_list = []
    for exp in cell_dataframe.experiment.unique():
        for pos in cell_dataframe[cell_dataframe.experiment==exp].position_int.unique():
            exp_pos_list.append((exp, pos))
            
    
    print('getting the medial axes...')
    medial_axis_dict = load_specific_medial_axes(results_path, exp_pos_list[index][0], exp_pos_list[index][1])
    #    bad_cells = pickle.load(open(results_path+'/bad_cell_axis', 'rb'))
    print('getting the cell labels...')
    cell_labels_dict = get_specific_cell_labels(results_path, exp_pos_list[index][0], exp_pos_list[index][1], curated=True)
    
    position_dataframe = cell_dataframe[(cell_dataframe.experiment==exp_pos_list[index][0])&
                                        (cell_dataframe.position_int==exp_pos_list[index][1])]
    oned_coords_dict = apply_oned_coordinates(position_dataframe, cell_labels_dict, medial_axis_dict, 
                                              images_path, results_path, fluorescent_channels, edge_width, half_window, every_nth)
    return oned_coords_dict


def get_oned_coords(results_path):
    """
    Returns a dictionary that includes all realtive pixels coordinates (Pandas dataframe) as values,
    for each corresponding cell ID (keys).
    """
    oned_coords_dict = {}
    dir_list = os.listdir(results_path)
    label_dir = [s for s in dir_list if 'oned_coords_dict' in s]
    for s in label_dir:
        print(s)
        pre_dict = pickle.load(open(results_path+'/'+s, 'rb'))
        oned_coords_dict = {**oned_coords_dict, **pre_dict}

    return oned_coords_dict       


def get_twod_dataframe(cell_dataframe, oned_coords_dict, sample_size):
    """
    Converts the dictionary returned by the get_oned_coords function into a Pandas dataframe,
    where the cell ID is a column.
    
    This dataframe was used for the reconstrunction of the average 2D fluorescence.
    """
    cell_dataframe = cell_dataframe[cell_dataframe.cell_id.isin(list(oned_coords_dict.keys()))]
    if sample_size > 0:    
        cell_dataframe = cell_dataframe.sample(sample_size, random_state=1)
    i=0
    for index, row in cell_dataframe.iterrows():
        if row.cell_id in oned_coords_dict:
            if i==0:
                oned_coords_df = oned_coords_dict[row.cell_id]
                oned_coords_df['cell_id'] = row.cell_id
                i+=1
            elif i>0:
                temp_df = oned_coords_dict[row.cell_id]
                temp_df['cell_id'] = row.cell_id
                oned_coords_df = pd.concat([oned_coords_df, temp_df])
            # print(row.cell_id)
    
    return oned_coords_df



def get_demo_graph_v2(results_path, mother_df, fluor_channels, bin_n, rol_win, sort_col):
    """
    This function was used to reconstruct the ensemble kymographs from the 1D fluorescent pixel projections.
    
    It returns two 2D arrays of the average or scaled cell fluorescence along the cell length (primary axis),
    for each single cell sorted by the sorting dimension (secondary axis).
    
    The sorting dimension can be any column included i nthe mother_df pandas dataframe, 
    such as the cell area, the cell division cycle phase or tge cell length. 
    For the ensemble kymographs the cell IDs were sorted first by the cell division cycle and then by cell length.
    """
    print('assembling the oneD pixel coordinates...')
    oned_coords_dict = get_oned_coords(results_path)
    
    fluor_arrays_dict = {}
    for ch in fluor_channels:
        fluor_arrays_dict[ch] = []
        
    if 'cell_length' not in mother_df:
        print('getting the medial axes...')
        medial_axis_dict = load_medial_axes(results_path)
        print('getting the cell length...')
        mother_df = get_cell_length(mother_df, medial_axis_dict)
    

    mother_dataframe = mother_df
    print('removing cells without estimated polarity...')
    mother_dataframe = mother_dataframe[~mother_dataframe.polarity.isnull()]
    print('sorting the cell cycle trajectories...')
    mother_dataframe = mother_dataframe.sort_values(sort_col)
    
    i=0
    for index, row in mother_dataframe.iterrows():
#        print(row.cell_id)
        if row.cell_id in oned_coords_dict:
            oned_coords = oned_coords_dict[row.cell_id]
            oned_coords['scaled_length_bin'] = pd.cut(oned_coords.scaled_length, bins=np.arange(-1,1+2/bin_n,2/bin_n), right=True, include_lowest=True, labels=False)
            mean_df = oned_coords.groupby('scaled_length_bin').mean().reindex(list(range(0,bin_n)))
            mean_df = mean_df.interpolate()
            smooth_df = mean_df.rolling(rol_win, min_periods=1, center=True).mean()
            
            for ch in fluor_channels:
                fluor_arrays_dict[ch].append(smooth_df['fluor_pixels_'+str(ch)].tolist())
        i += 1
        if i in [int(mother_dataframe.shape[0]/4), int(mother_dataframe.shape[0]/3), int(mother_dataframe.shape[0]/2), int(mother_dataframe.shape[0]/1.5), int(mother_dataframe.shape[0]/1.2)]:
            print('analyzed', i, 'out of', mother_dataframe.shape[0], 'cells in complete cell cycles...')
    
    print('estimating the z-score...')
    fluor_arrays_dict_z = {}
    for ch in fluor_arrays_dict:
        fluor_arrays_dict_z[ch] = []
        for cl in fluor_arrays_dict[ch]:
            cl = np.array(cl)
            fluor_arrays_dict_z[ch].append(list((cl-np.mean(cl))/np.std(cl)))
    
    return fluor_arrays_dict, fluor_arrays_dict_z


def map_particle_positions(cell_df, particle_df, half_window, results_path):
    """
    This function is used to map the relative cell coordinates of the Gaussian center of detected particles.
    It adds the particle position statistics to the cell_df (or the mother_df).
    """
    def get_position_projection(pixel_x, pixel_y, medial_axis_df, half_window):
        
        medial_axis_df['pixel_distance'] = np.sqrt((medial_axis_df.x-pixel_x)**2+(medial_axis_df.y-pixel_y)**2)
        min_df = medial_axis_df[medial_axis_df.pixel_distance == medial_axis_df.pixel_distance.min()]
        min_arch_centered_length = min_df.arch_length_centered.values[0]
        min_arch_scaled_length =  min_df.arch_length_scaled.values[0]
        min_distance_abs = min_df.pixel_distance.values[0]
        min_index = min_df.index.values[0]
        medial_axis_coords = (min_df.x.values[0], min_df.y.values[0])
        
        def get_relative_distance(min_distance_abs, medial_axis_df, min_index, medial_axis_coords, pixel_x, pixel_y, half_window):
    
            if min_index>=half_window and min_index<medial_axis_df.index.max()-half_window:
                index_range = (min_index-half_window, min_index+half_window)
            elif min_index<half_window and min_index<medial_axis_df.index.max()-half_window:
                index_range = (0, min_index+half_window)
            elif min_index>=half_window and min_index>=medial_axis_df.index.max()-half_window:
                index_range = (min_index-half_window, medial_axis_df.index.max())
            
            delta_x = (medial_axis_df.iloc[index_range[1]].x -  medial_axis_df.iloc[index_range[0]].x)
            delta_y = (medial_axis_df.iloc[index_range[1]].y -  medial_axis_df.iloc[index_range[0]].y)
            medial_axis_vector = [delta_x, delta_y]
            
            delta_x = pixel_x - medial_axis_coords[0]
            delta_y = pixel_y - medial_axis_coords[1]
            pixel_vector = [delta_x, delta_y]
            
            cross_product = np.cross(medial_axis_vector, pixel_vector)
            if cross_product != 0:
                min_distance = np.sign(cross_product)*min_distance_abs
#                return min_distance
            elif cross_product == 0:
#                half_window+=1
#                get_relative_distance(min_distance_abs, medial_axis_df, min_index, medial_axis_coords, pixel_x, pixel_y, half_window)
                min_distance = 0
            return min_distance
        
        min_distance = get_relative_distance(min_distance_abs, medial_axis_df, min_index, medial_axis_coords, pixel_x, pixel_y, half_window)
    
        return min_arch_centered_length, min_arch_scaled_length, min_distance
    

    print('getting the medial axes...')
    medial_axis_dict = load_medial_axes(results_path)
#    print('getting the cell labels...')
#    cell_labels_dict = get_cell_labels(results_path, curated=True)
    
    vol_dict = {}
    arch_dict = {}
    rel_dict = {}
    dis_dict = {}
    pos_dict = {}
    for index, row in particle_df.iterrows():
        
        par_df = cell_df[(cell_df.experiment==row.experiment)&(cell_df.position_int==row.position)&(cell_df.frame==row.frame)&(cell_df.label==row.label)]
        if par_df.shape[0] == 1:
            cl = par_df.cell_id.values[0]
            if cl in medial_axis_dict and type(medial_axis_dict[cl])!=float:
                if cl not in pos_dict:
                    print(cl)
                    par_stats = get_position_projection(row.gaussian_center[0], row.gaussian_center[1], medial_axis_dict[cl], half_window)
                    pos_dict[cl] = row.gaussian_center
                    vol_dict[cl] = row.gaussian_volume
                    arch_dict[cl] = par_stats[0] * par_df.polarity.values[0]
                    rel_dict[cl] = par_stats[1] * par_df.polarity.values[0]
                    dis_dict[cl] = par_stats[2] * par_df.polarity.values[0]
                else:
                    if row.gaussian_volume > vol_dict[cl]:
                        print('found bigger particle. Recalculating:', cl)
                        par_stats = get_position_projection(row.gaussian_center[0], row.gaussian_center[1], medial_axis_dict[cl], half_window)
                        pos_dict[cl] = row.gaussian_center
                        vol_dict[cl] = row.gaussian_volume
                        arch_dict[cl] = par_stats[0] * par_df.polarity.values[0]
                        rel_dict[cl] = par_stats[1] * par_df.polarity.values[0]
                        dis_dict[cl] = par_stats[2] * par_df.polarity.values[0]
                        
    cell_df['particle_pos'] = cell_df.cell_id.map(pos_dict)
    cell_df['particle_vol'] = cell_df.cell_id.map(vol_dict)
    cell_df['particle_arch_length'] = cell_df.cell_id.map(arch_dict)
    cell_df['particle_relative_length'] = cell_df.cell_id.map(rel_dict)
    cell_df['particle_distance'] = cell_df.cell_id.map(dis_dict)
    
    return cell_df


def get_object_relative_coords(medial_axis_df, polarity, object_coord_list):
    """
    This function returns the cell length (absolute and relative) coordinates 
    of the centroid of a masked object.
    """
    if type(medial_axis_df) != float:
        medial_axis_df['arch_length_centered'] = medial_axis_df.arch_length_centered*polarity
        medial_axis_df['arch_length_scaled'] = medial_axis_df.arch_length_scaled*polarity
        
        relative_coords_list = []
        arch_coords_list = []
        
        if type(object_coord_list)!=float:
            for obj in object_coord_list:
                medial_axis_df['distance'] = np.sqrt((medial_axis_df.x - obj[0])**2 + (medial_axis_df.y - obj[1])**2)
                min_list = medial_axis_df[medial_axis_df.distance == medial_axis_df['distance'].min()]
                relative_coords_list.append(min_list.arch_length_scaled.values[0])
                arch_coords_list.append(min_list.arch_length_centered.values[0])
            return arch_coords_list, relative_coords_list, len(object_coord_list)
        else:
            return (np.nan, np.nan, np.nan)
    else:
        return (np.nan, np.nan, np.nan)


# ------------------   8. CELL AND DETECTED OBJECT (e.g. NUCLEOID) FLUROESCENCE STATISTICS  --------------------- #
def get_mean_fluorescence(fluorescence_images_dict, cell_labels_dict, 
                          experiment, position, frame, label, every_nth):
    """
    This function returns the mean cell fluorescence and its standard deviation.
    The mean fluorescence corresponds to the average concentration of the reported protein.
    """
    if frame%every_nth != 0:
        return (np.nan, np.nan)
    else:
        xy_pos = experiment[experiment.find('xy'):]
        try:
            fluorescence_image = fluorescence_images_dict[xy_pos][position][frame]
        except KeyError:
            xy_pos_2 = 'xy'+str(int(xy_pos[2:]))
            fluorescence_image = fluorescence_images_dict[xy_pos_2][position][frame]
        cell_masks = cell_labels_dict[experiment][position][frame]
        pixel_fluor = fluorescence_image[cell_masks==label]
        mean_fluor = pixel_fluor.mean()
        std_fluor = pixel_fluor.std()
        
        return (mean_fluor, std_fluor)



def apply_fluorescence_stats(cell_dataframe, images_path, results_path, channel, edge_width, every_nth):
    """
    This function applies the get_mean_fluorescence function to all segmented and tracked cells and
    for one specific fluroescence channel.
    """
    print('getting fluorescence images...')
    frame_range = (cell_dataframe.frame.min(), cell_dataframe.frame.max()+1)
    fluorescence_images_dict = get_fluorescence_images(images_path, channel, frame_range, edge_width)
    print('getting cell labels...')
    cell_labels_dict = get_cell_labels(results_path, curated=True)
    
    if 'cc_phase' not in cell_dataframe:
        print('estimating the cell cycle phases...')
        cell_dataframe = get_cell_cycle_phases(cell_dataframe)
    
#    print('estimating cell polarity...')
#    cell_dataframe = get_polarity(cell_dataframe)
    
    print('applying fluorescence statistics...')
    cell_dataframe['fluor_analysis'] = cell_dataframe.apply(lambda x: get_mean_fluorescence(fluorescence_images_dict, cell_labels_dict, x.experiment, x.position_int, x.frame, x.label, every_nth), axis=1)
    cell_dataframe[['mean_fluor_'+str(channel), 'std_fluor_'+str(channel), ]] = pd.DataFrame(cell_dataframe.fluor_analysis.tolist(), index=cell_dataframe.index)
    cell_dataframe['total_fluor_'+str(channel)] = cell_dataframe.area_px * cell_dataframe['mean_fluor_'+str(channel)]

    cc_mean_fluor_dict = cell_dataframe.groupby('cell_trajectory_id').mean()['mean_fluor_'+str(channel)].to_dict()
    cc_std_fluor_dict = cell_dataframe.groupby('cell_trajectory_id').std()['mean_fluor_'+str(channel)].to_dict()

    cell_dataframe['cc_mean_fluor_'+str(channel)] = cell_dataframe.cell_trajectory_id.map(cc_mean_fluor_dict)
    cell_dataframe['cc_std_fluor_'+str(channel)] = cell_dataframe.cell_trajectory_id.map(cc_std_fluor_dict)
    cell_dataframe['z_score_fluor_'+str(channel)] = (cell_dataframe['mean_fluor_'+str(channel)]-cell_dataframe['cc_mean_fluor_'+str(channel)])/cell_dataframe['cc_std_fluor_'+str(channel)]
    
    cell_dataframe = cell_dataframe.drop(['fluor_analysis'], axis=1)
    cell_dataframe = cell_dataframe.sort_values(['cell_trajectory_id', 'frame'])
    
    return cell_dataframe


def apply_all_fluorescent_stats(images_path, results_path, channels, every_nth, edge_width):
    """
    This function applies the apply_fluorescence_stats function to all specified channels.
    See example in the Jupyter notebook: Microfluidics segmentation and tracking
    """
    print('getting the cell dataframe...')
    cell_dataframe = pd.read_pickle(results_path+'/cell_dataframe', compression='zip')
    
    for ch in channels:
        cell_dataframe = apply_fluorescence_stats(cell_dataframe, images_path, results_path, ch, edge_width, every_nth)
    
    cell_dataframe.to_pickle(results_path+'/cell_dataframe_fluor', compression='zip')
    
    return cell_dataframe


def get_masked_fluorescence(fluorescence_images_dict, cell_labels_dict, object_masks_dict,
                            experiment, position, frame, label, every_nth, min_mask_size):
    """
    This function is used to apply the fluorescence statistics to a masked cell object (e.h. the nucleoid).
    """
    if frame%every_nth != 0:
        return (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)
    else:
        xy_pos = experiment[experiment.find('xy'):]
        try:
            fluor_image = fluorescence_images_dict[xy_pos][position][frame]
        except KeyError:
            xy_pos_2 = 'xy'+str(int(xy_pos[2:]))
            fluor_image = fluorescence_images_dict[xy_pos_2][position][frame]
        cell_masks = cell_labels_dict[experiment][position][frame]
        try:
            object_mask = object_masks_dict[xy_pos][position][frame]
#            print('Using xy position as object mask key')
        except KeyError:
#            print('Using experiment as object mask key')
            object_mask = object_masks_dict[experiment][position][frame]
        
        cell_mask = cell_masks==label
        object_mask = (object_mask*cell_masks)==label
        non_object_mask = (cell_mask.astype(int)+object_mask.astype(int))==1
        
        all_pixel_fluor = fluor_image[cell_mask]
        object_pixel_fluor = fluor_image[object_mask]
        non_object_pixel_fluor = fluor_image[non_object_mask]
        
        all_mean_fluor = all_pixel_fluor.mean()
        all_std_fluor = all_pixel_fluor.std()
        obj_mean_fluor = object_pixel_fluor.mean()
        obj_std_fluor = object_pixel_fluor.std()
        non_obj_mean_fluor = non_object_pixel_fluor.mean()
        non_obj_std_fluor = non_object_pixel_fluor.std()
        #cell_area = all_pixel_fluor.shape[0]
        obj_area = object_pixel_fluor.shape[0]
        non_obj_area = non_object_pixel_fluor.shape[0]
        
        obj_xy = []
        
        for obj_lbl in measure.regionprops(measure.label(object_mask)):
            if obj_lbl.area >= min_mask_size:
                obj_xy.append(obj_lbl.centroid[::-1])   
        
        return all_mean_fluor, all_std_fluor, obj_area, obj_mean_fluor, obj_std_fluor, non_obj_area, non_obj_mean_fluor, non_obj_std_fluor, obj_xy
                

def apply_object_fluorescence(results_path, images_path, channels, every_nth, edge_width,
                              min_mask_size=5, object_mask_dict_suffix='nucleoid_masks_dict'):
    """
    Applies the get_masked_fluorescence and the get_object_relative_coords functions.
    See example in the Jupyter notebook: Microfluidics segmentation and tracking
    """
    fluorescence_images_dict = {}
    print('getting the cell_dataframe')
    if os.path.exists(results_path+'/cell_dataframe_fluor_object_asymmetry'):
        open_path = results_path+'/cell_dataframe_fluor_object_asymmetry'
    elif os.path.exists(results_path+'/cell_dataframe_fluor_object'):
        open_path = results_path+'/cell_dataframe_fluor_object'
    elif  os.path.exists(results_path+'/cell_dataframe_fluor'):
        open_path = results_path+'/cell_dataframe_fluor'
    print('loading dataframe:', open_path)
    cell_dataframe = pd.read_pickle(open_path, compression='zip')
    print('getting fluorescence images...')
    frame_range = (cell_dataframe.frame.min(), cell_dataframe.frame.max()+1)
    for channel in channels:
        fluorescence_images_dict[channel] = get_fluorescence_images(images_path, channel, frame_range, edge_width)
    print('getting cell labels...')
    cell_labels_dict = get_cell_labels(results_path, curated=True)
    print('getting the medial axes...')
    medial_axis_dict = load_medial_axes(results_path)

    print('getting object masks...')
    if os.path.exists(results_path+'/'+object_mask_dict_suffix):
        object_masks_dict = pickle.load(open(results_path+'/'+object_mask_dict_suffix, 'rb'))
    else:
        print('No object segmentation file was detected. Perform object detection first...')
        return ()
    
    print('applying object statistics...')
    for channel in channels:
        print('...for channel',str(channel))
        cell_dataframe['fluor_analysis'] = cell_dataframe.apply(lambda x: get_masked_fluorescence(fluorescence_images_dict[channel], cell_labels_dict, object_masks_dict, x.experiment, x.position_int, x.frame, x.label, every_nth, min_mask_size), axis=1)
        cell_dataframe[['cell_mean_fluor_'+str(channel), 'cell_std_fluor_'+str(channel), 
                    'object_area', 'object_mean_fluor_'+str(channel), 'object_std_fluor_'+str(channel), 
                   'non_object_area', 'non_object_mean_fluor_'+str(channel), 
                    'non_object_std_fluor_'+str(channel), 'object_positions']] = pd.DataFrame(cell_dataframe.fluor_analysis.tolist(), index=cell_dataframe.index)
    print('mapping the nucleoid position on the central line...')
    cell_dataframe['relative_object_coords'] =  cell_dataframe.apply(lambda x: get_object_relative_coords(medial_axis_dict[x.cell_id], x.polarity, x.object_positions), axis=1)
    cell_dataframe[['object_centered_length', 'object_scaled_length', 'number_of_objects']] = pd.DataFrame(cell_dataframe.relative_object_coords.tolist(), index=cell_dataframe.index)
    print('cleaning the dataframe...')
    cell_dataframe = cell_dataframe.drop(['fluor_analysis', 'relative_object_coords'], axis=1)
    cell_dataframe = cell_dataframe.sort_values(['cell_trajectory_id', 'frame'])
    print('saving data...')
    cell_dataframe.to_pickle(open_path, compression='zip')
    
    return cell_dataframe


def get_nucleoid_asymmetry(array_of_centers):
    """
    Returns the mean cell length position of segmented objects in a cell.
    The cell length coordinates towards the old pole are negative,
    towards the new pole are negative,
    and the cell center is set to zero.
    
    This function can be applied to any cell length statistic (relative or absolute).
    """
    cnt = np.array(array_of_centers)
    cnt = np.array([np.nanmean(cnt[cnt<0]), np.nanmean(cnt[cnt>=0])])
    
    if cnt[~np.isnan(cnt)].shape[0] == 0:
        return np.nan
    elif cnt[~np.isnan(cnt)].shape[0] > 0:
        return np.nanmean(cnt)
    

def get_nucleoid_distance(array_of_centers):
    """
    Return the distance between two segmented objects.
    
    This function can be applied to any cell length statistic (relative or absolute).
    """
    cnt = np.array(array_of_centers)
    cnt = np.array([np.nanmean(cnt[cnt<0]), np.nanmean(cnt[cnt>=0])])
        
    if cnt[~np.isnan(cnt)].shape[0] == 0:
        return np.nan
    elif cnt[~np.isnan(cnt)].shape[0] == 1:
        return 0
    elif cnt[~np.isnan(cnt)].shape[0] > 1:
        return np.nansum(np.abs(cnt))


def get_fluorescence_asymmetry(oned_coords_dict, cell_id, cell_range, channel):
    """
    Returns the mean fluorescence of two mirrored cell length ranges around the cell center.
    Two fluorescence values are returned, one for each mirrored range.
    The cell length coordinates towards the old pole are negative,
    towards the new pole are negative,
    and the cell center is set to zero.
    """
    if cell_id in oned_coords_dict:
        oned_fluor_df = oned_coords_dict[cell_id]
        new_pole = oned_fluor_df[oned_fluor_df.scaled_length.between(cell_range[0],cell_range[1])]['fluor_pixels_'+str(channel)].mean()
        new_area = oned_fluor_df[oned_fluor_df.scaled_length.between(cell_range[0],cell_range[1])]['fluor_pixels_'+str(channel)].count()
        old_pole = oned_fluor_df[oned_fluor_df.scaled_length.between(-cell_range[1],-cell_range[0])]['fluor_pixels_'+str(channel)].mean()
        old_area = oned_fluor_df[oned_fluor_df.scaled_length.between(-cell_range[1],-cell_range[0])]['fluor_pixels_'+str(channel)].count()
    else:
        new_pole = np.nan
        new_area = np.nan
        old_pole = np.nan
        old_area = np.nan
    
    return new_pole, old_pole, new_area, old_area


def get_symmetric_mean(oned_coords_dict, cell_id, cell_range, channel):
    """
    Returns the mean fluorescence of two mirrored cell length ranges around the cell center.
    One value is returned, which corresponds to the average fluorescence of the pixels included,
    in the mirrored cell length sectors.
    """
    if cell_id in oned_coords_dict:
        oned_fluor_df = oned_coords_dict[cell_id]
        if type(cell_range)==float:
            sector_mean = oned_fluor_df[oned_fluor_df.scaled_length.between(-cell_range,cell_range)]['fluor_pixels_'+str(channel)].mean()
            sector_area = oned_fluor_df[oned_fluor_df.scaled_length.between(-cell_range,cell_range)]['fluor_pixels_'+str(channel)].count()
        elif type(cell_range)==tuple:
            sector_mean = oned_fluor_df[oned_fluor_df.scaled_length.between(cell_range[0],cell_range[1])|oned_fluor_df.scaled_length.between(-cell_range[1],-cell_range[0])]['fluor_pixels_'+str(channel)].mean()
            sector_area = oned_fluor_df[oned_fluor_df.scaled_length.between(cell_range[0],cell_range[1])|oned_fluor_df.scaled_length.between(-cell_range[1],-cell_range[0])]['fluor_pixels_'+str(channel)].count()
    else:
        sector_mean = np.nan
        sector_area = np.nan
    
    return sector_mean, sector_area


# The kernels below are used to determine the average timepoint of consistent nucleoid splitting,
# for nucleoids with consecutive splitting merging events of the dynamic nucleoid polymer. 
def distance_kernel_1(data_array):
       
    return data_array.rolling(window=7, center=True, min_periods=1).sum() 

def distance_kernel_2(data_array):
    bool_array = data_array>0
    bool_array = bool_array.astype(int)
  
    return bool_array.expanding().sum() 

def distance_kernel_3(data_array): # this kernel is not being used currently

    return data_array.rolling(window=3, center=True, min_periods=1).median()


def get_single_cell_asymmetries(cell_dataframe, results_path, 
                                asymmetry_list=[(0,0.25), (0.25,0.5), (0.5,0.75), (0.75,1), (0.5,1), (0,0.5), (0,1)],
                                symmetry_list=[(0.05),(0.10),(0.25),(0.5),(0.75),(0.25,0.5),(0.5,0.75),(0.75,1),(0.5,1),(0.25,1),(0.25,0.75)]):
    """
    Applies the functions:
        get_nucleoid_asymmetry
        get_nucleoid_distance
        get_nucleoid_distance
        get_fluorescence_asymmetry
        get_symmetric_mean
        distance_kernel_1
        distance_kernel_3
    """
    oned_coords_dict = get_oned_coords(results_path)
    
    cell_dataframe['object_asymmetry'] = cell_dataframe.object_scaled_length.apply(get_nucleoid_asymmetry)
    cell_dataframe['object_distance'] = cell_dataframe.object_centered_length.apply(get_nucleoid_distance)
    # apply the first kernel
    cell_dataframe['curated_object_distance'] = cell_dataframe.object_distance.copy()
    cell_dataframe['kernel_1'] = cell_dataframe.groupby('cell_trajectory_id').curated_object_distance.apply(distance_kernel_1)
    cell_dataframe.loc[cell_dataframe.curated_object_distance==cell_dataframe.kernel_1, 'curated_object_distance']=0
    # apply the second kernel
    cell_dataframe['kernel_2'] = cell_dataframe.groupby('cell_trajectory_id').curated_object_distance.apply(distance_kernel_2)
    cell_dataframe.loc[(cell_dataframe.curated_object_distance==0)&(cell_dataframe.kernel_2>0), 'curated_object_distance']=np.nan
    # remove kernel columns
    cell_dataframe = cell_dataframe.drop(['kernel_1', 'kernel_2'], axis=1)
    
    channel = 2
    for cell_range in asymmetry_list:
        print('asymmetric range:', cell_range)
        cell_dataframe['ribo_asymmetry'] = cell_dataframe.apply(lambda x: get_fluorescence_asymmetry(oned_coords_dict, x.cell_id, cell_range, channel), axis=1)
        cell_dataframe[['newp_ribo_'+str(cell_range), 'oldp_ribo_'+str(cell_range), 'newp_area_'+str(cell_range), 'oldp_area_'+str(cell_range)]] = pd.DataFrame(cell_dataframe.ribo_asymmetry.tolist(), index=cell_dataframe.index)
        cell_dataframe = cell_dataframe.drop(['ribo_asymmetry'], axis=1)
    
    for cell_range in symmetry_list:
#        if type(cell_range)==float:
#            cl_rg = (0,cell_range)
#        elif type(cell_range)==tuple:
#            cl_rg = cell_range
        print('symmetric range:', cell_range)
        cell_dataframe['ribo_symmetry'] = cell_dataframe.apply(lambda x: get_symmetric_mean(oned_coords_dict, x.cell_id, cell_range, channel), axis=1)
        cell_dataframe[['sym_ribo_'+str(cell_range), 'sym_area_'+str(cell_range)]] = pd.DataFrame(cell_dataframe.ribo_symmetry.tolist(), index=cell_dataframe.index)
        cell_dataframe = cell_dataframe.drop(['ribo_symmetry'], axis=1)
    
    channel = 3
    for cell_range in symmetry_list:
        print('symmetric range:', cell_range)
        cell_dataframe['nuc_symmetry'] = cell_dataframe.apply(lambda x: get_symmetric_mean(oned_coords_dict, x.cell_id, cell_range, channel), axis=1)
        cell_dataframe[['sym_nuc_'+str(cell_range), 'sym_area_2_'+str(cell_range)]] = pd.DataFrame(cell_dataframe.nuc_symmetry.tolist(), index=cell_dataframe.index)
        cell_dataframe = cell_dataframe.drop(['nuc_symmetry'], axis=1)

    return cell_dataframe   
   
    
        
# ------------------   9. NUCLEOID TRACKING  --------------------- #
def correct_nucleoid_positions(nucleoid_positions):
    """
    This function returns the average nucleoid position for cells with two nucleoid objects.
    The mid-point of the two nucleoid centroids is returned.
    """
    nuc_pos = np.array(nucleoid_positions)
    if nuc_pos[~np.isnan(nuc_pos)].shape[0] == 0:
        return nuc_pos
    elif nuc_pos[~np.isnan(nuc_pos)].shape[0] == 1:
        return nuc_pos
    elif nuc_pos[~np.isnan(nuc_pos)].shape[0] > 1:
        return np.array([np.nanmean(nuc_pos[nuc_pos<0]), np.nanmean(nuc_pos[nuc_pos>=0])])

def follow_the_nucleoid(oned_coords_dict, cell_id, nucleoid_positions, channel, px_wid):
    """
    This function links the nuceoids in the new or the old cell half.
    """
    if cell_id in oned_coords_dict:
        oned_df = oned_coords_dict[cell_id]
        nuc_pos = correct_nucleoid_positions(nucleoid_positions)
            
        if nuc_pos[~np.isnan(nuc_pos)].shape[0] == 0:
            return [np.nan, np.nan, np.nan]
        elif nuc_pos[~np.isnan(nuc_pos)].shape[0] == 1:
            try:
                nuc_df = oned_df[oned_df.arch_length.between(nuc_pos[0]-px_wid, nuc_pos[0]+px_wid)]
            except IndexError:
                nuc_df = oned_df[oned_df.arch_length.between(nuc_pos-px_wid, nuc_pos+px_wid)]
            return [nuc_df['fluor_pixels_'+str(channel)].mean(), np.nan, np.nan]
        elif nuc_pos[~np.isnan(nuc_pos)].shape[0] == 2:
            nuc_df = oned_df[oned_df.arch_length.between(nuc_pos[0]-px_wid, nuc_pos[0]+px_wid)|oned_df.arch_length.between(nuc_pos[1]-px_wid, nuc_pos[1]+px_wid)]
            mean_fluor = nuc_df['fluor_pixels_'+str(channel)].mean()
            oldp_fluor = oned_df[oned_df.arch_length.between(nuc_pos[0]-px_wid, nuc_pos[0]+px_wid)]['fluor_pixels_'+str(channel)].mean()
            newp_fluor = oned_df[oned_df.arch_length.between(nuc_pos[1]-px_wid, nuc_pos[1]+px_wid)]['fluor_pixels_'+str(channel)].mean()
            return [mean_fluor, newp_fluor, oldp_fluor]
    else:
        return [np.nan, np.nan, np.nan]


def apply_follow_the_nucleoid(oned_coords_dict, cell_dataframe, px_wid, nucleoid_center_column):
    # nucleoid_center_column = 'object_centered_length'
    # nucleoid_center_column = 'tracked_nucleoid_center'
    """
    This function applies the follow_the_nucleoid function.
    """
    for channel in range(2,4):
        cell_dataframe['nuc_follow_'+str(channel)] =  cell_dataframe.apply(lambda x: follow_the_nucleoid(oned_coords_dict, x.cell_id, x[nucleoid_center_column], channel, px_wid), axis=1)
        cell_dataframe[['nuc_follow_mean_'+str(channel), 'nuc_follow_new_'+str(channel), 'nuc_follow_old_'+str(channel)]] = pd.DataFrame(cell_dataframe['nuc_follow_'+str(channel)].tolist(), index=cell_dataframe.index)
        cell_dataframe = cell_dataframe.drop(['nuc_follow_'+str(channel)], axis=1)
       
    return cell_dataframe


def track_nucleoids(mother_df, cell_df, results_path, buffer_zone=4, kernel_width=11, kernel_rounds=2, center_smooth=11, show=False):
    """
    This function uses the linked nucleoid objects to track the nucleoids. 
    """
    def get_mother_trajectory(cell_trajectory, cell_dataframe):
        mother_cell = cell_dataframe[cell_dataframe.cell_trajectory_id==cell_trajectory].mother_cell.values[0]
        if mother_cell in cell_dataframe.cell_id.unique():
            return cell_dataframe[cell_dataframe.cell_id==mother_cell].cell_trajectory_id.values[0]
        else:
            return 'no cell here'
    
    def tracked_positions_in_mother(list_of_positions, moth_polarity, traj_polarity, buffer):
            nuc_pos = np.array(list_of_positions)
            if moth_polarity == -1:
                if traj_polarity == 1:
                    return nuc_pos[nuc_pos>buffer]
                elif traj_polarity == -1:
                    return nuc_pos[nuc_pos<-buffer]
            elif moth_polarity == 1:
                if traj_polarity == 1:
                    return nuc_pos[nuc_pos<-buffer]
                elif traj_polarity == -1:
                    return nuc_pos[nuc_pos>buffer]
        
    def correct_traj_positions(list_of_positions):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
        nuc_pos = np.array(list_of_positions)
        if nuc_pos[~np.isnan(nuc_pos)].shape[0] == 0:
            return nuc_pos[0]
        elif nuc_pos[~np.isnan(nuc_pos)].shape[0] == 1:
            return nuc_pos[0]
        elif nuc_pos[~np.isnan(nuc_pos)].shape[0] > 1:
            nuc_pos = np.array([np.nanmean(nuc_pos[nuc_pos<0]), np.nanmean(nuc_pos[nuc_pos>=0])])
            return np.nanmean(nuc_pos)
    
    def count_objects(list_of_objects):
        nuc_pos = np.array(list_of_objects)
        return nuc_pos.shape[0]
    
    f = 0
    i = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for traj in mother_df.cell_trajectory_id.unique():
            traj_df = mother_df[mother_df.cell_trajectory_id==traj]
            moth_df = cell_df[cell_df.cell_trajectory_id==get_mother_trajectory(traj, cell_df)]
            traj_df = traj_df[~traj_df.number_of_objects.isnull()]
            moth_df = moth_df[~moth_df.number_of_objects.isnull()]
            
            if moth_df.shape[0]>0 and traj_df.shape[0]>0:
                traj_polarity = traj_df.polarity.values[0]
                moth_polarity = moth_df.polarity.values[0]
                
                polarity_group = traj_polarity + moth_polarity
                if polarity_group == 0:
                    polarity_group = traj_polarity.copy()
                
                if moth_polarity in [1,-1] and traj_polarity in [1,-1]:
                    traj_df['tracked_nucleoid_positions'] = traj_df.object_centered_length
                    moth_df['tracked_nucleoid_positions'] = moth_df.apply(lambda x: tracked_positions_in_mother(x.object_centered_length, moth_polarity, traj_polarity, buffer_zone), axis=1)
                    traj_df['tracked_nucleoid_center'] = traj_df.object_centered_length.apply(correct_traj_positions)
                    moth_df['tracked_nucleoid_center'] = moth_df.tracked_nucleoid_positions.apply(np.nanmean)
                    
                    traj_df['tracked_number_of_objects'] = traj_df.number_of_objects
                    moth_df['tracked_number_of_objects'] = moth_df.tracked_nucleoid_positions.apply(count_objects).astype(int)
                    
                    for rnds in range(kernel_rounds):
                        traj_df['tracked_number_of_objects'] = traj_df.tracked_number_of_objects.rolling(kernel_width, center=True, min_periods=1).apply(lambda x: x.mode()[0])
#                        traj_df['tracked_number_of_objects'] = traj_df.tracked_number_of_objects.rolling(kernel_width, center=True, min_periods=1).apply(lambda x: x.mode()[0])
                        moth_df['tracked_number_of_objects'] = moth_df.tracked_number_of_objects.rolling(kernel_width, center=True, min_periods=1).apply(lambda x: x.mode()[0])
#                        moth_df['tracked_number_of_objects'] = moth_df.tracked_number_of_objects.rolling(kernel_width, center=True, min_periods=1).apply(lambda x: x.mode()[0])
                    
                    traj_df['tracked_nucleoid_center'] = traj_df['tracked_nucleoid_center'].rolling(center_smooth, center=True, min_periods=1).mean()
                    moth_df['tracked_nucleoid_center'] = moth_df['tracked_nucleoid_center'].rolling(center_smooth, center=True, min_periods=1).mean()
                    
                    min_nucleoid_div = moth_df[moth_df.tracked_number_of_objects==1].frame.min()
                    max_nucleoid_abs = moth_df[moth_df.tracked_number_of_objects==0].frame.max()
                    min_nucleoid_div_zer = moth_df[moth_df.tracked_number_of_objects==1].zeroed_frame.min()
                    if max_nucleoid_abs > min_nucleoid_div:
                        min_nucleoid_div = max_nucleoid_abs+1
                    min_nucleoid_div_2 = traj_df[traj_df.tracked_number_of_objects==2].frame.min()
                    min_nucleoid_div_zer_2 = traj_df[traj_df.tracked_number_of_objects==2].zeroed_frame.min()
    #                moth_df = moth_df[moth_df.zeroed_frame>=min_nucleoid_div]
    #                traj_df = traj_df[traj_df.zeroed_frame<min_nucleoid_div_2]
                    
                    nuc_df = pd.concat([moth_df, traj_df])[['cell_id', 'cell_trajectory_id', 
                                      'experiment', 'position_int', 'label',
                                      'frame', 'cc_phase', 'zeroed_frame', 'polarity', 
                                      'mean_fluor_2', 'std_fluor_2', 'cc_mean_fluor_2', 'cc_std_fluor_2',
                                      'mean_fluor_3', 'std_fluor_3', 'cc_mean_fluor_3', 'cc_std_fluor_3',
                                      'number_of_objects', 'object_asymmetry', 'object_distance', 'curated_object_distance',
                                      'tracked_nucleoid_positions', 'tracked_nucleoid_center', 'tracked_number_of_objects']]
    
                    nuc_df = nuc_df[nuc_df.frame.between(min_nucleoid_div, min_nucleoid_div_2-1)]
                    nuc_df['nucleoid_id']= moth_df.cell_trajectory_id.values[0]+'+'+traj_df.cell_trajectory_id.values[0]
                    nuc_df['zeroed_frame_2'] = nuc_df.frame - nuc_df.frame.min()
                    nuc_df['nuc_phase'] = nuc_df.zeroed_frame_2/nuc_df.zeroed_frame_2.max()
                    nuc_df['polarity_group'] = polarity_group
                    
                    if min_nucleoid_div_zer > 0 and nuc_df.shape[0]>0:
                        if show == True:
                            print(moth_polarity, traj_polarity, min_nucleoid_div_zer, min_nucleoid_div_zer_2)
                            ax1 = plt.subplot()
                            ax2 = ax1.twinx()
                            ax1.plot(nuc_df.frame, nuc_df.tracked_nucleoid_center)
                            ax2.plot(nuc_df.frame, nuc_df.tracked_number_of_objects)
                            plt.show()
                            input()
                        if i==0:
                            final_df = nuc_df
                        elif i > 0:
                            final_df = pd.concat([final_df, nuc_df])
    #                    print(i, f, final_df.shape[0])
                        i+=1
        #                input()
            f+=1
    
    if os.path.isdir(results_path):
        final_df.to_pickle(results_path+'/nucleoid_tracking_df')
    print(i/f*100, '% of cells with good nucleoid trajectories:',final_df.shape[0],'cell instances from',i,'nucleoid cycles.')
    return final_df
    


# ------------------   10. CALCULATION OF CELL DIVISION CYCLE AND INSTANTANEOUS RATES  --------------------- #
# Only works for 1min interval in the phase channel! 
def get_instantaneous_growth_rate(area_array):
    """
    This function returns the instantaneous relative grwoth rate from an array of cell areas.
    Exponential growth is assumed
    """
    log_area = np.log(area_array)
    time_h = np.arange(0, log_area.shape[0])/60 #interval is a global variable
    # interval is a gloval variable
    return np.polyfit(time_h, log_area, 1)[0]

def get_instantaneous_growth_rate_lin(area_array):
    """
    This function returns the instantaneous linear grwoth rate from an array of cell areas.
    """
#    log_area = np.log(area_array)
    time_h = np.arange(0, area_array.shape[0])/60 #interval is a global variable
    # interval is a gloval variable
    return np.polyfit(time_h, area_array, 1)[0]

def get_mean_area(area_array):
    """
    Returns the average area from an array of cell areas.
    """
    # interval is a gloval variable
    return np.mean(area_array)

def get_mean_cc(cc_array):
    """
    Returns the average area from an array of cell areas across the entire cell division cycle.
    """
    return np.mean(cc_array)

def get_area_increase(area_array):
    """
    Returns the added area from a sorted array.
    """
    area_list = list(area_array)
    return area_list[-1]-area_list[0]

def apply_instantaneous_rate(cell_dataframe, rolling_window, px_scale):
    """
    Estimates the instantaneous growth rates along a specified rolling window.
    """
    cell_dataframe = cell_dataframe.sort_values(['cell_trajectory_id', 'frame'])
    cell_dataframe['area_um'] = cell_dataframe.area_px*px_scale**2
    cell_dataframe['instant_log_gr'] = cell_dataframe.area_um.rolling(rolling_window, center=True).apply(get_instantaneous_growth_rate)
    cell_dataframe['instant_lin_gr'] = cell_dataframe.area_um.rolling(rolling_window, center=True).apply(get_instantaneous_growth_rate_lin)
    cell_dataframe['mean_area_um'] = cell_dataframe.area_um.rolling(rolling_window, center=True).apply(get_mean_area)
    cell_dataframe['mean_cc_phase'] = cell_dataframe.cc_phase.rolling(rolling_window, center=True).apply(get_mean_cc)

    cell_dataframe['cell_trajectory_code'] = cell_dataframe.cell_trajectory_id.astype('category').cat.codes
    cell_dataframe['cell_trajectory_av'] = cell_dataframe.cell_trajectory_code.rolling(rolling_window, center=True).mean()
    cell_dataframe.loc[cell_dataframe.cell_trajectory_av!=cell_dataframe.cell_trajectory_code, 'instant_log_gr']=np.nan
    cell_dataframe.loc[cell_dataframe.cell_trajectory_av!=cell_dataframe.cell_trajectory_code, 'instant_lin_gr']=np.nan
    cell_dataframe.loc[cell_dataframe.cell_trajectory_av!=cell_dataframe.cell_trajectory_code, 'mean_area_um']=np.nan
    cell_dataframe.loc[cell_dataframe.cell_trajectory_av!=cell_dataframe.cell_trajectory_code, 'mean_cc_phase']=np.nan

    cell_dataframe = cell_dataframe.drop(['cell_trajectory_av', 'cell_trajectory_code'], axis=1)
    return cell_dataframe


def get_cell_cycle_growth_rate(area_px, frames, px_scale, area_scale):
    """
    Estimates the cell division cycle averaged growth rates.
    """
    # ensure that the data are sorted
    gr_df = pd.DataFrame()
    gr_df['area'] = area_px*px_scale**2
    gr_df['time'] = frames/60
    gr_df = gr_df.sort_values('time')
    gr_df['area_scaled'] = gr_df.area/gr_df.area.tolist()[0]
    try:
        if area_scale == False:
            gr = np.polyfit(gr_df.time, np.log(gr_df.area), 1)[0]
        elif area_scale ==True:
            gr = np.polyfit(gr_df.time, np.log(gr_df.area_scaled), 1)[0]
    except:
        gr = np.nan
    return gr


def get_cc_gr_dict(cell_dataframe, px_scale):
    """
    Returns dictionaries of the cell division average growth rates (values),
    for each cell division cycle trajectory ID (keys).
    """
    cell_dataframe = cell_dataframe.sort_values(['cell_trajectory_id', 'frame'])
    gr_dict = cell_dataframe.groupby('cell_trajectory_id').apply(lambda x: get_cell_cycle_growth_rate(x.area_px, x.zeroed_frame, px_scale, False)).to_dict()
    gr_dict_scaled = cell_dataframe.groupby('cell_trajectory_id').apply(lambda x: get_cell_cycle_growth_rate(x.area_px, x.zeroed_frame, px_scale, True)).to_dict()
    return gr_dict, gr_dict_scaled


def apply_growth_rates(cell_dataframe, px_scale, rolling_window):
    """
    Applies the calculations of the instantaneous and cell division cycle average growth rates.
    """
    cell_dataframe = apply_instantaneous_rate(cell_dataframe, rolling_window, px_scale)
    gr_dict, gr_dict_scaled = get_cc_gr_dict(cell_dataframe, px_scale)
    cell_dataframe['cc_gr'] = cell_dataframe.cell_trajectory_id.map(gr_dict)
    cell_dataframe['cc_gr_scaled'] = cell_dataframe.cell_trajectory_id.map(gr_dict_scaled)
    return cell_dataframe


# The functions below can be applied for the estimation of instantaneous 
# and cell division cycle average growth rates for any quantified cell feature and not just the cell area.
def get_z_score(cell_dataframe, feature):
    
    mean_dict = cell_dataframe.groupby('cell_trajectory_id').mean()[feature]
    std_dict = cell_dataframe.groupby('cell_trajectory_id').std()[feature]
    
    cell_dataframe['mean_feat'] = cell_dataframe.cell_trajectory_id.map(mean_dict)
    cell_dataframe['std_feat'] = cell_dataframe.cell_trajectory_id.map(std_dict)
    cell_dataframe[feature+'_zscore'] = (cell_dataframe[feature]-cell_dataframe['mean_feat'])/cell_dataframe['std_feat']
    cell_dataframe = cell_dataframe.drop(['mean_feat','std_feat'], axis=1)
    
    return cell_dataframe

def get_instantaneous_rate(cell_dataframe, columns, rol_window, object_distance):
    
    if object_distance==True:
        cell_dataframe = cell_dataframe[cell_dataframe.curated_object_distance>0]
    
    for col in columns:
        cell_dataframe[col+'_rolwin'] = cell_dataframe[col].rolling(rol_window, center=True).mean()
    cell_dataframe['frame_rolwin'] = cell_dataframe['frame'].rolling(rol_window, center=True).mean()
    cell_dataframe['cell_trajectory_code'] = cell_dataframe.cell_trajectory_id.astype('category').cat.codes
    cell_dataframe['cell_trajectory_av'] = cell_dataframe.cell_trajectory_code.rolling(rol_window, center=True).mean()
    for col in columns:
        cell_dataframe.loc[cell_dataframe.cell_trajectory_av!=cell_dataframe.cell_trajectory_code, col+'_rolwin']=np.nan
    cell_dataframe.loc[cell_dataframe.cell_trajectory_av!=cell_dataframe.cell_trajectory_code, 'frame_rolwin']=np.nan
    cell_dataframe = cell_dataframe.drop(['cell_trajectory_av'], axis=1)
    
    for col in columns:
        cell_dataframe[col+'_dif'] = cell_dataframe[col+'_rolwin']-cell_dataframe[col+'_rolwin'].shift(1)
    cell_dataframe['frame_dif'] = cell_dataframe['frame_rolwin']-cell_dataframe['frame_rolwin'].shift(1)
    cell_dataframe['cell_trajectory_dif'] = cell_dataframe['cell_trajectory_code']-cell_dataframe['cell_trajectory_code'].shift(1)
    for col in columns:
        cell_dataframe.loc[cell_dataframe.cell_trajectory_dif!=0, col+'_dif']=np.nan
    cell_dataframe.loc[cell_dataframe.cell_trajectory_dif!=0, 'frame_dif']=np.nan
    cell_dataframe = cell_dataframe.drop(['cell_trajectory_dif', 'cell_trajectory_code'], axis=1)
    
    for col in columns:
        cell_dataframe[col+'_rate'] = cell_dataframe[col+'_dif']/cell_dataframe['frame_dif']/60
    
    return cell_dataframe

def get_cell_cycle_rate(variable_array, frames):
    
    gr_df = pd.DataFrame()
    gr_df['vari'] = variable_array*1
    gr_df['time'] = frames/60
    gr_df = gr_df.dropna()
    gr_df = gr_df.sort_values('time')
    try:
        gr = np.polyfit(gr_df.time, gr_df.vari, 1)[0]
    except:
        gr = np.nan
    return gr

def get_cc_rate_dict(cell_dataframe, variable):
    cell_dataframe = cell_dataframe.sort_values(['cell_trajectory_id', 'frame'])
    rate_dict = cell_dataframe.groupby('cell_trajectory_id').apply(lambda x: get_cell_cycle_rate(x[variable], x.zeroed_frame)).to_dict()
    return rate_dict

def get_nuc_rate_dict(cell_dataframe, variable):
    cell_dataframe = cell_dataframe.sort_values(['nucleoid_id', 'frame'])
    rate_dict = cell_dataframe.groupby('nucleoid_id').apply(lambda x: get_cell_cycle_rate(x[variable], x.zeroed_frame_2)).to_dict()
    return rate_dict
        


#-----------  11. DEFINITION OF CELL AGES  ----------#
def get_mother_cells_from_mother_trajectories(cell_df, mother_trajectories):
    """
    Returs the cell IDs of the pre-divisional mother cells of a list of trajectories.
    """
    list_of_mother_cells = []
    for traj in mother_trajectories:
        traj_df = cell_df[cell_df.cell_trajectory_id==traj]
        list_of_mother_cells.append(traj_df.sort_values('frame').cell_id.tolist()[-1])
    return list_of_mother_cells

def get_next_age(cell_df, mother_cells, mother_polarity, age):
    """
    Returns a list of the younger age.
    """
    if age%2==0:
        mother_polarity = mother_polarity * -1
    daughter_df = cell_df[(cell_df.polarity==mother_polarity)&(cell_df.mother_cell.isin(mother_cells))]
    return daughter_df.cell_trajectory_id.unique().tolist()

def get_ages(mother_df, initial_threshold, mother_polarity):
    """
    This function determines the ages of the tracked cell lineages.
    """
    age = 1
    
    mean_df = mother_df.groupby('cell_trajectory_id').mean()
    thres = filters.threshold_otsu(mean_df[mean_df.y>initial_threshold].y)
    print(thres)
    plt.hist(mean_df.y, bins=np.arange(200,500,2))
    plt.axvline(thres)
    plt.xlabel('Mean trajectory centroid (px)')
    plt.ylabel('Number of cell cycles')
    plt.show()
    
    if mother_polarity == -1: 
        initial_mothers_df = mean_df[(mean_df.y>thres)&(mean_df.polarity==mother_polarity)]
    elif mother_polarity == 1:
        initial_mothers_df = mean_df[(mean_df.y<thres)&(mean_df.polarity==mother_polarity)]
    
    # age 1
    list_of_mother_traj = initial_mothers_df.index.tolist()
    list_of_mother_cells = get_mother_cells_from_mother_trajectories(mother_df, list_of_mother_traj)
    more_mothers_list = mother_df[(mother_df.polarity==mother_polarity)&(mother_df.mother_cell.isin(list_of_mother_cells))].cell_trajectory_id.unique().tolist()
    list_of_mother_traj = list(set(list_of_mother_traj) | set(more_mothers_list))
    print(len(list_of_mother_traj), 'mother cells of age:',age)
    
    age+=1 # age 2
    list_of_d2_traj = get_next_age(mother_df, list_of_mother_cells, -1, age)
    list_of_d2_cells = get_mother_cells_from_mother_trajectories(mother_df, list_of_d2_traj)
    print(len(list_of_d2_traj), 'daughter cells of age:',age)
    d2_df = mother_df[mother_df.cell_trajectory_id.isin(list_of_d2_traj)]

    age+=1 # age 3
    list_of_d4_traj = get_next_age(mother_df, list_of_d2_cells, 1, age)
#    list_of_d4_cells = get_mother_cells_from_mother_trajectories(mother_df, list_of_d4_traj)
    print(len(list_of_d4_traj), 'daughter cells of age:',age)
    d4_df = mother_df[mother_df.cell_trajectory_id.isin(list_of_d4_traj)]
    
    age+=1 # age 4
    list_of_d3_traj = get_next_age(mother_df, list_of_d2_cells, 1, age)
    list_of_d3_cells = get_mother_cells_from_mother_trajectories(mother_df, list_of_d3_traj)
    print(len(list_of_d3_traj), 'daughter cells of age:',age)
    d3_df = mother_df[mother_df.cell_trajectory_id.isin(list_of_d3_traj)]
    
    age+=1 # age 5
    list_of_d5_traj = get_next_age(mother_df, list_of_d3_cells, -1, age)
#     list_of_d4_cells = get_mother_cells_from_mother_trajectories(mother_df, list_of_d4_traj)
    print(len(list_of_d5_traj), 'daughter cells of age:',age)
    d5_df = mother_df[mother_df.cell_trajectory_id.isin(list_of_d5_traj)]
    
    age+=1 # age 6
    list_of_d6_traj = get_next_age(mother_df, list_of_d3_cells, -1, age)
#     list_of_d4_cells = get_mother_cells_from_mother_trajectories(mother_df, list_of_d4_traj)
    print(len(list_of_d6_traj), 'daughter cells of age:',age)
    d6_df = mother_df[mother_df.cell_trajectory_id.isin(list_of_d6_traj)]
    
    
    plt.hist(initial_mothers_df.groupby('cell_trajectory_id').mean().y, bins=np.arange(200,500,4), label='age 1', density=True)
    plt.hist(d2_df.groupby('cell_trajectory_id').mean().y, bins=np.arange(200,500,4), label='age 2', density=True)
    plt.hist(d3_df.groupby('cell_trajectory_id').mean().y, bins=np.arange(200,500,4), label='age 3', density=True)
    plt.hist(d4_df.groupby('cell_trajectory_id').mean().y, bins=np.arange(200,500,4), label='age 4', density=True)
    plt.hist(d5_df.groupby('cell_trajectory_id').mean().y, bins=np.arange(200,500,4), label='age 5', density=True)
    plt.hist(d6_df.groupby('cell_trajectory_id').mean().y, bins=np.arange(200,500,4), label='age 6', density=True)
    plt.xlabel('Mean trajectory centroid (px)')
    plt.ylabel('Frequency of cell cycles')
    
    
    mean_df = initial_mothers_df.groupby('cell_trajectory_id').mean()
    mean_df['age']=1
    age1_dict = mean_df.age.to_dict()
    
    mean_df = d2_df.groupby('cell_trajectory_id').mean()
    mean_df['age']=2
    age2_dict = mean_df.age.to_dict()
    
    mean_df = d3_df.groupby('cell_trajectory_id').mean()
    mean_df['age']=3
    age3_dict = mean_df.age.to_dict()
    
    mean_df = d4_df.groupby('cell_trajectory_id').mean()
    mean_df['age']=4
    age4_dict = mean_df.age.to_dict()
    
    mean_df = d5_df.groupby('cell_trajectory_id').mean()
    mean_df['age']=5
    age5_dict = mean_df.age.to_dict()
    
    mean_df = d6_df.groupby('cell_trajectory_id').mean()
    mean_df['age']=6
    age6_dict = mean_df.age.to_dict()
    
    age_dict = {**age1_dict, **age2_dict, **age3_dict, **age4_dict, **age5_dict, **age6_dict}
    
    d7_df = mother_df[~mother_df.cell_trajectory_id.isin(age_dict)]
    mean_df = d7_df.groupby('cell_trajectory_id').mean()
    mean_df['age']=7
    plt.hist(d7_df.groupby('cell_trajectory_id').mean().y, bins=np.arange(200,500,4), label='age NA', density=True)
    age7_dict = mean_df.age.to_dict()
    age_dict = {**age_dict, **age7_dict}
    
    plt.legend()
    plt.show()
    
    return age_dict