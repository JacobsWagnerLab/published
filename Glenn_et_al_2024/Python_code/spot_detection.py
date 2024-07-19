# -*- coding: utf-8 -*-
"""
Code for detecting spot of MipZ in C.rescentus and extract properties of MipZ foci

@author: A.Fragasso
"""
# import modules
from scipy.signal import lfiltic, lfilter
import skimage
from scipy import io
import os
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import math
import tifffile as tf
import torch
import keyboard
import time
import cupy as cp
import json
import pickle
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import splprep, splev, splrep, interp1d, UnivariateSpline, BSpline
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import laplace
from scipy.signal import lfilter, savgol_filter
from scipy.stats import spearmanr
from skimage import graph, morphology, measure
import imageio.v3 as iio
import imageio
from scipy.ndimage import gaussian_filter, gaussian_filter1d, binary_fill_holes, binary_closing, binary_opening, binary_erosion, center_of_mass, label
from scipy import ndimage
from scipy.spatial import distance
from PIL import Image as im
from skimage.morphology import skeletonize
from skimage.morphology import binary_closing
from skimage.filters import threshold_otsu, threshold_local
from skimage.segmentation import find_boundaries
from skimage import measure, morphology
from numpy.lib.stride_tricks import as_strided
import warnings
warnings.filterwarnings('ignore')


runfile("path_to_analysis_functions_file", wdir="path_to_analysis_folder")

'''Write paths to folders wherer the data are stored'''
folder_path_list = [
    'path_to_experiment_1',
    'path_to_experiment_2',
    'path_to_experiment_3'
    ]


'''Assign names to experiments'''
exp_list = [
    'name_experiment_1',
    'name_experiment_2',
    'name_experiment_3'
        ]

'''Create list of positions to scan through'''
pos_list = create_pos_list(1, 10, double_digit = True) + create_pos_list(1, 10, double_digit = False)

count=0     # counter for number of cells

'''
Input global variables
'''
fr = 3     # [min] time between frames (write the specific one for each experiment)
px_size = 0.064 # px size
save = True
plot = True
random_cmap = generate_random_cmap(1000)


ch_list = ['phase', 'fluor1']  # input here the name of the fluorescence channel folders
ch_int = [1, 3]  # frame intervals [min]
ch_fluor_list = [s for s in ch_list if 'fluor' in s]


for s in range(len(folder_path_list)):
    experiment_path = folder_path_list[s]
    save_dict_path = experiment_path + '/output'
    save_dict_path_2 = experiment_path + '/output_new'
    file_names_list = os.listdir(save_dict_path)
    filename_0 = file_names_list[0]
    idx = filename_0.find('xy')
    exp_name = filename_0[:idx-1]
    print(exp_name)
    for pos in pos_list:
        folder_path = experiment_path + '/' +\
            pos  # path to position
        if not os.path.exists(folder_path):
            continue
        cells_path = folder_path+"/cell"
        mask_path = folder_path + '/masks'
        phase_path = folder_path + '/phase'
        fluor1_path = folder_path + '/fluor1'     # MipZ channel
        cells_list = os.listdir(cells_path)
        if len(cells_list) == 0:
            continue
        masks_list = os.listdir(mask_path)
        phase_list = os.listdir(phase_path)

        '''
        Import json and pickle files
        '''
        good_cells_list_name = [
            s for s in file_names_list if 'good_cell' in s and pos in s and exp_name in s]
        if len(good_cells_list_name) == 0:
            continue
        with open(save_dict_path+'/'+good_cells_list_name[0], 'r') as json_file:
            good_cells_list = json.load(json_file)
        with open(save_dict_path+'/'+exp_name+'_'+pos+'_cell_ch.pkl', 'rb') as pickle_file:
            cell_ch_dict = pickle.load(pickle_file)
        with open(save_dict_path+'/'+exp_name+'_'+pos+'_cell_info.pkl', 'rb') as pickle_file:
            cell_info_dict = pickle.load(pickle_file)
        files = os.listdir(save_dict_path_2)
        df_name = [s for s in files if 'screened' in s and pos in s]
        if len(df_name) == 0:
            continue
        with open(save_dict_path_2+'/'+df_name[0], 'rb') as f:
            df = pickle.load(f)
        cells_list = list(df['cell_id'].unique())
        paths_dict = {}
        file_list_dict = {}

        for ch in ch_list:
            paths_dict[ch] = folder_path + '/' + ch
            file_list_dict[ch] = os.listdir(paths_dict[ch])

        '''First create dictionary of images of particles, fluor_bkg, masks and phase'''

        fluor_bkg_sub_dict = {}
        mask_dict = {}
        phase_dict = {}
        particle_dict = {}
        centers_of_mass_dict = {}
        block_size = 501
        for ch in ch_fluor_list:
            print('..Subtract background for channel ' + ch +
                  ' at pos ' + pos + ' for experiment ' + exp_name)
            fluor_list = file_list_dict[ch]
            fluor_path = paths_dict[ch]
            for k in range(len(fluor_list)):
                mask_temp = iio.imread(mask_path+'/'+masks_list[k])
                phase_temp = iio.imread(phase_path+'/'+phase_list[k])
                mask_dict[k] = mask_temp
                phase_dict[k] = phase_temp
                if not k % ch_int[1] == 0:
                    continue
                fluor_temp = iio.imread(fluor_path + '/' + fluor_list[k])
                bs = int(fluor_temp.shape[0]/16)
                fluor_bkg_sub, fluor_bkg, fluor_median_bkg = local_bkg_sub_cp(
                    fluor_temp, mask_temp, pos, exp_name, frame=str(k), box_size=bs, dilations=15, sigma_=60, show=False)  # subtract background in fluorescence channel
                fluor_bkg_sub_dict[(ch, k)] = fluor_bkg_sub
                
                ''' Algorithm for spot detection'''
                blurred = gaussian_filter(fluor_bkg_sub, sigma=3)
                log_image = laplace(blurred)
                inverted_img = log_image - np.min(log_image)
                inverted_img = np.max(inverted_img)-inverted_img
                inverted_img = inverted_img/np.max(inverted_img)

                blurred_2 = gaussian_filter(inverted_img, sigma=4)
                otsu = threshold_otsu(blurred_2)
                # adaptive_thr = threshold_local(blurred_2, block_size, offset=-0.45)   # for 20240329
                adaptive_thr = threshold_local(
                    blurred_2, block_size, offset=-otsu/2.5)  # for 20240328
                binary_img = inverted_img > adaptive_thr
                mask_temp_cropped = mask_temp[:binary_img.shape[0],
                                              :binary_img.shape[1]]
                mask_labeled_img = (
                    binary_img*(mask_temp_cropped > 0)+mask_temp_cropped)
                mask_labeled_particle_img = mask_labeled_img*binary_img

                labeled_image, num_features = label(mask_labeled_particle_img)

                labeled_particles = (
                    labeled_image + mask_labeled_particle_img)*(labeled_image > 0)
                relabeled_image, num_features = measure.label(
                    labeled_particles, return_num=True, background=0)
                relabeled_particles = curate_small_objects(
                    relabeled_image, min_size=3)

                particle_dict[(ch, k)] = relabeled_particles
                centers_of_mass_dict[k] = center_of_mass(
                    relabeled_particles > 0, relabeled_particles, range(1, num_features + 1))

        df_new = pd.DataFrame()
        
        '''Scan through each cell'''
        for cell in cells_list:
            print('..Analyze cell ' + cell)
            birth, death, divide, motherID, sisterID, daughterID, gen, pole_1_label, pole_2_label, cell_type, p1_list, p2_list, xaxis_off_list = cell_info_dict[
                cell]
            cell_df = df[df['cell_id'] == cell]
            frame_list = cell_df['frame'].to_list()
            if plot:
                fig = plt.figure(figsize=(12, 9))
                fig.suptitle(cell + ' at position ' + pos + ', ' + exp_name)
            num_p_list = []
            dist_p_list = []
            mean_fluor1_particle_list = []
            max_fluor1_particle_list = []
            mean_log1_particle_list = []
            tot_particle1_list = []
            coord_particle1_list = []
            max_log_particle1_list = []
            area_particle1_list = []

            mean_fluor2_particle_list = []
            max_fluor2_particle_list = []
            mean_log2_particle_list = []
            tot_particle2_list = []
            coord_particle2_list = []
            max_log_particle2_list = []
            area_particle2_list = []
            edge_flag = False
            '''Scan through each frame of each cell'''
            for i in frame_list:
                '''Only consider frames were fluorescence image was acquired'''
                if not i % ch_int[1] == 0:
                    mean_fluor1_particle_list.append(np.nan)
                    max_fluor1_particle_list.append(np.nan)
                    mean_log1_particle_list.append(np.nan)
                    coord_particle1_list.append([np.nan, np.nan])
                    tot_particle1_list.append(np.median(np.nan))
                    max_log_particle1_list.append(np.nan)
                    area_particle1_list.append(np.nan)

                    mean_fluor2_particle_list.append(np.nan)
                    max_fluor2_particle_list.append(np.nan)
                    mean_log2_particle_list.append(np.nan)
                    coord_particle2_list.append([np.nan, np.nan])
                    tot_particle2_list.append(np.nan)
                    max_log_particle2_list.append(np.nan)
                    area_particle2_list.append(np.nan)
                    num_p_list.append(np.nan)
                    dist_p_list.append(np.nan)
                    continue
                else:
                    phase_temp = phase_dict[i]
                    mask_temp = mask_dict[i]
                    fluor1_temp = fluor_bkg_sub_dict[('fluor1', i)]

                    particles_img = particle_dict[('fluor1', i)]
                    found_p = False

                    cell_ch = cell_ch_dict[cell, i]
                    cell_stack_temp = cell_ch[2]
                    cropped_mask = cell_ch[0]

                    A, rcm, r_off, angle, x1, x2, y1, y2, edge_flag = get_cell_coordinates(
                        cell_stack_temp)
                    cell_mask = np.zeros(mask_temp.shape)
                    cell_mask[y1:y2, x1:x2] = cropped_mask

                    dims_fluor = particles_img.shape
                    if dims_fluor[1]-offs < x2 or dims_fluor[0]-offs < y2:
                        print(cell+' is on the edge and will be excluded')
                        edge_flag = True
                        break

                    particles_img_cropped = particles_img[y1:y2,
                                                          x1:x2]*cropped_mask
                    particles_img_cropped, num_features = measure.label(
                        particles_img_cropped, return_num=True, background=0)
                    fluor1_cropped = fluor1_temp[y1:y2, x1:x2]
                    
                    '''If no particles were found'''
                    if num_features == 0:
                        print('No particles at frame ' +
                              str(i) + ' for ' + cell)
                        num_p_list.append(0)
                        dist_p_list.append(np.nan)

                        mean_fluor1_particle_list.append(0)
                        max_fluor1_particle_list.append(0)
                        coord_particle1_list.append([np.nan, np.nan])
                        tot_particle1_list.append(0)
                        area_particle1_list.append(0)

                        mean_fluor2_particle_list.append(0)
                        max_fluor2_particle_list.append(0)
                        coord_particle2_list.append([np.nan, np.nan])
                        tot_particle2_list.append(0)
                        area_particle2_list.append(0)
                    
                    '''If one particle was found'''
                    elif num_features == 1:
                        found_p = True
                        regions = measure.regionprops(
                            particles_img_cropped, intensity_image=fluor1_cropped)
                        region1 = extract_region_properties(regions[0])
                        num_p_list.append(1)
                        dist_p_list.append(0)

                        mean_fluor1_particle_list.append(
                            region1['mean_intensity'])
                        max_fluor1_particle_list.append(
                            region1['max_intensity'])
                        area_particle1_list.append(region1['area'])
                        tot_particle1_list.append(
                            region1['mean_intensity']*region1['area'])
                        coord_particle1_list.append(region1['centroid'])

                        mean_fluor2_particle_list.append(0)
                        max_fluor2_particle_list.append(0)
                        coord_particle2_list.append([np.nan, np.nan])
                        tot_particle2_list.append(0)
                        area_particle2_list.append(0)

                    '''If more than one particle was found'''
                    elif num_features > 1:
                        found_p = True
                        '''If more than two particle were found, only keep the two with the largest area'''
                        if num_features > 2:

                            regions = measure.regionprops(
                                particles_img_cropped, intensity_image=fluor1_cropped)
                            sorted_regions = sorted(
                                regions, key=lambda region: region.area, reverse=True)
                            top_two_labels = [
                                sorted_regions[0].label, sorted_regions[1].label]
                            particles_img_cropped = np.where(
                                np.isin(particles_img_cropped, top_two_labels), particles_img_cropped, 0)

                        regions = measure.regionprops(
                            particles_img_cropped, intensity_image=fluor1_cropped)
                        region1 = extract_region_properties(regions[0])
                        region2 = extract_region_properties(regions[1])
                        centroid1 = np.array(region1['centroid'])
                        centroid2 = np.array(region2['centroid'])
                        relative_distance = distance.euclidean(
                            centroid1, centroid2)

                        num_p_list.append(2)
                        dist_p_list.append(relative_distance)

                        mean_fluor1_particle_list.append(
                            region1['mean_intensity'])
                        max_fluor1_particle_list.append(
                            region1['max_intensity'])
                        area_particle1_list.append(region1['area'])
                        tot_particle1_list.append(
                            region1['mean_intensity']*region1['area'])
                        coord_particle1_list.append(region1['centroid'])

                        mean_fluor2_particle_list.append(
                            region2['mean_intensity'])
                        max_fluor2_particle_list.append(
                            region2['max_intensity'])
                        area_particle2_list.append(region2['area'])
                        tot_particle2_list.append(
                            region2['mean_intensity']*region2['area'])
                        coord_particle2_list.append(region2['centroid'])

                    y1_n = np.max([y1-offs, 0])
                    y2_n = np.min([y2+offs, mask_temp.shape[0]])
                    x1_n = np.max([x1-offs, 0])
                    x2_n = np.min([x2+offs, mask_temp.shape[1]])

                    if plot:
                        if not i % ch_int[1] == 0:
                            continue
                        else:
                            contours = measure.find_contours(
                                cropped_mask, level=0.5)
                            contours_global = measure.find_contours(
                                cell_mask, level=0.5)

                            fig = plt.figure(figsize=(30, 20))
                            fig.suptitle(cell + ', ' + str(i) + ', ' +
                                         pos + ', ' + exp_name, fontsize=24)
                            ax1 = fig.add_subplot(2, 3, 1)
                            plt.title('phase')
                            plt.imshow(phase_temp, cmap='grey')
                            if contours_global != False:
                                for contour in contours_global:
                                    ax1.plot(contour[:, 1], contour[:, 0],
                                             linewidth=2, color='yellow')
                            plt.scatter(rcm[0], rcm[1],
                                        marker='x', color='red', s=104)
                            plt.xlim(x1_n, x2_n)
                            plt.ylim(y2_n, y1_n)

                            fig.add_subplot(2, 3, 2)
                            plt.title('masks')
                            plt.imshow(cropped_mask, cmap=random_cmap,
                                       interpolation='nearest')
                            if found_p:
                                for region in regions:
                                    plt.scatter(
                                        region['centroid'][1], region['centroid'][0], marker='x', color='red', s=180)
                            ax1 = fig.add_subplot(2, 3, 3)
                            plt.title('MipZ_fluor_bkg_sub')
                            plt.imshow(fluor1_temp)
                            if contours_global != False:
                                for contour in contours_global:
                                    ax1.plot(contour[:, 1], contour[:, 0],
                                             linewidth=2, color='yellow')
                            plt.xlim(x1_n, x2_n)
                            plt.ylim(y2_n, y1_n)

                            ax1 = fig.add_subplot(2, 3, 4)
                            plt.title('MipZ_fluor_bkg_sub cropped')
                            plt.imshow(fluor1_cropped)
                            if contours != False:
                                for contour in contours:
                                    ax1.plot(contour[:, 1], contour[:, 0],
                                             linewidth=2, color='yellow')

                            fig.add_subplot(2, 3, 5)
                            plt.title('Cell and particles masks')
                            plt.imshow(cropped_mask+particles_img_cropped,
                                       cmap=random_cmap, interpolation='nearest')

                            ax1 = fig.add_subplot(2, 3, 6)
                            plt.title('Particles masks')
                            plt.imshow(particles_img_cropped,
                                       cmap=random_cmap, interpolation='nearest')
                            if found_p:
                                for region in regions:
                                    plt.scatter(
                                        region['centroid'][1], region['centroid'][0], marker='x', color='red', s=180)
                            if contours != False:
                                for contour in contours:
                                    ax1.plot(contour[:, 1], contour[:, 0],
                                             linewidth=2, color='yellow')

                            plt.show()

            if edge_flag:    ## if cell was found on the edge --> exclude
                continue
            cell_df['num_p'] = num_p_list
            cell_df['dist_p'] = dist_p_list

            cell_df['mean_fluor_particle_1'] = mean_fluor1_particle_list
            cell_df['max_fluor_particle_1'] = max_fluor1_particle_list
            cell_df['area_particle_1'] = area_particle1_list
            cell_df['tot_fluor_particle_1'] = tot_particle1_list
            cell_df['centroid_particle_1'] = coord_particle1_list

            cell_df['mean_fluor_particle_2'] = mean_fluor2_particle_list
            cell_df['max_fluor_particle_2'] = max_fluor2_particle_list
            cell_df['area_particle_2'] = area_particle2_list
            cell_df['tot_fluor_particle_2'] = tot_particle2_list
            cell_df['centroid_particle_2'] = coord_particle2_list

            df_new = pd.concat([df_new, cell_df])

        if save:
            df_new_name = exp_name+'_'+pos+'_cell_features_particle_df'
            df_new.to_pickle(save_dict_path_2+'/'+df_new_name+'.pkl')
