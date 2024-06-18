# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 15:15:08 2024

@author: Alexandros Papagiannakis, Christine Jacobs-Wagner lab, Stanford University, 2024
"""

import os
import scipy as sp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import Biviriate_medial_axis_estimation as bma
import LoG_adaptive_image_filter as logad
import pickle
import skimage
import warnings
warnings.filterwarnings("ignore")
# import LoG_adaptive_image_filter as logad



class otsu_segmentation_class(object):
    
    def __init__(self, tif_path, experiment, channels, frame_interval, frame_range, every_nth, save_path):
        
        start = time.time()
        
        self.save_path = save_path
        self.experiment = experiment
        self.frame_interval = frame_interval
        
        xy_position_folders = os.listdir(tif_path)
        xy_position_folders = [x for x in xy_position_folders if experiment in x]
        xy_positions = []
        xy_position_number = len(xy_position_folders)
        xy_position_folder_dict = {}
        for pos in range(xy_position_number):
            xy_position_folder_dict[pos] = tif_path+'/'+[x for x in xy_position_folders if 'xy'+str(pos+1) in x][0]
            xy_positions.append('xy'+str(pos+1))
        
        self.xy_positions = xy_positions
        self.xy_position_number = xy_position_number
        self.xy_position_folder_dict = xy_position_folder_dict
        self.channels = channels
        self.channel_number = len(channels)
        
        image_arrays = {}
        
        for xy_pos in range(self.xy_position_number):
            image_arrays[xy_pos] = {}
            for ch in range(self.channel_number):
                image_arrays[xy_pos][channels[ch]] = {}
                channel_pos_directory = [x for x in os.listdir(self.xy_position_folder_dict[xy_pos]) if x[-5]==str(ch+1)]
                for t in range(len(channel_pos_directory)):
                    image_arrays[xy_pos][channels[ch]][t] = skimage.io.imread(self.xy_position_folder_dict[xy_pos]+'/'+channel_pos_directory[t])
        
        self.image_arrays = image_arrays
        
        end = time.time() 
        print(round(end-start,1), 'seconds to load the cropped images and masks from', 
              self.xy_position_number, 'xy positions,', self.channel_number, 
              'channels and', len(list(range(frame_range[0], frame_range[1]))), 'timepoints')
    
    
    def otsu_filter(self, phase_image, gaussian_level, filter_percentage, dilation_rounds, fill_holes = True, binary_closing = True):
        
        image = 1/phase_image
        smoothed_image = skimage.filters.gaussian(image, gaussian_level)
                                                 
        otsu_thresh =  skimage.filters.threshold_otsu((smoothed_image).ravel())
        
        otsu_image = smoothed_image >= otsu_thresh*filter_percentage
        
        dil = 0
        while dil < dilation_rounds:
            otsu_image = skimage.morphology.binary_dilation(otsu_image)
            dil+=1
            print(dil, 'dilation rounds')
            
        if fill_holes:
            otsu_image = sp.ndimage.binary_fill_holes(otsu_image)
            print('filling holes in masks')
            
        if binary_closing:
            otsu_image = skimage.morphology.binary_closing(otsu_image)
            print('smoothing the masks')
        
        return otsu_image
    
    
    def segment_phase_images(self, gaussian_level, otsu_fraction, dilation_rounds,
                             hole_fill, closing):
        
        phase_otsu_masks = {}
        
        for xy_pos in self.image_arrays:
            phase_otsu_masks[xy_pos]= {}
            for tm in self.image_arrays[xy_pos]['Phase']:
                image_labels = skimage.measure.label(self.otsu_filter(self.image_arrays[xy_pos]['Phase'][tm], 
                                                                      gaussian_level, otsu_fraction, dilation_rounds, hole_fill, closing))
                phase_otsu_masks[xy_pos][tm] = image_labels
                print(str(np.max(image_labels.ravel()))+' cells in position xy'+str(xy_pos+1)+',  in timepoint '+str(tm))
        self.phase_otsu_masks = phase_otsu_masks
        
        with open(self.save_path+'/'+self.experiment+'_otsu_masks_dict', 'wb') as handle:
            pickle.dump(phase_otsu_masks, handle)

    
    def load_otsu_masks(self):
        with open(self.save_path+'/'+self.experiment+'_otsu_masks_dict', 'rb') as handle:
            self.phase_otsu_masks = pickle.load(handle)
    
    
    def get_sinuosity(self, medial_axis_df):
        
        min_df = medial_axis_df[medial_axis_df.arch_length_centered==medial_axis_df.arch_length_centered.min()]
        max_df = medial_axis_df[medial_axis_df.arch_length_centered==medial_axis_df.arch_length_centered.max()]
        length = min_df.arch_length_centered.abs().values[0]+max_df.arch_length_centered.values[0]
        
        eucli_distance = np.sqrt((max_df.cropped_x.values[0] - min_df.cropped_x.values[0])**2
                                 +(max_df.cropped_y.values[0]-min_df.cropped_y.values[0])**2)
        
        return length - eucli_distance
        
        
    def get_cell_id_string(self, position, timepoint, mask_label):
        
        return 'xy'+str(position+1)+'_time'+str(timepoint)+'_label'+str(mask_label)
        
    
    def get_cropping_pad(self, masked_image, pad_frame=5):
    
        nonzero_coords = np.nonzero(masked_image)
        min_x = np.min(nonzero_coords[1])-pad_frame
        max_x = np.max(nonzero_coords[1])+pad_frame
        min_y = np.min(nonzero_coords[0])-pad_frame
        max_y = np.max(nonzero_coords[0])+pad_frame    
        
        return min_x, min_y, max_x, max_y
    
    
    def label_curation(self):
        
        self.load_otsu_masks()
        
            
        pad_coords = {}
        cropped_centroid_coords = {}
        cropped_masks = {}
        medial_axis_dict = {}
        max_distance = {}
        sinuosity = {}
        cell_area = {}
        bad_cells = []
        
        for xy_pos in self.phase_otsu_masks:
            for lbl in np.unique(self.phase_otsu_masks[xy_pos][0])[1:]:
                
                masked_image = self.phase_otsu_masks[xy_pos][0]==lbl
                min_x, min_y, max_x, max_y = self.get_cropping_pad(masked_image)
                # print(xy_pos, lbl, min_x, min_y, max_x, max_y)
                if min_x >0 and min_y >0 and max_x<2048 and max_y<2048:
                    # cell_id = 'xy'+str(xy_pos+1)+'_label'+str(lbl)
                    cell_id = self.get_cell_id_string(xy_pos, 0, lbl)
                    print(cell_id)
                    pad_coords[cell_id] = (min_x, min_y, max_x, max_y)
                    
                    cropped_masked_image  = masked_image[min_y:max_y+1, min_x:max_x+1]
                    
                    cropped_masks[cell_id] = cropped_masked_image
                    cropped_centroid_coords[cell_id] = skimage.measure.centroid(cropped_masked_image)
                    medial_axis_dict[cell_id] = bma.get_medial_axis(cropped_masked_image, radius_px=8, half_angle=22, cap_knot=13, max_degree=60, verbose=True)
                    
                    medial_axis_df = medial_axis_dict[cell_id][0]
                    cell_mask_df = bma.get_oned_coordinates(cropped_masked_image, medial_axis_df, half_window=5)
                    
                    max_distance[cell_id] = cell_mask_df.width.abs().max()
                    sinuosity[cell_id]= self. get_sinuosity(medial_axis_df)
                    cell_area[cell_id] = cropped_masked_image[np.nonzero(cropped_masked_image)].shape[0]
                    print('max distance',max_distance[cell_id])
                    print('sinuosity',sinuosity[cell_id])
                    print('cell area', cell_area[cell_id])
                    
                    if max_distance[cell_id] > 10:
                        bad_cells.append(cell_id)
                        print(cell_id, 'is a bad cell segmentation')
                        # input()
                    if sinuosity[cell_id] > 1:
                        if cell_id not in bad_cells:
                            bad_cells.append(cell_id)
                            print(cell_id, 'is a bad cell segmentation')
                            # input()
            
        
        curation_df = pd.DataFrame()
        curation_df['cell_id'] = list(sinuosity.keys())
        curation_df['sinuosity'] = curation_df.cell_id.map(sinuosity)
        curation_df['max_distance'] = curation_df.cell_id.map(max_distance)
        curation_df['cell_area'] = curation_df.cell_id.map(cell_area)
        
        curation_df.to_pickle(self.save_path+'/'+self.experiment+'_curation_statistics_df', compression='zip')
        print(bad_cells, 'are possibly bad cells')
        
        
    def get_good_cells(self, sinuosity_range, max_distance_range, cell_area_range):
        
        curation_df = pd.read_pickle(self.save_path+'/'+self.experiment+'_curation_statistics_df', compression='zip')
        
        plt.figure(figsize=(5,5))
        plt.scatter(curation_df.sinuosity*0.066, curation_df.max_distance*0.066, c=curation_df.cell_area*0.066**2)
        cbar = plt.colorbar()
        plt.yscale('log')
        plt.xscale('log')
        
        good_df = curation_df.copy()
        
        if sinuosity_range != 'none':
            good_df = good_df[good_df.sinuosity.between(sinuosity_range[0], sinuosity_range[1])]
            
            plt.axvline(sinuosity_range[0]*0.066, color='black', linewidth=0.5)
            plt.axvline(sinuosity_range[1]*0.066, color='black', linewidth=0.5)
        
        if max_distance_range != 'none':
            good_df = good_df[good_df.max_distance.between(max_distance_range[0], max_distance_range[1])]
            
            plt.axhline(max_distance_range[0]*0.066, color='black', linewidth=0.5)
            plt.axhline(max_distance_range[1]*0.066, color='black', linewidth=0.5)
        
        if cell_area_range != 'none':
            good_df = good_df[good_df.cell_area.between(cell_area_range[0], cell_area_range[1])]
            cbar.set_ticks([2,4,6,8,10,12,14]+[cell_area_range[0]*0.066**2, cell_area_range[1]*0.066**2])
        
        cbar.set_label('Cell area (um$^2$)', rotation=90)
        
        plt.xlabel('Sinuosity (um)')
        plt.ylabel('Max half width (um)')
        plt.show()
        
        good_cells = good_df.cell_id.tolist()
        bad_cells = curation_df[~curation_df.cell_id.isin(good_df.cell_id)].cell_id.to_list()
            
        with open(self.save_path+'/'+self.experiment+'_good_cells_list', 'wb') as handle:
            pickle.dump(good_cells, handle)
        
        with open(self.save_path+'/'+self.experiment+'_bad_cells_list', 'wb') as handle:
            pickle.dump(bad_cells, handle)
            
    
    def load_good_cells(self):
        
        with open(self.save_path+'/'+self.experiment+'_good_cells_list', 'rb') as handle:
            self.good_cells = pickle.load(handle)
            
            
    def load_bad_cells(self):
        
        with open(self.save_path+'/'+self.experiment+'_bad_cells_list', 'rb') as handle:
            self.good_cells = pickle.load(handle)
    
    
    def get_centroid_dataframe(self):
        
        self.load_otsu_masks()
        
        position_list = []
        time_list = []
        label_list = []
        x_list = []
        y_list = []
        area_list = []
        cell_id_list = []
        
        for xy_pos in self.phase_otsu_masks:
            print('Getting centroids for xy position', xy_pos+1)
            for tm in self.phase_otsu_masks[xy_pos]:
                
                large_cell_masks = skimage.morphology.remove_small_objects(self.phase_otsu_masks[xy_pos][tm]>0, min_size=300)
                cell_labels = self.phase_otsu_masks[xy_pos][tm]*large_cell_masks
                
                for lbl in np.unique(cell_labels)[1:]:
                    
                    whole_mask = self.phase_otsu_masks[xy_pos][tm]==lbl
                    centroid = skimage.measure.centroid(whole_mask)
                    x_list.append(centroid[1])
                    y_list.append(centroid[0])
                    position_list.append(xy_pos)
                    time_list.append(tm)
                    label_list.append(lbl)
                    area_list.append(np.nonzero(whole_mask)[0].shape[0])
                    cell_id_list.append(self.get_cell_id_string(xy_pos, tm, lbl))
        
        mask_df = pd.DataFrame()
        mask_df['cell_id'] = cell_id_list
        mask_df['position'] = position_list
        mask_df['x'] = x_list
        mask_df['y'] = y_list
        mask_df['cell_area'] = area_list
        mask_df['cell_label'] = label_list
        mask_df['t'] = time_list
        
        mask_df.to_pickle(self.save_path+'/'+self.experiment+'_centroids_df', compression='zip')
       
        
    def load_centroids_dataframe(self):
        
        self.centroid_df = pd.read_pickle(self.save_path+'/'+self.experiment+'_centroids_df', compression='zip')
        
    
    def link_cells(self, max_radius, max_area_per):
        
        self.load_centroids_dataframe()
        
        cell_linkage_dict = {}
        
        for xy_pos in self.centroid_df.position.unique():
            print('Linking cells in position XY'+str(xy_pos+1))
            pos_df = self.centroid_df[self.centroid_df.position==xy_pos]
            
            for cell_index, cell_row in pos_df.iterrows():
                cell_x = cell_row.x
                cell_y = cell_row.y
                cell_frame = cell_row.t
                cell_id_string = cell_row.cell_id
                cell_area_px = cell_row.cell_area
                
                if cell_frame < pos_df.t.max():
                    
                    temporary_df = pos_df[pos_df.t==cell_frame+1]
                    # temporary_df = pos_df.copy()
                    # temporary_df = temporary_df.drop(index=cell_index)
                    temporary_df['distance'] = np.sqrt((temporary_df.x-cell_x)**2 + (temporary_df.y-cell_y)**2) 
                    temporary_df['area_difference'] = temporary_df.cell_area- cell_area_px
                    
                    minimum_distance_df = temporary_df[temporary_df['distance']<= max_radius]
                    minimum_distance_df = minimum_distance_df[minimum_distance_df.area_difference.between(max_area_per[0]*cell_area_px, max_area_per[1]*cell_area_px)]
                    
                    if  minimum_distance_df.shape[0] == 1: 
                        cell_linkage_dict[cell_id_string] = minimum_distance_df.cell_id.values[0]
                    elif minimum_distance_df.shape[0] > 1:
                        cell_linkage_dict[cell_id_string] = minimum_distance_df[minimum_distance_df.distance == minimum_distance_df.distance.min()].cell_id.values[0]
                    elif minimum_distance_df.shape[0] == 0:
                        cell_linkage_dict[cell_id_string] = 'none'
                    
        with open(self.save_path+'/'+self.experiment+'_cell_linkage_dict', 'wb') as handle:
            pickle.dump(cell_linkage_dict, handle)
        
    
    def load_cell_linkage_dictionary(self):
        
        with open(self.save_path+'/'+self.experiment+'_cell_linkage_dict', 'rb') as handle:
            self.cell_linkage_dict = pickle.load(handle)
        
    
    def connect_cells_into_lineages(self):
        
        self.load_cell_linkage_dictionary()
        
        cell_trajectories_dictionary = {}
        
        traj = 0
        
        for cl in self.cell_linkage_dict:
            if cl not in cell_trajectories_dictionary:
                cell_trajectories_dictionary[cl] = 'traj_'+str(traj)
                if self.cell_linkage_dict[cl]!='none':
                    cell_trajectories_dictionary[self.cell_linkage_dict[cl]] = cell_trajectories_dictionary[cl]
            else:
                if self.cell_linkage_dict[cl]!='none':
                    cell_trajectories_dictionary[self.cell_linkage_dict[cl]] = cell_trajectories_dictionary[cl]
            traj+=1
        
        with open(self.save_path+'/'+self.experiment+'_cell_lineages_dict', 'wb') as handle:
            pickle.dump(cell_trajectories_dictionary, handle)
      
    
    def load_cell_lineages(self):
        
        with open(self.save_path+'/'+self.experiment+'_cell_lineages_dict', 'rb') as handle:
            self.cell_trajectories_dictionary = pickle.load(handle)
            
    
    def apply_cell_trajectories_to_centroid_dataframe(self):
        
        self.load_centroids_dataframe()
        self.load_cell_lineages()
        
        self.centroid_df['trajectory'] = self.centroid_df.cell_id.map(self.cell_trajectories_dictionary)
    
    
    def get_long_trajectories(self, min_trajectory_length, centroid_df):
        
        trajectory_length_dict = centroid_df.groupby('trajectory').count().cell_id.to_dict()
        centroid_df['trajectory_length'] = centroid_df.trajectory.map(trajectory_length_dict)
        
        return centroid_df[centroid_df.trajectory_length>=min_trajectory_length].trajectory.to_list()
    
    
    def get_cell_contour(self, cropped_cell_mask):
        
        cell_mesh_y, cell_mesh_x = zip(*skimage.measure.find_contours(cropped_cell_mask, level=0.5)[0])
        
        return cell_mesh_x, cell_mesh_y
    
    
    def load_manually_curated_bad_cells(self, exclude_double_cells=True, exclude_long_cells=True):
        
        bad_cells = pd.read_excel(self.save_path+'/Manual_curation_trajectories.xlsx')
        bad_trajectories_1 = np.array(bad_cells.Bad_trajectories)[np.nonzero(np.array(bad_cells.Bad_trajectories))]
        
        if exclude_double_cells:    
            bad_trajectories_2 = np.array(bad_cells.Double_cells)[np.nonzero(np.array(bad_cells.Double_cells))]
        else:
            bad_trajectories_2 = np.array([])
        
        if exclude_long_cells == True:
            bad_trajectories_3 = np.array(bad_cells.Very_long)[np.nonzero(np.array(bad_cells.Very_long))]
        else:
            bad_trajectories_3 = np.array([])
        bad_trajectories = list(bad_trajectories_1)+list(bad_trajectories_2)+list(bad_trajectories_3)
        
        bad_trajectories_string = []
        for traj in bad_trajectories:
            bad_trajectories_string.append('traj_'+str(traj))
        
        return bad_trajectories_string
    
    
    def check_segmentation(self, min_trajectory_length, check_curation=False, save_folder='none'):
        
        
        if self.experiment+'_medial_axis_dict' in os.listdir(self.save_path):
            self.load_medial_axis()
            got_medial_axis = True
            
            if self.experiment+'_bad_medial_axis_list' in os.listdir(self.save_path):
                self.load_bad_medial_axis()
            else:
                self.bad_medial_axis_list = []
            
        else:
            got_medial_axis = False
            
        self.apply_cell_trajectories_to_centroid_dataframe()
        # self.load_good_cells()
        self.load_otsu_masks()
        
        if check_curation == True:
            bad_trajectories = self.load_manually_curated_bad_cells()
        else:
            bad_trajectories = []
        

        cen_df = self.centroid_df.copy()
        long_trajectories = self.get_long_trajectories(min_trajectory_length, cen_df)
        
        cen_df = cen_df[cen_df.trajectory.isin(long_trajectories)]
        cen_df = cen_df[~cen_df.trajectory.isin(bad_trajectories)]
        print(cen_df.trajectory.unique().shape[0], 'long trajectories')
        
        for traj in cen_df.trajectory.unique():
            traj_df = cen_df[cen_df.trajectory==traj]
            traj_df = traj_df.sort_values('t')
            print(traj)
            for cell_index, cell_row in traj_df.iterrows():
                masked_image = self.phase_otsu_masks[cell_row.position][cell_row.t]==cell_row.cell_label
                min_x, min_y, max_x, max_y = self.get_cropping_pad(masked_image)

                # print(xy_pos, lbl, min_x, min_y, max_x, max_y)
                if min_x >0 and min_y >0 and max_x<2048 and max_y<2048:
                    # cell_id = 'xy'+str(xy_pos+1)+'_label'+str(lbl)
                    cropped_masked_image  = masked_image[min_y:max_y+1, min_x:max_x+1]
                    cropped_phase_image = self.image_arrays[cell_row.position]['Phase'][cell_row.t][min_y:max_y+1, min_x:max_x+1]
                    plt.imshow(cropped_phase_image, cmap='gray')
                    cell_mesh_x, cell_mesh_y = self.get_cell_contour(cropped_masked_image)
                    plt.plot(cell_mesh_x, cell_mesh_y, color='yellow', linewidth=0.5)
                    plt.text(5,5, traj+' '+str(cell_row.t*self.frame_interval)+'min', fontsize=12, fontweight='bold')
                    plt.xticks([])
                    plt.yticks([])
                    
                    if got_medial_axis == True:
                        if cell_row.cell_id in self.medial_axis_dict and cell_row.cell_id not in self.bad_medial_axis_list:
                            medial_axis_df = self.medial_axis_dict[cell_row.cell_id][0]
                            plt.plot(medial_axis_df.cropped_x, medial_axis_df.cropped_y, color='red', linewidth=1)
                    
                    if os.path.isdir(save_folder):
                        plt.savefig(save_folder+'/'+traj+'_'+str(cell_row.t)+'.jpeg')
                    plt.close()
            # print('end of trajectory. Press ENTER to continue...')
            # input()
        
        
    def draw_medial_axes(self, min_trajectory_length=18, half_angle_px=12, curated=True):
        
        medial_axis_dict = {}
        oned_coords_dict = {}
        
        self.apply_cell_trajectories_to_centroid_dataframe()
        self.load_otsu_masks()
        
        if curated == True:
            bad_trajectories = self.load_manually_curated_bad_cells(False, False)
        else:
            bad_trajectories = []
        

        cen_df = self.centroid_df.copy()
        long_trajectories = self.get_long_trajectories(min_trajectory_length, cen_df)
        
        cen_df = cen_df[cen_df.trajectory.isin(long_trajectories)]
        cen_df = cen_df[~cen_df.trajectory.isin(bad_trajectories)]
        print(cen_df.trajectory.unique().shape[0], 'long trajectories')
        
        for traj in cen_df.trajectory.unique():
            traj_df = cen_df[cen_df.trajectory==traj]
            traj_df = traj_df.sort_values('t')
            print(traj)
            for cell_index, cell_row in traj_df.iterrows():
                masked_image = self.phase_otsu_masks[cell_row.position][cell_row.t]==cell_row.cell_label
                min_x, min_y, max_x, max_y = self.get_cropping_pad(masked_image)

                # print(xy_pos, lbl, min_x, min_y, max_x, max_y)
                if min_x >0 and min_y >0 and max_x<2048 and max_y<2048:
                    # cell_id = 'xy'+str(xy_pos+1)+'_label'+str(lbl)
                    cropped_masked_image  = masked_image[min_y:max_y+1, min_x:max_x+1]
                    medial_axis_dict[cell_row.cell_id] = bma.get_medial_axis(cropped_masked_image, radius_px=8, half_angle=half_angle_px, cap_knot=13, max_degree=60, verbose=False)
                    oned_coords_dict[cell_row.cell_id] = bma.get_oned_coordinates(cropped_masked_image, medial_axis_dict[cell_row.cell_id][0], half_window=5)
                                       
        with open(self.save_path+'/'+self.experiment+'_medial_axis_dict', 'wb') as handle:
            pickle.dump(medial_axis_dict, handle)
    
        with open(self.save_path+'/'+self.experiment+'_oned_coordinates_dict', 'wb') as handle:
            pickle.dump(oned_coords_dict, handle)
            
            
    def load_medial_axis(self):
        
        with open(self.save_path+'/'+self.experiment+'_medial_axis_dict', 'rb') as handle:
            self.medial_axis_dict = pickle.load(handle)
            
        with open(self.save_path+'/'+self.experiment+'_oned_coordinates_dict', 'rb') as handle:
            self.oned_coords_dict = pickle.load(handle)
        
    
    def get_cell_length(self, medial_axis_df):
        
        min_arch = float(medial_axis_df[medial_axis_df.arch_length_scaled==-1].arch_length.abs())
        max_arch = float(medial_axis_df[medial_axis_df.arch_length_scaled==1].arch_length.abs())
        
        return min_arch + max_arch
    
    
    def curate_medial_axes(self):
        
        self.apply_cell_trajectories_to_centroid_dataframe()
        cen_df = self.centroid_df.copy()
        self.load_medial_axis()
        
        cells_with_medial_axis = list(self.medial_axis_dict.keys())
        cen_df = cen_df[cen_df.cell_id.isin(cells_with_medial_axis)]
        
        time_list = []
        length_list = []
        area_list = []
        sinuosity_list = []
        max_distance_list = []
        cell_id_list = []
        trajectory_list = []
        
        for  cell_index, cell_row in cen_df.iterrows():
            medial_axis_length = self.get_cell_length(self.medial_axis_dict[cell_row.cell_id][0]) 
            time_list.append(cell_row.t)
            length_list.append(medial_axis_length)
            area_list.append(cell_row.cell_area)
            sinuosity_list.append(self.get_sinuosity(self.medial_axis_dict[cell_row.cell_id][0]))
            max_distance_list.append(self.oned_coords_dict[cell_row.cell_id].width.max())
            
            cell_id_list.append(cell_row.cell_id)
            trajectory_list.append(cell_row.trajectory)
        
        curation_df = pd.DataFrame()
        curation_df['trajectory'] = trajectory_list
        curation_df['cell_id'] = cell_id_list
        curation_df['time'] = time_list
        curation_df['cell_area'] = area_list
        curation_df['sinuosity'] = sinuosity_list
        curation_df['max_distance'] = max_distance_list
        curation_df['cell_length'] = length_list
        
        curation_df.to_pickle(self.save_path+'/'+self.experiment+'_medial_axis_statistics_df', compression='zip')
        
        bad_medial_axis = curation_df[curation_df.sinuosity>1].cell_id.tolist()
        
        with open(self.save_path+'/'+self.experiment+'_bad_medial_axis_list', 'wb') as handle:
            pickle.dump(bad_medial_axis, handle)
        
        print(len(bad_medial_axis), 'cells with a bad medial axis definition, out of',len(list(self.medial_axis_dict.keys())), 'cells')
        
        
        
    def load_bad_medial_axis(self):
        
        self.medial_axis_curation_df = pd.read_pickle(self.save_path+'/'+self.experiment+'_medial_axis_statistics_df', compression='zip')
        
        try:
            with open(self.save_path+'/'+self.experiment+'_bad_medial_axis_list', 'rb') as handle:
                self.bad_medial_axis_list = pickle.load(handle)
        except FileNotFoundError:
            self.bad_medial_axis_list = []
        
        print(len(self.bad_medial_axis_list), 'cells with a bad medial axis definition loaded...')

    
    def apply_oned_fluorescence(self):
        
        demograph_data = {}
        
        self.apply_cell_trajectories_to_centroid_dataframe()
        cen_df = self.centroid_df.copy()
        self.load_medial_axis()
        self.load_bad_medial_axis()
        self.load_otsu_masks()
        
        for cell_index, cell_row in cen_df.iterrows():
            
            if cell_row.cell_id in self.oned_coords_dict and cell_row.cell_id not in self.bad_medial_axis_list:
                
                xy_pos = cell_row.position
                tm = cell_row.t
                lbl = cell_row.cell_label
                
                if self.image_arrays[xy_pos][self.channels[1]][tm].max() > 0:
                    
                    cell_mask_df = self.oned_coords_dict[cell_row.cell_id]
                    masked_image = self.phase_otsu_masks[xy_pos][tm]==lbl
                    min_x, min_y, max_x, max_y = self.get_cropping_pad(masked_image)
                    
                    for ch in self.channels:

                        cropped_signal_image = self.image_arrays[xy_pos][ch][tm][min_y:max_y+1, min_x:max_x+1]
                        cropped_masked_image  = masked_image[min_y:max_y+1, min_x:max_x+1]
                        cell_mask_df[ch+'_fluor'] = cropped_signal_image[np.nonzero(cropped_masked_image)]
                        
                        demograph_data[cell_row.cell_id] = cell_mask_df
        
        
        with open(self.save_path+'/'+self.experiment+'_demograph_fluorescence_dict', 'wb') as handle:
            pickle.dump(demograph_data, handle)
    
    
    def load_oned_fluorescence(self):
        
        with open(self.save_path+'/'+self.experiment+'_demograph_fluorescence_dict', 'rb') as handle:
            self.oned_fluorescence_dict = pickle.load(handle)

    
    def segment_fluorescence_accumulation(self, signal_image, log_parameters, min_area):
        
        # [2,1000,75,4,7,-2,2] for RplA-GFP
        # [2,1000,75,4,7,-1,2] for HupA-mCherry
        
        # signal_image = signal_image - np.min(signal_image)
        # signal_image = signal_image / np.max(signal_image) *100
        
        masked_fluorescence = logad.log_adaptive_filter(signal_image, log_parameters)[2]
        masked_fluorescence = skimage.morphology.remove_small_objects(masked_fluorescence, min_area)
        masked_fluorescence = sp.ndimage.binary_fill_holes(masked_fluorescence)
        
        return masked_fluorescence
    
    
    def apply_fluorescence_segmentation(self, channel, log_parameters, min_area):
        
        segmentation_dict = {}
        print('segmenting objects in the',channel,'channel')
        for pos in self.image_arrays:
            segmentation_dict[pos] = {}
            for tm in self.image_arrays[pos][channel]:
                img = self.image_arrays[pos][channel][tm]
                if np.mean(img) > 0:
                    print('position:',pos,', timepoint:',tm)
                    segmentation_dict[pos][tm] = self.segment_fluorescence_accumulation(img, log_parameters, min_area)
        
        with open(self.save_path+'/'+self.experiment+'_'+channel+'_fluorescence_segmentation_dict', 'wb') as handle:
            pickle.dump(segmentation_dict, handle)
    
    
    def get_fluorescence_segmentation_dict(self, channel):
        
        with open(self.save_path+'/'+self.experiment+'_'+channel+'_fluorescence_segmentation_dict', 'rb') as handle:
            segmentation_dict = pickle.load(handle)
            
        return segmentation_dict
    
        
                
                
    
    
    
        
        
    