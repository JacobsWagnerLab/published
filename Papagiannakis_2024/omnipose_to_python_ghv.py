# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 16:25:32 2024

@author: Alexandros Papagiannakis, Christine Jacobs-Wagner lab, HHMI/Stanford University 2024
"""


import os
import scipy as sp
# import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
import time
import Biviriate_medial_axis_estimation as bma
import pickle


# test



class omnipose_to_python_timelapse(object):
    
    def __init__(self, omni_cell_path, experiment, fluorescent_channels, min_trajectory_length, frame_interval, every_nth, save_path):
        
        start = time.time()
        cropped_masks = {}
        cropped_fluor = {}
        bound_box = {}
        cell_frames = {}
        cell_centroid = {}
        divided_index = {}
        daughter_id = {}
        mother_id = {}
        cell_areas = {}
        
        cell_ids_list = os.listdir(omni_cell_path)
        
        for cl in cell_ids_list:
            
            cell_id = int(cl[4:-4])            
            cell_array = sp.io.loadmat(omni_cell_path+'/'+cl)
            birth_frame = cell_array['birth'][0][0]
            death_frame = cell_array['death'][0][0]
            divided = cell_array['divide'][0][0]
            frame_range = np.arange(birth_frame-1, death_frame)
            
            divided_index[cell_id] = divided
            if cell_array['daughterID'].shape[1]>0:
                daughter_id[cell_id] = cell_array['daughterID'][0]
            else:
                daughter_id[cell_id] = np.array([0,0])
            mother_id[cell_id] = cell_array['motherID'][0]
            
            if frame_range.shape[0] > min_trajectory_length:
                
                cell_frames[cell_id] = frame_range
                cropped_masks[cell_id] = {}
                bound_box[cell_id] = {}
                cell_centroid[cell_id] = {}
                cropped_fluor[cell_id] = {}
                cell_areas[cell_id] = []
                
                # print(cell_id, birth_frame, death_frame, divide_frame)
                # print( daughter_id[cell_id], mother_id[cell_id])
                
                for tm in range(len(cell_array['CellA'][0])):
                    # print(tm)
                    if frame_range.shape[0]>1:
                        cropped_masks[cell_id][frame_range[tm]] = sp.ndimage.binary_fill_holes(cell_array['CellA'][0][tm][0][0][3]).astype(int)
                        cell_areas[cell_id].append(np.nonzero(cell_array['CellA'][0][tm][0][0][3])[0].shape[0])
                        cropped_fluor[cell_id][frame_range[tm]] = {}
                        ch_index = 11
                        for ch in fluorescent_channels:
                            cropped_fluor[cell_id][frame_range[tm]][ch]= cell_array['CellA'][0][tm][0][0][ch_index]
                            ch_index+=3
                        cropped_fluor[cell_id][frame_range[tm]]['Phase']= cell_array['CellA'][0][tm][0][0][7]
                        bound_box[cell_id][frame_range[tm]]=cell_array['CellA'][0][tm][0][0][5]
                        cell_centroid[cell_id][frame_range[tm]]=cell_array['CellA'][0][tm][0][0][8][0][0][8][0]

        end = time.time() 
        print(round(end-start,1), 'seconds to load the cropped images and masks')
        
        self.cropped_masks = cropped_masks
        self.cropped_fluor = cropped_fluor
        self.bound_box = bound_box
        self.cell_frames = cell_frames
        self.cell_centroid = cell_centroid
        self.divided_index = divided_index
        self.daughter_id = daughter_id
        self.mother_id = mother_id
        self.cell_areas = cell_areas
        self.save_path = save_path
        self.experiment = experiment
        self.fluorescent_channels = fluorescent_channels+['Phase']
    
    
    def get_lineages_without_plasmid(self, thresholds=(100,300,300)):
        
        low_cell_list = []
        
        thresholded_range = list(np.arange(thresholds[0], thresholds[1]+1))
        
        for cl in self.cell_frames:
            if np.intersect1d(thresholded_range, self.cell_frames[cl]).shape[0]>0:
                bfp_trace = []
                for tm in self.cropped_fluor[cl]:
                    mask = self.cropped_masks[cl][tm]
                    bfp = self.cropped_fluor[cl][tm]['bfp']
                    mean_bfp = bfp[np.nonzero(mask)].mean()
                    bfp_trace.append(mean_bfp)
                bfp_trace = np.array(bfp_trace)[np.nonzero(bfp_trace)]
                if np.mean(bfp_trace)<=thresholds[2]:
                    low_cell_list.append(cl)
        
        for cl in low_cell_list:
            # print(cl)
            if cl in self.mother_id[cl] and self.mother_id[cl] > 0:
                low_cell_list.append(self.mother_id[cl][0])
            if cl in self.daughter_id[cl]:
                for dcl in self.daughter_id[cl]:
                    if dcl > 0:
                        low_cell_list.append(dcl)
                
        return list(np.unique(low_cell_list))
                    
        
    def get_cell_out_of_boundaries(self, limits):
        
        out_of_bound = []
        
        for cl in self.cell_centroid:
            for tm in self.cell_centroid[cl]:
                if int(self.cell_centroid[cl][tm][0]) not in range(limits[0],limits[1]+1) or int(self.cell_centroid[cl][tm][1]) not in range(limits[0],limits[1]+1):
                    out_of_bound.append(cl)
    
        return list(np.unique(out_of_bound))
    
    
    def get_mothers_without_daughters(self):
        
        mothers_list = []
        
        for cl_id in self.mother_id:
            if self.mother_id[cl_id][0]==0:
                mothers_list.append(cl_id)
        
        return mothers_list
    
    
    def get_medial_axes(self, bad_cells, verb=False):
        
        medial_axis_dict = {}
        
        for cl in self.cropped_masks:
            if cl not in bad_cells:
                medial_axis_dict[cl] = {}
                for tm in self.cropped_masks[cl]:
                    if np.mean(self.cropped_fluor[cl][tm]['hu']) != 0:
                        print(cl)
                        medial_axis_dict[cl][tm] = bma.get_medial_axis(self.cropped_masks[cl][tm], radius_px=8, half_angle=22, cap_knot=13, max_degree=60, verbose=verb)
                    else:
                        print('no fluorescence image taken for cell',cl,'at time point',tm)
        with open(self.save_path+'/'+self.experiment+'_medial_axis_dict', 'wb') as handle:
            pickle.dump(medial_axis_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
            
            
    def locate_cell_id(self, cell_position, frame, radius):
        
        centroid_distance = radius
        cell_id_found = 'none'
        
        for cl in self.cell_centroid:
            if frame in self.cell_centroid[cl]:
                centroid_coords = self.cell_centroid[cl][frame]
                distance = np.sqrt((cell_position[0]-centroid_coords[0])**2 + (cell_position[1]-centroid_coords[1])**2)
                if distance < centroid_distance:
                    print(cl)
                    cell_id_found = cl
                    centroid_distance = distance
        
        return cell_id_found
    
    
    def get_lineage_mother(self, single_cell_id):
        
        while self.mother_id[single_cell_id][0]>0:
            single_cell_id = self.mother_id[single_cell_id][0]
            print(single_cell_id)
        
        return single_cell_id
    
    
    def get_oned_fluorescence(self, single_cell_id):
        
        oned_coords_dict = {}
        cell_length_dict = {}
        
        if os.path.isdir(self.save_path+'/'+self.experiment+'_medial_axis_dict'):
            print('Loading the medial axis dataframes...')
            with open(self.save_path+'/'+self.experiment+'_medial_axis_dict', 'rb') as handle:
                medial_axis_dict = pickle.load(handle)
        else:
            print('Drawing the medial axes...')
            medial_axis_dict = {}
            for tm in self.cropped_masks[single_cell_id]:
                if np.mean(self.cropped_fluor[single_cell_id][tm]['hu']) != 0:
                    cell_mask = self.cropped_masks[single_cell_id][tm]
                    
                    medial_axis_dict[tm] = bma.get_medial_axis(cell_mask, 
                                                               radius_px=8, half_angle=22, cap_knot=13, 
                                                               max_degree=60, verbose=True)
                    
                    medial_axis_df = medial_axis_dict[tm][0]
                    
                    cell_mask_df = bma.get_oned_coordinates(cell_mask, medial_axis_df, half_window=5)
                    print('Applying 1D fluorescence...')
                    for ch in self.fluorescent_channels:
                        crop_signal_image = self.cropped_fluor[single_cell_id][tm][ch]
                        cell_mask_df[ch+'_fluor'] = crop_signal_image[np.nonzero(cell_mask)]
                    oned_coords_dict[tm] = cell_mask_df
                    cell_length_dict[tm] = medial_axis_df.arch_length_centered.max()*2
        
        return oned_coords_dict, cell_length_dict
    
    

        



