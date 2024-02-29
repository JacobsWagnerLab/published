#Python 3.7.0 (default, Jun 28 2018, 07:39:16)
#Type "copyright", "credits" or "license" for more information.
#
#IPython 7.8.0 -- An enhanced Interactive Python.

"""
Created on Tue Jan  30 14:32:50 2024

@author: Dr Alexandros Papagiannakis, Christine Jacobs-Wagner lab, Stanford University 2023
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
import warnings
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import skimage
from skimage import io, filters
from skimage.registration import phase_cross_correlation
from skimage.measure import label, regionprops
import scipy
from skimage.filters import threshold_otsu, threshold_local
from scipy import ndimage
from shapely.geometry import Point,Polygon,LineString
from PIL import Image
import pickle
from skimage.morphology import medial_axis, remove_small_objects
import os
from pims import ND2_Reader
from itertools import combinations


class unet_snapshots(object):
    """
    Developer: Alexandros Papagiannakis, Christine Jacobs-Wagner lab, Stanford University, 2023
    
    This class contains all the functions used for particle tracking, particle MSD and speed calculation, 
    correlation between macromolecular and the cell cycle, the ribosomal concentration or the N/C ratio.
    The disclacement of the tracked particles is estimated in the 2D cell projection, as well as relative
    to the centroid or across the medial axis (1D particle projection).
    
    Additional fucntions construct the medial axis of the single cells and the level of constriction.
    
    The code allows for making movies of the tracked particles. 
    
    Functions included in this class:
        __init__
            nd2_to_array
            unet_to_python
        show_unet_masks
        estimate_phase_image_drift
        get_medial_axis
            correct_angle_difference
            get_next_position
            recursive_medial_axis
        run_medial_axis
        check_cell_segmentation
        cell_free_bkg_estimation
        back_sub
    """
    def __init__(self, unet_path, snapshots_path, experiment, save_path):
        """
        The class is iniialized.
        
        Parameters
        ----------
        ___File paths___
        snapshots_paths: the path of the phase and signal images (e.g. HU-mCherry) snapshots. The path to the nd2 file
            A list of paths should be included with the right order (e.g. [phase_path, phase_after_path, signal_path] or [phase_path, signal_path])
            In the case of an empty lists, no snapshots are loaded
        unet_path: the path to the unet .tif cell labels
            If a non-valud file name is used no cell masks are loaded
            For instance, use 'none' to avoid loading cell masks
        save_path: path where the results are saved
        
        ___General parameters___
        experiment: a string describing the experiment (e.g. '07202919_seqA')
        The _init_ function contains a number of sub-functions which use the input parameters of the class
        
        Exceptions
        ----------
        raise ValueError if the XY positions do not include the same number of channels
        """
                   
        def nd2_to_array(images_path):
            """
            This function is used to convert .nd2 images to numpy arrays.
            It also returms the image metadata.
            
            Parameters
            ----------
            image_path - string: the path of the .nd2 file

            Returns
            -------
            [0] the iteration axis - string ('c', 't', 'mc' or 'mct')
            [1] the .nd2 metadata and images (from the ND2_Reader pims module). This is a class object and has multuiple functions
            [2] a dictionary which contains the images as numpy arrays organized by:
                iteration_axis 't' - For fast time lapse (stream acquisition), key2: frame
                iteration_axis 'c' - For snapshots of a single XY position, key1: channel
                iteration_axis 'mc' - For snapshots of multiple XY positions, key1: position, key2: channel
                iteration_axis = 'mct' - For time lapse across different channels, key1: position, key2: channel, key3: time-point
            [3] channels: list of strings - each string represents the channel (e.g. ['Phase', 'mCherry', 'GFP', 'Phase_after'])
                If a certain wavelength (lambda) is used two times in the ND acquisition, then the second channel instance is referred to as '_after'
                An empty list is returned if no channels are selected.
            [4] the number of time-points - positive integer or zero if the Time label was not selected in the ND acquisition
            [5] The number of XY positions - positive integer or zero if the XY label was not selected in the ND acquisition
            
            Notes
            -----
            This function was adapted to include all possible channel, time-point, xy-position permutations in our image acquisition protocols in NIS elements (including the JOBS module)
            New permutations may need to be included for new image permutations.
            The iteration axis determines how the image dimensions are iterated and stored into dictionaries
            """
            # The path of the .nd2 file 
            images = ND2_Reader(images_path)
            # "C:\Users\Alex\Anaconda3\Lib\site-packages\pims_nd2\nd2reader.py"
            # This path has been modified in lines 228 and 229 to accommodate the function.
            print('metadata:',images.metadata)
            print('dimensions:',images.sizes)
            
            scale = round(images.metadata['calibration_um'],3)  # μm/px scale
            sensor = (images.sizes['x'], images.sizes['y'])
            channels = []
            if 'c' in images.sizes:
                # get the channels and frames from the .nd2 metadata
                number_of_channels = images.sizes['c']
                
                for i in range(number_of_channels):
                    ch = images.metadata['plane_'+str(i)]['name']
                    if ch in channels:
                        channels.append(ch+'_after')
                    else:
                        channels.append(ch)   
            # number_of_frames = images.metadata['sequence_count']
            iteration_axis = ''
            if 'm' in images.sizes and images.sizes['m'] > 1:
                iteration_axis += 'm'
                number_of_positions = images.sizes['m']
            if 'c' in images.sizes and images.sizes['c'] > 1:
                iteration_axis += 'c'
            if 't' in images.sizes and images.sizes['t'] > 1:
                iteration_axis += 't'
                number_of_timepoints = images.sizes['t']
            # For a stream acquisition
            if iteration_axis == 't':
                image_arrays = {}
                number_of_positions = 0
                with images as frames:
                    t = 0 # time point
                    print(frames)
                    frames.iter_axes = iteration_axis
                    for frame in frames:
                        image_arrays[t] = np.array(frame)
                        t += 1
                frames.close()
            # For snapshots at different channels
            elif iteration_axis == 'c':
                image_arrays = {}
                number_of_timepoints = 0
                number_of_positions = 0
                with images as frames:
                    i = 0
                    print(frames)
                    frames.iter_axes = iteration_axis
                    for frame in frames:
                        image_arrays[channels[i]] = np.array(frame)
                        i += 1
                frames.close()
            # For snapshots at different XY positions for a single channel (this is how JOBS extracts the snapshots)
            elif iteration_axis == 'm':      
                image_arrays = {}
                number_of_timepoints = 0
                number_of_channels = 1
                with images as frames:
                    i = 0
                    print(frames)
                    frames.iter_axes = iteration_axis
                    for frame in frames:
                        image_arrays[i] = np.array(frame)
                        i += 1
                frames.close()
            # For snapshots at different channels and XY positions
            elif iteration_axis == 'mc':
                image_arrays = {}
                number_of_timepoints = 0
                with images as frames:
                    print(frames)
                    frames.iter_axes = iteration_axis
                    pos = 0
                    ch = 0
                    image_arrays[pos] = {}
                    for frame in frames:
                        if ch < number_of_channels:
                            if pos < number_of_positions:
                                image_arrays[pos][channels[ch]] = np.array(frame)
                                ch+=1
                        elif ch == number_of_channels:
                            pos += 1
                            image_arrays[pos] = {}
                            ch = 0
                            image_arrays[pos][channels[ch]] = np.array(frame)
                            ch+=1
                frames.close()
            # For snapshots at different channels and XY positions and timepoints
            elif iteration_axis == 'mt':
                image_arrays = {}
                with images as frames:
                    print(frames)
                    frames.iter_axes = iteration_axis
                    pos = 0
                    tm = 0
                    image_arrays[pos] = {}
                    for frame in frames:
                        if tm < number_of_timepoints:
                            image_arrays[pos][tm] = np.array(frame)
                            tm+=1
                        elif tm == number_of_timepoints:
                            tm = 0
                            if pos < number_of_positions-1:
                                pos += 1
                                image_arrays[pos] = {}
                                image_arrays[pos][tm] = np.array(frame)
                                tm+=1             
                frames.close()
            # For snapshots at different channels and XY positions and timepoints
            elif iteration_axis == 'mct':
                image_arrays = {}
                with images as frames:
                    print(frames)
                    frames.iter_axes = iteration_axis
                    pos = 0
                    ch = 0
                    tm = 0
                    image_arrays[pos] = {}
                    image_arrays[pos][channels[ch]] = {}
                    for frame in frames:
                        if tm < number_of_timepoints:
                            image_arrays[pos][channels[ch]][tm] = np.array(frame)
                            tm+=1
                        elif tm == number_of_timepoints:
                            tm = 0
                            if ch < number_of_channels-1:
                                ch += 1
                                image_arrays[pos][channels[ch]] = {}
                                image_arrays[pos][channels[ch]][tm] = np.array(frame)
                                tm+=1
                            elif ch == number_of_channels-1:
                                ch = 0
                                pos+=1
                                image_arrays[pos] = {}
                                image_arrays[pos][channels[ch]] = {}
                                image_arrays[pos][channels[ch]][tm] = np.array(frame)
                                tm+=1
                frames.close()
            # if no channels or time points are specified there should be only one image
            elif iteration_axis == '':
                number_of_timepoints = 0
                number_of_positions = 0
                with images as frames:
                    for frame in frames:
                        image_arrays = np.array(frame)
            
            return iteration_axis, images, image_arrays, channels, number_of_timepoints, number_of_positions, scale, sensor


        def unet_to_python(unet_path, experiment, position, pad=5):
            """
            This function incoorporates the masks from Unet to python
            
            Parameters
            ----------
            unet_path: string - the path of the tif image of the cell masks returned by the Unet
            experiment: the experiment ID string
            position: integer corresponding to the XY position (starting at zero)
            pad: integer - the size of the frame around the cells used to crop the cell masks
            
            Returns
            -------
            [0] cropped_cell_masks: a dictionary which inlcudes the cropped cell masks. A cropping pad of 3 pixels is used
            [1] pads: a tuple of 4 coordinates corresponding to the croping rectangle around the cell mask -> (miny, maxy, minx, maxx)
            """
            # check for bad cells in the results destination
            cell_count = 0
            
            if position < 9:
                position_string = 'xy0'+str(position+1)
            elif position >= 9:
                position_string = 'xy'+str(position+1)

            mask_array = io.imread(unet_path)
            cropped_cell_masks = {}
            pads = {}
            for cell in range(1, mask_array.max()+1):
                cell_id = experiment+'_'+position_string+'_'+str(cell)
                cell_mask = np.zeros((mask_array.shape[0], mask_array.shape[1]))
                cell_mask[mask_array==cell]=1
                cell_mask = ndimage.morphology.binary_fill_holes(cell_mask).astype(int)
                y_mask_coords, x_mask_coords = np.where(cell_mask==1)
                minx,miny,maxx,maxy = x_mask_coords.min(), y_mask_coords.min(), x_mask_coords.max(), y_mask_coords.max()
                # remove the masks at the edge of the sensor
                if maxy < mask_array.shape[0]-pad and maxx < mask_array.shape[1]-pad and miny > pad and minx > pad:
                    cropped_cell_masks[cell_id] = cell_mask[(miny-pad):(maxy+pad), (minx-pad):(maxx+pad)]
#                    plt.imshow(cropped_cell_masks[cell_id])
#                    plt.show()
                    pads[cell_id] = (minx-pad, miny-pad, maxx+pad, maxy+pad)
                    cell_count +=1
            print(cell_count, 'cells.')
            return cropped_cell_masks, pads
        
       
        image_arrays = nd2_to_array(snapshots_path)
        channels = image_arrays[3]
        
        if experiment+'_cropped_masks' not in os.listdir(save_path):
            mask_arrays = {}
            mask_dir = os.listdir(unet_path)
            mask_dir = [s for s in mask_dir if '_mask' in s]
            
            for pos in image_arrays[2]:
                pos_mask_str = [s for s in mask_dir if 'xy'+str(pos+1)+'c1_mask' in s][0]
                mask_arrays[pos] = unet_to_python(unet_path+'/'+pos_mask_str, experiment, pos, pad=5)
                    
            with open(save_path+'/'+experiment+'_cropped_masks', 'wb') as handle:
                pickle.dump(mask_arrays, handle)
        else:
            print('found the cropped masks...loading')
            with open(save_path+'/'+experiment+'_cropped_masks', 'rb') as handle:
                mask_arrays = pickle.load(handle)
     

        self.image_arrays = image_arrays
        self.mask_arrays = mask_arrays
        self.experiment = experiment
        self.save_path = save_path
        self.channels = channels
        self.sensor =  self.image_arrays[7] 
        self.unet_path = unet_path
        
        def count_cells(mask_arrays):
            cell_count = 0
            for pos in mask_arrays:
                cell_count+=len(mask_arrays[pos][0])
            return cell_count
    
        cell_count = count_cells(self.mask_arrays)
        self.cell_count = cell_count
        print(self.cell_count, 'total segmented cells...')
            

    # CHECK THE UNET SEGMENTATION 
    def show_unet_masks(self, save=False, curate=False):
    # def show_oufti_meshes(self, channel_offsets):
        """
        This function can be used to plot the perimeter of the Unet masks after processing.
        The perimeter of the masks is plotted over the phase image.
        """
        
        try:
            print('found bad cell coordinates. Getting bad cells...')
            bad_cells, bad_cell_coordinates = self.get_bad_cells()
        except:
            print('bad cells have not been curated. Generating data structures...')
            bad_cells = []
            bad_cell_coordinates = {}
            for pos in self.mask_arrays:
                bad_cell_coordinates[pos] = []
                
        for pos in self.mask_arrays:
            print('printing the cell boundaries on the phase contrast image...')
            plt.figure(figsize=(15,15))
            plt.imshow(self.image_arrays[2][pos][self.channels[0]])
            for test_cell in self.mask_arrays[pos][0]:
                if test_cell not in bad_cells:
                    cropped_mask = self.mask_arrays[pos][0][test_cell]
                    dilated_cell_mask = scipy.ndimage.morphology.binary_dilation(cropped_mask, iterations=1)
                    cell_mesh_y, cell_mesh_x = zip(*skimage.measure.find_contours(dilated_cell_mask, level=0.5)[0])
                    cell_pad = self.mask_arrays[pos][1][test_cell]
                    plt.plot(cell_mesh_x+cell_pad[0], cell_mesh_y+cell_pad[1], color='white')
            if save==True:
                plt.savefig(self.save_path+'/segmentation_'+self.experiment+'_xy'+str(pos+1)+'.jpeg')
            if curate==True:
                bad_coords = plt.ginput(n=0, timeout=0)
                print(bad_coords)
                if bad_cell_coordinates[pos] == 'none':
                    bad_cell_coordinates[pos] = []
                bad_cell_coordinates[pos] += bad_coords
                print('bad cells in locations:',bad_cell_coordinates[pos])
            plt.show()
            plt.close()
            
        with open(self.save_path+'/'+self.experiment+'_bad_cell_coordinates', 'wb') as handle:
            pickle.dump(bad_cell_coordinates, handle, protocol=pickle.HIGHEST_PROTOCOL)
            
            
    def get_unet_mask(self, position):
        
        position_string = 'xy'+str(position+1)+'c1'
        
        unet_dir = os.listdir(self.unet_path)
        cor_pos = position*len(self.channels)
        cor_pos = position_string
        image_string =  [s for s in unet_dir if cor_pos in s][0]
        single_unet_path = self.unet_path+'/'+image_string[0:-4]+'_mask.tif'
        mask_labels = io.imread(single_unet_path)
        return mask_labels
    
    
    def get_cell_id_from_label(self, mask_label, position_string):
        
        return self.experiment+'_'+position_string+'_'+str(mask_label)
    
    
    def show_selected_cells(self, position, channel, crop_square):
        
        position_string = 'xy'+str(position+1)+'c1'
        
        try:
            phase_image = self.image_arrays[2][position]['Trans']
        except KeyError:
            phase_image = self.image_arrays[2][position]['Phase']
        
        signal_image = self.image_arrays[2][position][channel]
        bkg_cor_image = self.back_sub(phase_image, signal_image, 25, 128, 60, False)[2]
        
        fluor_image = bkg_cor_image[crop_square[2]:crop_square[3],
                                    crop_square[0]:crop_square[1]] 
        
        mask_dir = os.listdir(self.unet_path)
        mask_dir = [s for s in mask_dir if '_mask' in s]
        pos_mask_str = [s for s in mask_dir if 'xy'+str(position+1)+'c1_mask' in s][0]
                        
        cell_masks = io.imread(self.unet_path+'/'+pos_mask_str)
        
        cropped_labels = cell_masks[crop_square[2]:crop_square[3],crop_square[0]:crop_square[1]]
        cropped_phase = phase_image[crop_square[2]:crop_square[3],crop_square[0]:crop_square[1]]
        
        cell_labels = []
        
        for lbl in np.unique(cropped_labels):
            if lbl > 0:
                cell_labels.append(self.get_cell_id_from_label(lbl, position_string))
                
        return fluor_image, cropped_phase, cropped_labels, cell_labels
        
            
    
    def get_bad_cells(self, show=False):
        
        bad_cells = []
        
        with open (self.save_path+'/'+self.experiment+'_bad_cell_coordinates', 'rb') as handle:
            bad_cell_coordinates = pickle.load(handle)
        print('getting bad cell IDs:')
        
        for pos in bad_cell_coordinates:
            
            if pos < 9:
                position_string = 'xy0'+str(pos+1)
            elif pos >= 9:
                position_string = 'xy'+str(pos+1)
            
            mask_labels = self.get_unet_mask(pos)  
          
            for bd_crd in bad_cell_coordinates[pos]:
                try:
                    bad_label = mask_labels[int(bd_crd[1]),int(bd_crd[0])]
                except IndexError:
                    bad_label = 0
                if bad_label>0:
                    # bad_cell_id = self.experiment+'_'+position_string+'_'+str(bad_label)
                    bad_cell_id = self.get_cell_id_from_label(bad_label, position_string)
                    print(bad_cell_id)
                    if show==True:
                        if bad_cell_id in self.mask_arrays[pos][0]:
                            plt.imshow(self.mask_arrays[pos][0][bad_cell_id])
                            plt.show()
                    bad_cells.append(bad_cell_id)
        print(bad_cells)
        with open(self.save_path+'/'+self.experiment+'_bad_cells', 'wb') as handle:
            pickle.dump(bad_cells, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        return bad_cells, bad_cell_coordinates
                
                    

     # MEDIAL AXIS ESTIMATION #
    def get_angle_from_slope(self, displacements, slope='none'):
        """
        This function is used to estimate the angle of each microtubule displacement from two adjacent time-points.
        The angle is estmated using the slope of the microtubule displacement and the orientation (negative or positive dx and dy)
        
        Specifically, the slope is converted to an angle.
        Using the dx and dy displacements it is decided to which of the quartiles in a circle the angle belongs and the appropriate adjustments are made.
        
        Parameters
        ----------
        displacements: tuple - (dx, dy). The displacement between two linked spots in the x and y direction
        slope: float - the slope of the displacement
               if slope == 'none' the slope is estimated as dy/dx
        
        Returns
        -------
        angle: float - the angle of the displacement in degrees
        
        Notes
        -----
        This function is used for the medial axis extimation: self.get_medial_axis()
        """
        dx = displacements[0]
        dy = displacements[1]
        if dx != 0:
            if slope =='none':
                slope = dy/dx
            angle = np.rad2deg(np.arctan(slope))
            if angle >= 0:
                angle = angle
            elif angle < 0:
                angle  = 360 + angle
            if slope >= 0:
                if dx >= 0:
                    angle = angle
                elif dx < 0:
                    angle = 180 + angle
            elif slope < 0:
                if dx >= 0:
                    angle = angle
                elif dx< 0:
                    angle = angle - 180
        elif dx == 0:
            if dy > 0:
                angle = 90
            elif dy < 0:
                angle = 270
            elif dy == 0:
                angle = 0
        return angle


    def get_medial_axis(self, cropped_cell_mask, cell_pad, radius_px=8, half_angle=22, cap_knot=13, max_degree=60):
        """
        This function construct the medial axis of a signle cell, 
        as well as the relative coordinates of the cell from one pole to the other.
        
        Parameters
        ----------
        the cell ID (date_xyPosition_cellNumber produced after initializing the class)
        radius_px: positive integer - the radius beyond which the algorithm searches for the next anchor point in the circle sector
        half_angle: positive_float - the half angle of the circle sector within which the next anchor point is searched for.
            e.g. an angle of 22 degrees indicates that the next angle will be searched for in the secotr of 44 degrees (quartile), 22 degrees adjuscent the previous orientation
        cap_px: positive_interger - the number of knots excluded at the poles. This region will be extended using the anlge from the previous 10 anchors. 
        max_degree: positive_integer - the maximum degree of the fitted polynomials
        
        Returns
        -------
        [0] A pandas DataFrame including the absolute, scaled and relative coordinates of the medial axis
            Columns:
                'cropped_x': the cropped x coordinates of the medial axis (cropped by the cell pad)
                'cropped_y': the croppedyx coordinates of the medial axis (cropped by the cell pad)
                'arch_length': the arch length of the medial axis along the cell length
                'arch_length_centered': the arch length of the medial axis scaled by the centroid
                'arch_length_scaled': the relative arch length from -1 to 1
        [1] The x,y coordinates of the cell centroid
        [2] The croped x,y coordinates of the cell centroid, cropped by the cell pad
        """

        # GET THE CELL MASK AND THE DISTANCE TRANSFORMATION OF THE CELL MASK
        cell_mask = cropped_cell_mask
        # get the cropped mask of the single cell
#        dilated_cell_mask =  scipy.ndimage.morphology.binary_dilation(cell_mask, iterations=2)
        # get the resized cell mask
        resized_cell_mask = np.array(Image.fromarray(cell_mask).resize((cell_mask.shape[1]*10, cell_mask.shape[0]*10), resample=Image.NEAREST))
        skel, dist = medial_axis(resized_cell_mask, return_distance=True)
#        dist = dist * (dist>5) # to set a distance threshold in the distance transformed mask
    
#        plt.figure(figsize=(10,8))
#        plt.imshow(dist)
#        plt.plot(np.nonzero(skel)[1], np.nonzero(skel)[0], 'o', markersize=0.2)
#        plt.show()
        
        # GET THE FIRST ANCHOR POINT AT THE CENTER OF THE MAX LENGTH DIMENSION
        # CORRESPONDING TO THE POINT WITH THE MAXIMUM DISTANCE FROM THE CELL EDGES
        # THE ANGLE OF THE CELL AT THE FIRST ANCHOR POINT IS ALSO ESTIMATED
        len_y = len(dist)
        len_x = len(dist[0])
        length = (len_x, len_y)
        # the index of the longest coordinate (x or y)
        max_index = length.index(max(length))
        if max_index == 0:
            half_x = int(round(len_x/2 ,0))
            half_y = np.argmax(np.transpose(dist)[half_x])
        elif max_index == 1:
            half_y = int(round(len_y/2, 0))
            half_x = np.argmax(dist[half_y])
            
        start_x = half_x
        start_y = half_y
        cropped_window = dist[(start_y-10):(start_y+11), (start_x-10):(start_x+11)]
#        plt.imshow(cropped_window)
#        plt.show()
        window_df = pd.DataFrame()
        window_df['x'] = np.nonzero(cropped_window)[1]
        window_df['y'] = np.nonzero(cropped_window)[0]
        window_df['fluor'] = cropped_window[np.nonzero(cropped_window)].ravel()
        window_df['distance'] = np.sqrt((window_df.x-10)**2 + (window_df.y-10)**2)
        window_df = window_df[window_df.distance>5]
        window_df = window_df[window_df.fluor == window_df.fluor.max()]
        if window_df.shape[0] == 1:
            start_angle = self.get_angle_from_slope((window_df.x.values[0]-10, window_df.y.values[0]-10))
        elif window_df.shape[0] > 1:
            window_df = window_df[window_df.distance == window_df.distance.max()]
            start_angle = self.get_angle_from_slope((window_df.x.values[0]-10, window_df.y.values[0]-10))
#        print(start_angle)
        
        # THIS CODE CORRECTS THE ANGLE DIFFERENCE
        def correct_angle_difference(source_angle, destination_angle):
            """
            This function is used to correct the difference between two angles.
            It returns a positive angle smaller than 180 degrees.
            """
            a = destination_angle - source_angle
            if a >= 180:
                return 360-a
            elif a <= -180:
                return 360+a
            else:
                return abs(a)
         
        # THIS CODE ESTIMATES THE NEXT ANCHOR POINT USING THE PREVIOUS ONE AND THE ANGLE
        def get_next_position(dist, x, y, angle, list_of_knots):
            """
            This function scans the cell mask distance transformation for the max distance knots
            that will be used to fit the medial axis.
            
            The knots have a resolution of 0.1 pixels.
            """
            dist_temp = dist.copy()
            dist_temp[y][x]=0
            # radius_px = 8 # good parameter
            crop_dist = dist_temp[(y-radius_px):(y+radius_px+1), (x-radius_px):(x+radius_px+1)]
#            col = np.argmax(np.amax(crop_dist, axis=1))
#            row = np.argmax(np.amax(crop_dist, axis=0))
            y_coords, x_coords = np.nonzero(crop_dist)
            intensities = crop_dist[np.nonzero(crop_dist)]
            
            intensity_df = pd.DataFrame()
            intensity_df['x'] = x_coords + x -radius_px
            intensity_df['y'] = y_coords + y -radius_px
            intensity_df['fluor'] = intensities
            
            intensity_df['dx'] =  intensity_df['x']-x
            intensity_df['dy'] =  intensity_df['y']-y
            intensity_df['distance'] = np.sqrt(intensity_df.dx**2 + intensity_df.dy**2)
            intensity_df['angle'] = intensity_df.apply(lambda row: self.get_angle_from_slope((row.dx, row.dy)), axis=1)
            intensity_df['angle_dif'] = intensity_df.apply(lambda row: correct_angle_difference(angle, row.angle), axis=1)
            
#            intensity_df = intensity_df[intensity_df.angle_dif<=45] 
            intensity_df = intensity_df[(intensity_df.angle_dif<=half_angle) & (intensity_df.distance>6)] 
            
            if intensity_df.shape[0]>0:
                max_df = intensity_df[intensity_df.fluor == intensity_df.fluor.max()] #new
#                max_df = intensity_df[intensity_df.angle_dif == intensity_df.angle_dif.min()]
                if max_df.shape[0] > 0:
                    if max_df.shape[0] > 1:
    #                    max_df['distance'] = np.sqrt(max_df.dx**2 + max_df.dy**2)
    #                    max_df = max_df[max_df.distance==max_df.distance.min()]
                        max_df = max_df[max_df.angle_dif==max_df.angle_dif.min()] #new
#                        max_df = max_df[max_df.fluor==max_df.fluor.max()] # old
    #                max_index = max_df.index[0]
                    new_x = max_df.x.values[0]
                    new_y = max_df.y.values[0]
                    new_angle = max_df.angle.values[0]
                    max_fluor = max_df.fluor.values[0]
    #                print(max_fluor)
                
                    if (new_x, new_y) not in list_of_knots:
                        if max_fluor >= 3:
    #                        print(new_x,  new_y, new_angle)
                            return new_x,  new_y, new_angle
                        elif max_fluor < 3:
                            print('This is the end of the cell:', index_increment)
                            return False
                    elif (new_x, new_y) in list_of_knots:
                        print('This is the end of the cell and a loop is formed:', index_increment)
                        return 'loop'
                else:
                    print('This is the end of the cell:', index_increment)
                    return False
            elif intensity_df.shape[0]==0:
                print('This is the end of the cell:', index_increment)
                return False
        
        
        # RECURSIVE ALGORITHM TO GET ALL ANCHOR POINTS
        def recursive_medial_axis(index_increment, dist, x, y, angle, index, list_of_knots):
            """
            This is a function that runs the "get_next_position" finction recursively.
            """
            new_knot =get_next_position(dist, x, y, angle, list_of_knots)
            if new_knot != False:
                if new_knot != 'loop':
                    new_x,new_y,new_angle = new_knot
                    list_of_knots.append((new_x, new_y))
    #                print(new_x, new_y, new_angle)
                    if index_increment == 1:
                        index += 1
                    if index_increment == -1:
                        index -= 1
                    xyz_coord_list.append((new_x, new_y, index))
                    
                    x_list, y_list, z_list = list(zip(*xyz_coord_list))
                    pre_df = pd.DataFrame()
                    pre_df['x'] = np.array(x_list)/10
                    pre_df['y'] = np.array(y_list)/10
                    pre_df['z'] = z_list
                    pre_df = pre_df.sort_values(['z'])
        
                    line_coords = list(map(lambda x, y:(x,y), pre_df.x, pre_df.y))
                    line = LineString(line_coords)
    #                input()
                    if line.is_simple == True:
                        recursive_medial_axis(index_increment, dist, new_x, new_y, new_angle, index, list_of_knots)
                    elif line.is_simple == False:
                        print('This is the end of the cell and a loop is formed...')
                # remove the loops
                elif new_knot == 'loop':
                    for i in range(20):
                         xyz_coord_list.pop()
                
        # Run the recursive algorithms to get the anchor points for the central line fit
        index_increment = 1 # run towards one dimension
        index = 0
        list_of_knots = [(start_x, start_y)]
        xyz_coord_list = [(start_x, start_y, index)]
        x = start_x
        y = start_y
        angle = start_angle
        recursive_medial_axis(index_increment, dist, x, y, angle, index, list_of_knots)
        xyz_coord_list_1 = xyz_coord_list.copy()
        index_increment = -1 # run towards the opposite dimension
        index = 0
        list_of_knots = [(start_x, start_y)]
        xyz_coord_list = [(start_x, start_y, index)]
        x = start_x
        y = start_y
        angle = start_angle + 180
        if angle >= 360:
            angle = angle-360
        recursive_medial_axis(index_increment, dist, x, y, angle, index, list_of_knots)
        
        xyz_coord_list = xyz_coord_list_1 + xyz_coord_list[1:] # combine the two lists of coordinates
        
        # GETTING THE XY COORDINATES OF THE ANCHOR POINTS
        # getting the x,y and z coordinates of the knots
        x_list, y_list, z_list = list(zip(*xyz_coord_list))
#        plt.figure(figsize=(10,10))
#        plt.imshow(resized_cell_mask)
#        plt.plot(x_list, y_list, 'o', markersize=1)
#        plt.show()
#        plt.figure(figsize=(10,10))
#        plt.scatter(x_list, y_list, c=z_list)
#        plt.show()
        # CHECKING IF THE CENTRAL LINE INTERSECTS ITSELF AND REMOVING THE LOOPS
        # rescaling and sorting the coordinates of the knots
        pre_df = pd.DataFrame()
        pre_df['x'] = np.array(x_list)/10
        pre_df['y'] = np.array(y_list)/10
        pre_df['z'] = z_list
        pre_df = pre_df.sort_values(['z'])
        line_coords = list(map(lambda x, y:(x,y), pre_df.x, pre_df.y))
        line_coords = line_coords[::2]
        line = LineString(line_coords)
        positive_intersection = False
        negative_intersection = False
        if line.is_simple == False:
            print('removing the loop...')
            list_of_lines = []
            line_intersections = []
            for pnt in range(len(line_coords)-1):
                list_of_lines.append(LineString(line_coords[pnt:pnt+2]))
            for l1, l2 in combinations(list_of_lines,2): #For all combinations of segments
#                print(l1.intersection(l2).coords[:])
                if l1.crosses(l2) == True: #Find crossings
#                    print('cross')
                    line_intersections.append(l1.intersection(l2).coords[:][0])
            intersection_points_positive = []
            intersection_points_negative = []
            for intersect in line_intersections:
                distance_df = pre_df.copy()
                distance_df['inter_distance'] = np.sqrt((distance_df.x - intersect[0])**2+(distance_df.y - intersect[1])**2)
                distance_df = distance_df.sort_values(['inter_distance'])
                distance_df = distance_df[0:4]
                intersection_point = distance_df[distance_df.z.abs()==distance_df.z.abs().min()].z.values[0]
                if intersection_point > 0:
                    intersection_points_positive.append(intersection_point)
                elif intersection_point < 0:
                    intersection_points_negative.append(intersection_point)
            if len(intersection_points_positive) > 0:
                pre_df = pre_df[pre_df.z < intersection_points_positive[np.argmax(intersection_points_positive)]]
                positive_intersection = True
            elif len(intersection_points_negative) > 0:
                pre_df = pre_df[pre_df.z > intersection_points_negative[np.argmin(intersection_points_negative)]]
                negative_intersection = True
#        plt.imshow(resized_cell_mask)
#        plt.plot(pre_df.x*10, pre_df.y*10, 'o')
#        plt.show()
                
        # TRUNCATE THE MEDIAL AXIS COORDINATES FROM THE EDGES
        if pre_df.shape[0]>2*cap_knot+5:
            pre_df_max_index = cap_knot
        elif pre_df.shape[0]<=2*cap_knot+5:
            pre_df_max_index = pre_df.shape[0]/2 - 5

        if positive_intersection == False and negative_intersection == False:
            truncated_df = pre_df[pre_df_max_index:-pre_df_max_index]
        elif positive_intersection == False and negative_intersection == True:
            truncated_df = pre_df[0:-pre_df_max_index]
        elif positive_intersection == True and negative_intersection == False:
            truncated_df = pre_df[pre_df_max_index:]
        elif positive_intersection == True and negative_intersection == True:
            truncated_df = pre_df
#        plt.imshow(cell_mask)
#        plt.plot(truncated_df.x, truncated_df.y, 'o', markersize=0.2)
#        plt.show()
        
        # EXTENDING THE CENTRAL LINE AT THE EDGES USING THE AVERAGE ANGLE FROM THE 10 PREVIOUS ANCHORS
        # For the negative side
        if truncated_df.shape[0] >= 10:
            trunc_index = 10
        elif truncated_df.shape[0] < 10:
            trunc_index = truncated_df.shape[0]
        
        slope_1_df = truncated_df[0:trunc_index]
        slope_1 = np.polyfit(slope_1_df.x, slope_1_df.y, 1)[0]
        x_1 = np.array(slope_1_df.x)
        y_1 = np.array(slope_1_df.y)
        dx_1 = round(x_1[0] - x_1[trunc_index-1], 0)
        dy_1 = round(y_1[0] - y_1[trunc_index-1],0)

        angle_1 = self.get_angle_from_slope((dx_1, dy_1), slope_1)
        
        # For the positive side
        slope_2_df = truncated_df[-trunc_index:]
        slope_2 = np.polyfit(slope_2_df.x, slope_2_df.y, 1)[0]
        x_2 = np.array(slope_2_df.x)
        y_2 = np.array(slope_2_df.y)
        dx_2 = round(x_2[trunc_index-1] - x_2[0],0)
        dy_2 = round(y_2[trunc_index-1] - y_2[0],0)
        angle_2 = self.get_angle_from_slope((dx_2, dy_2), slope_2)
        
        min_z = truncated_df.z.min()
        max_z = truncated_df.z.max()
        
        x_list_1 = [x_1[0]]
        y_list_1 = [y_1[0]]
        z_list_1 = [min_z]
        
        # extend towards the negative side using the average angle at the edge of the central line
        for i in range(55):
            x_list_1.append(x_list_1[-1]+0.5*np.cos(angle_1*np.pi/180))
            y_list_1.append(y_list_1[-1]+0.5*np.sin(angle_1*np.pi/180))
            z_list_1.append(z_list_1[-1]-1)
        
        x_list_2 = [x_2[-1]]
        y_list_2 = [y_2[-1]]
        z_list_2 = [max_z]
        
        # extend towards the positive side using the average angle at the edge of the central line
        for i in range(55):
            x_list_2.append(x_list_2[-1]+0.5*np.cos(angle_2*np.pi/180))
            y_list_2.append(y_list_2[-1]+0.5*np.sin(angle_2*np.pi/180))
            z_list_2.append(z_list_2[-1]+1)
            
        x_list_final = x_list_1[1:]+x_list_2[1:]
        y_list_final = y_list_1[1:]+y_list_2[1:]
        z_list_final = z_list_1[1:]+z_list_2[1:]
        
        pre_df_2 = pd.DataFrame()
        pre_df_2['x'] = x_list_final
        pre_df_2['y'] = y_list_final
        pre_df_2['z'] = z_list_final
        pre_df = pd.concat([truncated_df, pre_df_2])
        pre_df = pre_df.sort_values(['z'])
#        plt.imshow(cell_mask)
#        plt.plot(pre_df.x, pre_df.y)
#        pre_df = pre_df[50:-50]
        # a bivariate spline is fitted to the data
#        tck, u = interpolate.splprep([pre_df.x, pre_df.y, pre_df.z], k=1, s=50)
##        tck, u = interpolate.splprep([pre_df.x, pre_df.y, pre_df.z], k=1, s=30)
#        x_knots, y_knots, z_knots = interpolate.splev(tck[0], tck)
##        u_fine = np.linspace(-1,2,pre_df.shape[0]*100+1000)
#        u_fine = np.linspace(-2,3,pre_df.shape[0]*300)
#        x_hat, y_hat, z_hat = interpolate.splev(u_fine, tck)
        
        # FIT A nth DEGREE POLYNOMIAL TO THE EXTENDED CENTRAL LINES
        # use this code if a polynomial fit is preferred versus a bivariate spline
        polynomial_degree = int(pre_df.shape[0]/10-5)
        if polynomial_degree > max_degree:
            polynomial_degree = max_degree
        print('fitting polynomial fucntions of degree:', polynomial_degree)
        extended_z = np.arange(pre_df.z.min(), pre_df.z.max(), 0.01)
        fit_x = np.polyfit(pre_df.z, pre_df.x, polynomial_degree) #nth degree polynomial fit long the X axis
        x_hat = np.polyval(fit_x, extended_z)
#        plt.plot(pre_df.z, pre_df.x, 'o')
#        plt.plot(extended_z, x_hat)
#        plt.show()
        fit_y = np.polyfit(pre_df.z, pre_df.y, polynomial_degree) #nth degree polynomail fit along the Y axis
        y_hat = np.polyval(fit_y, extended_z)
#        plt.plot(pre_df.z, pre_df.y, 'o')
#        plt.plot(extended_z, y_hat)
#        plt.show()
        # REMOVE THE CENTRAL LINE COORDINATES THAT DO NOT FALL INTO THE CELL MASK
        # getting only those coordinates of the medial axis that fall into the original cell mask
        x_fine_round = np.around(x_hat,0).astype(int)
        y_fine_round = np.around(y_hat,0).astype(int)
        good_indexes = (x_fine_round<cell_mask.shape[1])*(y_fine_round<cell_mask.shape[0])
        good_indexes_2 = (x_fine_round>=0)*(y_fine_round>=0)
        good_indexes= good_indexes * good_indexes_2
        x_fine = x_hat[good_indexes]
        y_fine = y_hat[good_indexes]
        x_fine_round = x_fine_round[good_indexes]
        y_fine_round = y_fine_round[good_indexes]
        nonzero_medial_indexes = np.nonzero(cell_mask[y_fine_round, x_fine_round])     
        x_fine_good = x_fine[nonzero_medial_indexes]
        y_fine_good = y_fine[nonzero_medial_indexes]
#        plt.imshow(cell_mask)
#        plt.plot(x_fine_good, y_fine_good)
#        plt.show()
        # GENERATE THE RELATIVE CELL COORDINATES AND THE CENTROID
        # generate the medial axis dataframe
        medial_axis_df = pd.DataFrame()
        medial_axis_df['cropped_x'] = x_fine_good
        medial_axis_df['cropped_y'] = y_fine_good
        # get the arch length of the medial axis
        delta_x_sqr = (x_fine_good[1:] - x_fine_good[0:-1])**2
        delta_y_sqr = (y_fine_good[1:] - y_fine_good[0:-1])**2
        disp_array = np.sqrt(delta_x_sqr + delta_y_sqr)
        disp_list = [0]
        for disp in disp_array:
            disp_list.append(disp_list[-1]+disp)
        medial_axis_df['arch_length'] = disp_list 
        medial_axis_df['arch_length_centered'] = disp_list - np.max(disp_list)/2
        medial_axis_df['arch_length_scaled'] = medial_axis_df['arch_length_centered'] / medial_axis_df['arch_length_centered'].max()
        # get the cropped centroid of the medial axis
        center_df = medial_axis_df[medial_axis_df.arch_length_centered.abs()==medial_axis_df.arch_length_centered.abs().min()]
        cropped_centroid = (center_df.cropped_x.mean(), center_df.cropped_y.mean())
        # PLOT THE CELL MASK WITH THE CENTRAL LINE AND THE CENTROID
        # plot and medial axis on the cell mask and the centroid
        plt.imshow(cell_mask)
        plt.plot(medial_axis_df.cropped_x, medial_axis_df.cropped_y, color='red')
        plt.plot(cropped_centroid[0], cropped_centroid[1], 'o')
        plt.show()
        # get the original centroid coordinates from the cropped centroid
        centroid = (cropped_centroid[0]+cell_pad[0], cropped_centroid[1]+cell_pad[1])
        
        return medial_axis_df, centroid, cropped_centroid
    
    
    def run_medial_axis(self):
        """
        This function gets the medial axis for all cells in the image,
        assembles the results in a dictionary, which is then saved in the specified folder.
        
        Exception
        ---------
        If the medial axis cannot be estimated, the cell is not included in the medial axis dictionary.
        """
        medial_axis_dict = {}
        for pos in self.mask_arrays:
            medial_axis_dict[pos]={}
            for cell in self.mask_arrays[pos][0]:
                cropped_cell_mask = self.mask_arrays[pos][0][cell]
                cell_pad = self.mask_arrays[pos][1][cell]
                try:
                    medial_axis_dict[pos][cell] = self.get_medial_axis(cropped_cell_mask, cell_pad)
                except:
                    plt.imshow(cropped_cell_mask)
                    plt.show()
                    print('This spec is probably not a cell')
            
        with open(self.save_path+'/'+self.experiment+'_medial_axis_dict', 'wb') as handle:
            pickle.dump(medial_axis_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    

    # BACKGROUND ESTIMATION AND SUBTRACTION FUNCTIONS 
    def cell_free_bkg_estimation(self, masked_signal_image, step):
        """
        This function scans the image using squared regions of specified size (step) 
        and applies the average cell-free background fluorescence per region.
        This function is used in the self.back_sub() function.
        
        Parameters
        ----------
        masked_signal_image: 2D numpy array - the signal image were the cell pixels are annotated as 0 
                             and the non-cell pixels maintain their original grayscale values
        step: integer (should be a divisor or the square image dimensions) - the dimensions of the squared region where 
              the cell-free background fluorescence is averaged
                example: for an 2048x2048 image, 128 is a divisor and can be used as the size of the edge of the square 
    
        Returns
        -------
        A 2D numpy array where the cell-free average background is stored for each square region with specified step-size
        """
        zero_image = np.zeros(self.sensor) # initiated an empty image to store the average cell-free background
        
        for y in range(0, self.sensor[1], step):
            for x in range(0, self.sensor[0], step):
                # cropped_image = img_bkg_sig[y:(y+step), x:(x+step)]
                cropped_mask = masked_signal_image[y:(y+step), x:(x+step)]
#                mean_bkg = np.mean(cropped_mask[np.nonzero(cropped_mask)].ravel()) # get the mean of the non-zero pixels
#                mean_bkg = scipy.stats.mode(cropped_mask[cropped_mask!=0].ravel())[0][0] # get the mode of the non-zero pixels
                mean_bkg = np.median(cropped_mask[np.nonzero(cropped_mask)].ravel()) # get the mean of the non-zero pixels
                zero_image[y:(y+step), x:(x+step)] = mean_bkg # apply this mean fluorescence to the original empty image
                       
        return zero_image
    
    
    def back_sub(self, phase_image, signal_image, dilation, estimation_step, smoothing_sigma, show):
        """
        Subtracts an n_order second degree polynomial fitted to the non-cell pixels.
        The 2D polynomial surface is fitted to the non-cell pixels only.
            The order of the polynomial depends on whether there is uneven illumination or not
        The non-cell pixels are masked as thos below the otsu threshold estimated on the basis of the inverted phase image.
        
        Parameters
        ----------
        signal_image: numpy.array - the image to be corrected
        dilation: non-negative integer - the number of dilation rounds for the cell mask
        estimation_step: positive_integer - the size of the square edge used for average background estimation
        smoothing_sigma: non-negative integer - the smoothing factor of the cell free background
        show: binary - True if the user wants to visualize the 2D surface fit
        
        Returns
        -------
        [0] The average background
        [1] The background corrected image (after subtracting the 2D polynomial surface)
        """
        print('Subtracting background...')
        # invert the image and apply an otsu threshold to separate the dimmest 
        # (or inversely brightest pixels) which correspond to the cells
        inverted_phase_image = 1/phase_image
        inverted_threshold = threshold_otsu(inverted_phase_image.ravel())
        phase_mask = inverted_phase_image > inverted_threshold
        # dilate the masked phase images
        threshold_masks_dil = ndimage.binary_dilation(phase_mask, iterations=dilation)
        threshold_masks_dil = np.array(threshold_masks_dil)
        # mask the signal image, excluding the dilated cell pixels
        masked_signal_image = signal_image * ~threshold_masks_dil
        if show == True:
            plt.figure(figsize=(10,10))
            plt.imshow(threshold_masks_dil)
            plt.show()
        # The dimensions of the averaging square
        step = estimation_step
        img_bkg_sig = self.cell_free_bkg_estimation(masked_signal_image, step)
        if show == True:
            plt.figure(figsize=(20,20))
            plt.imshow(img_bkg_sig, cmap='Blues')
            plt.clim(np.mean(img_bkg_sig.ravel())-5*np.std(img_bkg_sig.ravel()),np.mean(img_bkg_sig.ravel())+2.5*np.std(img_bkg_sig.ravel()))
            plt.colorbar()
            plt.show()
        # Smooth the reconstructed background image, with the filled cell pixels.
        img_bkg_sig = img_bkg_sig.astype(np.int16)
        img_bkg_sig = ndimage.gaussian_filter(img_bkg_sig, sigma=smoothing_sigma)
        norm_img_bkg_sig = img_bkg_sig/np.max(img_bkg_sig.ravel())
        if show == True:
            plt.figure(figsize=(20,20))
            plt.imshow(img_bkg_sig, cmap='Blues')
            # plt.clim(0,25*np.std(bkg_cor.ravel()))
            plt.colorbar()
            plt.show()
        # subtract the reconstructed background from the original signal image
        bkg_cor = (signal_image - img_bkg_sig)/norm_img_bkg_sig
        bkg_cor_2 = signal_image - img_bkg_sig
        # use this line if you want to convert negative pixels to zero
        # bkg_cor[bkg_cor<0]=0
        if show == True:
            plt.figure(figsize=(20,20))
            plt.imshow(bkg_cor, cmap='Blues')
            plt.clim(0,25*np.std(bkg_cor.ravel()))
            plt.colorbar()
            plt.show()
            plt.figure(figsize=(20,20))
            plt.imshow(img_bkg_sig*threshold_masks_dil, cmap='Blues')
            plt.clim(np.mean(img_bkg_sig.ravel())-5*np.std(img_bkg_sig.ravel()),np.mean(img_bkg_sig.ravel())+2.5*np.std(img_bkg_sig.ravel()))
            plt.colorbar()
            plt.show()
        
        return bkg_cor, np.mean(img_bkg_sig.ravel()), bkg_cor_2
    

    def get_background_corrected_images(self):
#        plt.imshow(snap.back_sub(phase_image, signal_image, 25, 128, 60, True)[2])
#        plt.colorbar()
        bkg_cor_arrays = {}
        fluor_channels = self.channels.copy()
        try:
            fluor_channels.remove('Trans')
        except ValueError:
            fluor_channels.remove('Phase')
        for ch in fluor_channels:
            bkg_cor_arrays[ch]={}
            for pos in self.image_arrays[2]:
                print(self.experiment,', channel:',ch,', position:',pos, ', correcting background')
                try:
                    phase_image = self.image_arrays[2][pos]['Trans']
                except KeyError:
                    phase_image = self.image_arrays[2][pos]['Phase']
                signal_image = self.image_arrays[2][pos][ch]
                bkg_cor_arrays[ch][pos] = self.back_sub(phase_image, signal_image, 25, 128, 60, False)[2]
        return bkg_cor_arrays
                
    
    def get_cell_mean_stats(self):
        cor_image_arrays = self.get_background_corrected_images()
        cell_list = []
        for pos in self.mask_arrays.keys():
            cell_list += list(self.mask_arrays[pos][0].keys())
        mean_df = pd.DataFrame()
        mean_df['cell_id'] = cell_list
        for ch in cor_image_arrays:
            area_dict  = {}
            mean_fluor_dict = {}
            std_fluor_dict = {}
            per_fluor_dict = {}
            total_fluor_dict = {}
            position_dict = {}
            for pos in cor_image_arrays[ch]:
                for cl in self.mask_arrays[pos][0]:
                    cell_mask = self.mask_arrays[pos][0][cl]
                    cell_pad = self.mask_arrays[pos][1][cl]
                    crop_signal_image =  cor_image_arrays[ch][pos][cell_pad[1]:cell_pad[3], cell_pad[0]:cell_pad[2]]
#                    plt.imshow(crop_signal_image)
#                    plt.show()
#                    plt.imshow(cell_mask)
#                    plt.show()
#                    input()
                    if 'cell_area_px' not in mean_df:
                        area_dict[cl] = np.nonzero(cell_mask)[0].shape[0]
                    mean_fluor_dict[cl] = crop_signal_image[np.nonzero(cell_mask)].mean()
                    total_fluor_dict[cl] = crop_signal_image[np.nonzero(cell_mask)].sum()
                    std_fluor_dict[cl] = crop_signal_image[np.nonzero(cell_mask)].std()
                    per_fluor_dict[cl] = np.percentile(crop_signal_image[np.nonzero(cell_mask)], 90)
                    position_dict[cl] = pos
            mean_df[ch+'_mean'] = mean_df.cell_id.map(mean_fluor_dict)
            mean_df[ch+'_total'] = mean_df.cell_id.map(total_fluor_dict)
            mean_df[ch+'_std'] = mean_df.cell_id.map(std_fluor_dict)
            mean_df[ch+'_90'] = mean_df.cell_id.map(per_fluor_dict)
            mean_df['experiment'] = self.experiment
            mean_df['position'] = mean_df.cell_id.map(position_dict)
            if 'cell_area_px' not in mean_df:
                mean_df['cell_area_px'] = mean_df.cell_id.map(area_dict)
        mean_df.to_pickle(self.save_path+'/'+self.experiment+'_mean_df', compression='zip')
        return mean_df
        
    
    def segment_nucleoids(self, channel, test_segmentation=False, hard_threshold = 10000, parameters=[1,1000,90,2,9,-2,0], min_nucleoid_size=10):
        """
        This function can be used to segment the nucleloids in a heterogeneous population.
        
        Parameters:
            channel: string - the channel of the HupA or DAPI fluorescence
            test_segmentation: binary - set to True to inspet the nucleoid segregation
            hard_threshold: integer - a hard fluorescence threshold above which the brighter nucleoid objects are considered
            parameters: list - the parameters of the LoG-adaptive filter
        Reutnrs:
            A dictionary with the cellID as key (string) and the nucleoid mask as value (2D binary numpy array)
        """
        nucleoid_mask = {}
        number_of_nucleoids = {}
        nucleoid_area = {}
        nc_ratio = {}
        
        # parameters = [2,1000,90,2,3,-9,0]
        cor_image_arrays = self.get_background_corrected_images()
        cell_list = []
        for pos in self.mask_arrays.keys():
            cell_list += list(self.mask_arrays[pos][0].keys())
        
        for pos in cor_image_arrays[channel]:
            log_masked_image = self.log_adaptive_filter(cor_image_arrays[channel][pos], parameters)[2]
            hard_masked_image = self.hard_threshold(cor_image_arrays[channel][pos], hard_threshold)
            
            masked_image = log_masked_image * hard_masked_image
            
            for cl in self.mask_arrays[pos][0]:
                # print(cl)
                cell_mask = self.mask_arrays[pos][0][cl]
                cell_pad = self.mask_arrays[pos][1][cl]
                crop_mask_image = masked_image[cell_pad[1]:cell_pad[3], cell_pad[0]:cell_pad[2]]
                crop_signal_image =  cor_image_arrays[channel][pos][cell_pad[1]:cell_pad[3], cell_pad[0]:cell_pad[2]]
                
                cell_pixels = crop_signal_image[np.nonzero(cell_mask)]
                otsu_mask = crop_signal_image>threshold_otsu(cell_pixels)
        
                nuc_mask = crop_mask_image*cell_mask.astype(int)
                nuc_mask = nuc_mask * otsu_mask
                nuc_mask = remove_small_objects(nuc_mask.astype(bool), min_nucleoid_size).astype(int)
                
                nucleoid_mask[cl] = nuc_mask
                number_of_nucleoids[cl] = label(nuc_mask).max()
                nucleoid_area[cl] = nuc_mask[np.nonzero(nuc_mask)].shape[0]
                nc_ratio[cl] = nucleoid_area[cl] / cell_mask[np.nonzero(cell_mask)].shape[0]
                
                if test_segmentation ==True:
                    plt.imshow(crop_signal_image*cell_mask)
                    plt.show()
                    plt.imshow(nuc_mask)
                    plt.show()
                    print(nc_ratio[cl], number_of_nucleoids[cl])
                    input()
        
        nucleoid_df = pd.DataFrame()
        nucleoid_df['cell_id'] = list(nucleoid_mask.keys())
        nucleoid_df['nucleoid_number'] = nucleoid_df.cell_id.map(number_of_nucleoids)
        nucleoid_df['nucleoid_area'] = nucleoid_df.cell_id.map(nucleoid_area)
        nucleoid_df['nc_ratio'] = nucleoid_df.cell_id.map(nc_ratio)
        
        nucleoid_df.to_pickle(self.save_path+'/'+self.experiment+'_nucleoid_df', compression='zip')
        with open(self.save_path+'/'+self.experiment+'_nucleoid_masks', 'wb') as handle:
            pickle.dump(nucleoid_mask, handle)
        
        return nucleoid_df
                

    
    def get_oned_coordinates(self, cell_mask, medial_axis_df, half_window): 

        cell_mask_df = pd.DataFrame()
        cell_mask_df['x'] = np.nonzero(cell_mask)[1]
        cell_mask_df['y'] = np.nonzero(cell_mask)[0]
    #    cell_mask_df['z'] = fluor_image[np.nonzero(cell_mask)]
    
        def get_pixel_projection(pixel_x, pixel_y, medial_axis_df, half_window):
            
            medial_axis_df['pixel_distance'] = np.sqrt((medial_axis_df.cropped_x-pixel_x)**2+(medial_axis_df.cropped_y-pixel_y)**2)
            min_df = medial_axis_df[medial_axis_df.pixel_distance == medial_axis_df.pixel_distance.min()]
            min_arch_centered_length = min_df.arch_length_centered.values[0]
            min_arch_scaled_length =  min_df.arch_length_scaled.values[0]
            min_distance_abs = min_df.pixel_distance.values[0]
            min_index = min_df.index.values[0]
            medial_axis_coords = (min_df.cropped_x.values[0], min_df.cropped_y.values[0])
            
            def get_relative_distance(min_distance_abs, medial_axis_df, min_index, medial_axis_coords, pixel_x, pixel_y, half_window):
        
                if min_index>=half_window and min_index<medial_axis_df.index.max()-half_window:
                    index_range = (min_index-half_window, min_index+half_window)
                elif min_index<half_window and min_index<medial_axis_df.index.max()-half_window:
                    index_range = (0, min_index+half_window)
                elif min_index>=half_window and min_index>=medial_axis_df.index.max()-half_window:
                    index_range = (min_index-half_window, medial_axis_df.index.max())
                
                delta_x = (medial_axis_df.iloc[index_range[1]].cropped_x -  medial_axis_df.iloc[index_range[0]].cropped_x)
                delta_y = (medial_axis_df.iloc[index_range[1]].cropped_y -  medial_axis_df.iloc[index_range[0]].cropped_y)
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
        
    
    def apply_oned_coordinates(self):
        
        
        with open(self.save_path+'/'+self.experiment+'_medial_axis_dict', 'rb') as handle:
            medial_axis_dict = pickle.load(handle)
        
        oned_coords_dict = {}
        cell_length_dict = {}
        
        cor_image_arrays = self.get_background_corrected_images()
        
        for pos in self.mask_arrays:
            for cl in self.mask_arrays[pos][0]:
                if cl in medial_axis_dict[pos]:
                    print(cl, 'mapping 1D pixel coordinates')
                    cell_mask = self.mask_arrays[pos][0][cl]
                    cell_pad = self.mask_arrays[pos][1][cl]
                    medial_axis_df = medial_axis_dict[pos][cl][0]
                    cell_mask_df = self.get_oned_coordinates(cell_mask, medial_axis_df, half_window=10)
                    for ch in cor_image_arrays:
                        crop_signal_image =  cor_image_arrays[ch][pos][cell_pad[1]:cell_pad[3], cell_pad[0]:cell_pad[2]]
                        cell_mask_df[ch+'_fluor'] = crop_signal_image[np.nonzero(cell_mask)]
                    oned_coords_dict[cl] = cell_mask_df
                    cell_length_dict[cl] = medial_axis_df.arch_length_centered.max()*2

        pickle.dump(oned_coords_dict, open(self.save_path+'/'+self.experiment+'_oned_coords_dict', 'wb'))
        pickle.dump(cell_length_dict, open(self.save_path+'/'+self.experiment+'_cell_length_dict', 'wb'))
        
        return oned_coords_dict
    
    
    #--------- IMAGE FILTERS ---------#
    # Customized image filters used to sharpen or smooth, as well as threshold images.
    def hard_threshold(self, image, threshold):
        
        return image>threshold
    
    
    def log_adaptive_filter(self, image, parameters): 
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
        
        # image = particle.back_sub(particle.particle_images[0], show=True)[0]
        # parameters = [4, 1000, 99.5] for testing
        image[image<0]=0
        # LoG filter with hard threshold
        # image_smoothed = filters.gaussian(image, parameters[0])
        image_smoothed = ndimage.gaussian_filter(image, sigma = parameters[0])
        image_laplace = filters.laplace(image_smoothed, parameters[1])
        image_pixel_intensities = image_laplace.ravel()
        sorted_pixel_intensities = np.sort(image_pixel_intensities)
        pixel_intensity_threshold = sorted_pixel_intensities[int((parameters[2]/100)*len(sorted_pixel_intensities))]
        log_image = image_laplace > pixel_intensity_threshold
        # plt.figure(figsize=(20,20))
        # plt.imshow(masked_image)
        # # plt.colorbar()
        # plt.show()
        # Adaptice threshold
        # image_smoothed_2 = filters.gaussian(image, parameters[3])
        image_smoothed_2 = ndimage.gaussian_filter(image, sigma=parameters[3])
        adaptive_threshold = threshold_local(image_smoothed_2, block_size=parameters[4], offset=parameters[5])
        # adaptive_threshold = threshold_local(image, block_size=parameters[4], offset=parameters[5])
        adaptively_masked_image = image_smoothed_2 > adaptive_threshold
        # plt.figure(figsize=(20,20))
        # plt.imshow(adaptively_masked_image)
        # # plt.colorbar()
        # plt.show()
        masked_image = adaptively_masked_image * log_image
            # erode the image if the erosion iterations are higher than 1
        if parameters[6] > 0:
            adaptively_masked_image = scipy.ndimage.morphology.binary_erosion(masked_image, iterations = parameters[6])
            final_image = adaptively_masked_image
        elif parameters[6] == 0:
            final_image = masked_image

        return log_image, adaptively_masked_image, final_image
    

    def log_filter(self, image, parameters): 
        """
        This fucntion constructs an LoG filer (Laplace of Gaussian) as well as an adaptive filter to segment particles
        
        Parameters
        ----------
        image: numpy array - the image to be filtered and thresholded (usually this is the background subtracted image)
        particle_detection_parameters: list - the parameters for the LoG/adaptive filter
            [0] Smoothing factor for the Gaussian filter (https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.gaussian)
            [1] Laplace threshold (https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.laplace)
            [2] Hard threshold on the LoG image as a percentile of the brightest pixels
        
        Returns
        -------
        The thresholded image after the LoG filter
        """
        
        # image = particle.back_sub(particle.particle_images[0], show=True)[0]
        # parameters = [4, 1000, 99.5] for testing
        image[image<0]=0
        # LoG filter with hard threshold
        # image_smoothed = filters.gaussian(image, parameters[0])
        image_smoothed = ndimage.gaussian_filter(image, sigma = parameters[0])
        image_laplace = filters.laplace(image_smoothed, parameters[1])
        image_pixel_intensities = image_laplace.ravel()
        sorted_pixel_intensities = np.sort(image_pixel_intensities)
        pixel_intensity_threshold = sorted_pixel_intensities[int((parameters[2]/100)*len(sorted_pixel_intensities))]
        masked_image = image_laplace > pixel_intensity_threshold
        return masked_image
    
    