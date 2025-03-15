# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 13:24:09 2022

@author: Alexandros Papagiannakis, Christine Jacobs-Wagner lab, Sarafan ChEM-H, Stanford University, 2022
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skimage import io
from itertools import combinations
from skimage.filters import threshold_otsu
from scipy import ndimage
from shapely.geometry import LineString
from PIL import Image, ImageDraw
import pickle
from skimage.morphology import medial_axis
import os
from scipy.io import loadmat


class oufti_snapshots_GrayGovers(object):
    """
    Developer: Alexandros Papagiannakis, Christine Jacobs-Wagner lab, Stanford University, 2022
    
    This class contains all the functions for analyzing the polysome and nucleoid fluorescence statistics using Oufti cell meshes.
    
    Functions included in this class:
        __init__
            tiff_to_arrays
            cell_list_to_python (conversion from Oufti)
            oufti_meshes_to_masks
        check_cell_meshes
        angle_from_slope
        get_medial_axis
            correct_angle_difference
            get_next_position
            recursive_medial_axis
        run_medial_axis
        check_cell_segmentation
        cell_free_bkg_estimation
        back_sub
        get_background_corrected_images
        get_cell_mean_stats
        get_oned_coordinates
            get_pixel_projection
                get_relative_distance
        apply_oned_coordinates
    """
    def __init__(self, experiment_path, save_path):
        """
        The class is iniialized.
        
        Parameters
        ----------
        ___File paths___
        experiment_path: The path of all the tiff images and the oufti paths
        save_path: path where the results are saved
        
        """
        
        def tiff_to_arrays(images_path):
            """
            This function is used to read all the tiff images in the specified folder.
            The tiff images are saved in the dictionary as numpy arrays.
            """
            
            images_list = os.listdir(images_path)
            ch = int(images_list[0][-5])
            self.image_arrays_dict[ch]={}
            for img in images_list:
                if 'xy' in img:
                    try:
                        pos = int(img[img.find('xy')+2:img.find('xy')+4])-1
                    except ValueError:
                        xy_index = img.find('xy', img.find('xy')+1)
                        pos = int(img[xy_index+2:xy_index+4])-1
                        
                    print(ch, pos)
                    self.image_arrays_dict[ch][pos]=io.imread(images_path+'/'+img)
                
            
        def cell_list_to_python(oufti_path, position):
            """
            This function incoorporates the cell list data from oufti to python.
            It searches into the cell list for data fields which include 4 columns, corresponding to the cell mesh data.
            
            The cell meshes are broken down into x,y coordinates and the medial axis is calculated.
            
            Input:
                cell_list_path: string - the path of the oufti cell list .mat file
                position: non-negative integer - the position to be analyzed (this is useful if batches of images are analyzed)
            
            Returns:
                cell_mesh_dict: a dictionary which inlcudes the mesh x,y coordinates (value - tuple) for each cell (key - integer)
            """
            
            
            oufti = loadmat(oufti_path)
            cell_list = oufti['cellList']
#            cell_list_N = oufti['cellListN']
            
            # get the cell IDs from oufti
            cell_ID_oufti = cell_list[0][0][1][0][position][0]
            
#            n = int(cell_list_N[0][position])
            n = len(cell_ID_oufti)
        #    n_2 = 10  ## used to test the exeption
            
            cell_mesh_dict = {}
#            cell_medial_axis_dict = {}
            
            for cell in range(n):
                index = 0
                shape_1 = 0
                # looking for a field with 4 columns. This field corresponds to the cell meshes.
                while shape_1 != 4:
                    cell_spine = cell_list[0][0][0][0][position][0][cell][0][0][index]
                    shape_1 = cell_spine.shape[1]
                    index += 1
#                    if shape_1 == 4:
#                            print(cell_spine, cell, position)
                # Create empty lists to append the cell mesh and medial axis coordinates
                polygon_oufti_x = []
                polygon_oufti_y = []
#                medial_axis_oufti_x = []
#                medial_axis_oufti_y = []
                for i in range(cell_spine.shape[0]):
                    # WARNING
                    # subtract 1 from all the coordinates to match the zero index of the images
                    polygon_oufti_x.append(cell_spine[i][0]-1)
                    polygon_oufti_y.append(cell_spine[i][1]-1)
#                    medial_axis_oufti_x.append((cell_spine[i][0]-1+cell_spine[i][2]-1)/2)
#                    medial_axis_oufti_y.append((cell_spine[i][1]-1+cell_spine[i][3]-1)/2)
                for i in range(1,cell_spine.shape[0]+1):
                    polygon_oufti_x.append(cell_spine[-i][2]-1)
                    polygon_oufti_y.append(cell_spine[-i][3]-1)
                
                
                # create a new cell ID which includes the experiment and the XY position, as well as the original oufit cell ID.
                if position < 9:
                    cell_ID = self.experiment+'_xy0'+str(position+1)+'_cell_'+str(int(cell_ID_oufti[cell]))
                elif position >= 9:
                    cell_ID = self.experiment+'_xy'+str(position+1)+'_cell_'+str(int(cell_ID_oufti[cell]))
                
                # add the cell specific coordinates into the dictionary
                cell_mesh_dict[cell_ID] = (polygon_oufti_x, polygon_oufti_y)
#                cell_medial_axis_dict[cell_ID] = (medial_axis_oufti_x, medial_axis_oufti_y)
            
            return cell_mesh_dict



        def oufti_meshes_to_masks(meshes_dictionary, pad):
            """
            This function will be used to create cell masks and crop the images around the cell masks
            
            Input:
                meshes_dictionary: dictionary with the cell meshes from the cell_list_to_python[0] function
                pad: non-negative integer - the number of pixels around the max/min coordinates of the cell polygon surrounding the cropped image.
            Returns:
                [0] cropped_cell_masks: dictionary linking the cropped cell masks (values) to each cell ID (keys)
                [1] pads: dictionary linking the cropping pads (values) to each cell ID (keys)
            """
            
            # start empty dictionaries to store the masked and cropped images
            cropped_cell_masks = {}
            pads = {}
            
            # create a masked image
            masked_image = Image.new('L', self.sensor, 0)
            
            for cell in meshes_dictionary:    
                try:      
                    polygon_coords = []
#                    cell_mesh = np.around(meshes_dictionary[cell]).astype(int)
                    cell_mesh = np.array(meshes_dictionary[cell]).astype(int)
                    # iterate over the two columns of x and y mesh coordinates and create one column of (x,y) tuples.
                    for coord_index in range(len(cell_mesh[0])):
                        polygon_coords.append((cell_mesh[0][coord_index], cell_mesh[1][coord_index]))

                    # get the extremites of the polygon
                    minx,miny,maxx,maxy = np.min(cell_mesh[0]),np.min(cell_mesh[1]),np.max(cell_mesh[0]),np.max(cell_mesh[1])
                    # get the crop pad around the cell
                    pad_list = np.array([int(minx-pad), int(miny-pad), int(maxx+pad), int(maxy+pad)])
                    # covert negative values in a pad into zero
                    if all(it>=0 for it in pad_list) and all(it<2048 for it in pad_list):
                        
                        pads[cell] = (pad_list[0], pad_list[1], pad_list[2], pad_list[3])
                        
                        # Create a dark image with 0 pixel intensities and size equal to the original camera sensor
                        draw_image = Image.new('L', self.sensor, 0)
                        # Draw the mesh in the dark image
                        ImageDraw.Draw(draw_image).polygon(polygon_coords, outline=0, fill=1)
                        # Link each mask to the cell in the cell_masks dictionary
                        cell_mask = np.array(draw_image)
                        # Create a dictionary which links each cell with the cropped mask
                        cropped_cell_masks[cell] = cell_mask[pads[cell][1]:pads[cell][3], pads[cell][0]:pads[cell][2]]
                        # Add the cell masks into the image
                        masked_image+=cell_mask
                    else:
                        print('This cell is at the edge of the field of view and the cropping pad is out of bounds.')
                        print(cell, 'is excluded from the analysis')
            
                
                except ValueError:
                    print('Cell mesh:', cell, 'is not correct and was not included in the cell dictionaries')
                    print('Please correct this cell ID in oufti and the respective XY position')
                    print('This cell ID is deleted')
                    


            plt.figure(figsize=(20,20))
            plt.imshow(masked_image,cmap='gray')
            plt.show()
            
            
            return cropped_cell_masks, pads
        
        
        self.experiment_path = experiment_path
        self.image_arrays_dict = {}
        
        cell_list = [s for s in os.listdir(self.experiment_path) if "cellList.mat" in s][0]
        self.oufti_path = self.experiment_path + '/' + cell_list
        self.experiment = cell_list[0:cell_list.find('_cellList.mat')]
        phase_images =  [s for s in os.listdir(self.experiment_path) if "c1" in s][0]
        ribo_images =  [s for s in os.listdir(self.experiment_path) if "c2" in s][0]
        dna_images =  [s for s in os.listdir(self.experiment_path) if "c3" in s][0]
            
        tiff_to_arrays(self.experiment_path+'/'+phase_images)
        tiff_to_arrays(self.experiment_path+'/'+ribo_images)
        tiff_to_arrays(self.experiment_path+'/'+dna_images)
        
        self.sensor = self.image_arrays_dict[1][0].shape
        
        positions = []
        for img in os.listdir(self.experiment_path+'/'+phase_images):
            try:
                positions.append(int(img[img.find('xy')+2:img.find('xy')+4])-1)
            except ValueError:
                xy_index = img.find('xy', img.find('xy')+1)
                positions.append(int(img[xy_index+2:xy_index+4])-1)
        
        cell_mesh_dict = {}
        cropped_cell_masks = {}
        pads = {}
        for position in positions:
            cell_mesh_dict[position] = cell_list_to_python(self.oufti_path, position)
            cropped_cell_masks[position], pads[position] = oufti_meshes_to_masks(cell_mesh_dict[position], 5)
        
        self.cell_mesh_dict = cell_mesh_dict
        self.cropped_cell_masks = cropped_cell_masks
        self.pads = pads
        self.save_path = save_path
        
        
    def check_cell_meshes(self):
        
        for pos in self.cell_mesh_dict:
            plt.figure(figsize=(15,15))
            plt.imshow(self.image_arrays_dict[1][pos], cmap='gray')
            for cl in self.cell_mesh_dict[pos]:
                plt.plot(self.cell_mesh_dict[pos][cl][0], self.cell_mesh_dict[pos][cl][1], c='red')
            plt.show()
        

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
        for pos in self.cropped_cell_masks:
            medial_axis_dict[pos]={}
            for cell in self.cropped_cell_masks[pos]:
                cropped_cell_mask = self.cropped_cell_masks[pos][cell]
                cell_pad = self.pads[pos][cell]
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
        Local background estimation outside the cell regions (outside the otsu thresholded dilated cell area).
        The average background within a square of size equal to the 'estimation_step' is used to fill in the cellular background
        (within the dilated cell mask) and then a Gaussian smoothing is applied.
        
        Parameters
        ----------
        phsae_image: numpy.array - the phase contrast image
        signal_image: numpy.array - the image to be corrected
        dilation: non-negative integer - the number of dilation rounds for the cell mask
        estimation_step: positive_integer - the size of the square edge used for average background estimation
        smoothing_sigma: non-negative integer - the smoothing factor of the cell free background
        show: binary - True if the user wants to visualize the 2D surface fit
        
        Returns
        -------
        [0] The average background
        [1] The background corrected image (after subtracting the smoothed cell free signal)
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
#        fluor_channels = self.channels.copy()
#        fluor_channels.remove('Trans')
        for ch in range(2,4):
            bkg_cor_arrays[ch]={}
            for pos in self.image_arrays_dict[ch].keys():
                print(self.experiment,', channel:',ch,', position:',pos, ', correcting background')
                phase_image = self.image_arrays_dict[1][pos]
                signal_image = self.image_arrays_dict[ch][pos]
                bkg_cor_arrays[ch][pos] = self.back_sub(phase_image, signal_image, 25, 128, 60, False)[2]
        return bkg_cor_arrays
                
    
    def get_cell_mean_stats(self):
        cor_image_arrays = self.get_background_corrected_images()
        cell_list = []
        for pos in self.cropped_cell_masks.keys():
            cell_list += list(self.cropped_cell_masks[pos].keys())
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
                for cl in self.cropped_cell_masks[pos]:
                    cell_mask = self.cropped_cell_masks[pos][cl]
                    cell_pad = self.pads[pos][cl]
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
            mean_df['ch'+str(ch)+'_mean'] = mean_df.cell_id.map(mean_fluor_dict)
            mean_df['ch'+str(ch)+'_total'] = mean_df.cell_id.map(total_fluor_dict)
            mean_df['ch'+str(ch)+'_std'] = mean_df.cell_id.map(std_fluor_dict)
            mean_df['ch'+str(ch)+'_90'] = mean_df.cell_id.map(per_fluor_dict)
            mean_df['experiment'] = self.experiment
            mean_df['position'] = mean_df.cell_id.map(position_dict)
            if 'cell_area_px' not in mean_df:
                mean_df['cell_area_px'] = mean_df.cell_id.map(area_dict)
        mean_df.to_pickle(self.save_path+'/'+self.experiment+'_mean_df', compression='zip')
        return mean_df
    
    
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
        
        for pos in self.cropped_cell_masks:
            for cl in self.cropped_cell_masks[pos]:
                if cl in medial_axis_dict[pos]:
                    print(cl, 'mapping 1D pixel coordinates')
                    cell_mask = self.cropped_cell_masks[pos][cl]
                    cell_pad = self.pads[pos][cl]
                    medial_axis_df = medial_axis_dict[pos][cl][0]
                    cell_mask_df = self.get_oned_coordinates(cell_mask, medial_axis_df, half_window=10)
                    for ch in cor_image_arrays:
                        crop_signal_image =  cor_image_arrays[ch][pos][cell_pad[1]:cell_pad[3], cell_pad[0]:cell_pad[2]]
                        cell_mask_df['ch'+str(ch)+'_fluor'] = crop_signal_image[np.nonzero(cell_mask)]
                    oned_coords_dict[cl] = cell_mask_df
                    cell_length_dict[cl] = medial_axis_df.arch_length_centered.max()*2

        pickle.dump(oned_coords_dict, open(self.save_path+'/'+self.experiment+'_oned_coords_dict', 'wb'))
        pickle.dump(cell_length_dict, open(self.save_path+'/'+self.experiment+'_cell_length_dict', 'wb'))
        
        return oned_coords_dict
    
    
    
    
        
            
                
        
        
        
        