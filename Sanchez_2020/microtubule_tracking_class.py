

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 19:27:10 2020
@author: Alexandros Papagiannakis, Christine-Jacobs Wagner lab, Stanford University, 2020
"""

# libraries
from nd2reader import ND2Reader # library to open .nd2 images
# https://rbnvrw.github.io/nd2reader/index.html  
import matplotlib.pyplot as plt # library to plot
import os # library to manipulate and explore directories
import numpy as np # numpy arrays - data object
from skimage.filters import threshold_otsu, threshold_local, gaussian, laplace
from skimage import morphology
import scipy
from skimage.measure import label, regionprops
import matplotlib
import pandas as pd
import pickle
from numpy import arccos, array, dot, pi, cross
from numpy.linalg import det, norm
from scipy.interpolate import UnivariateSpline
#matplotlib.use('Qt5Agg')






class microtub_tracking(object):
    """
    A class that is used to track microtubules in C elegans embryos.
    Developer: Alexandros Papagiannakis, Christine-Jacobs Wagner group, Stanford University 2020.
    
    This class was developed for and used in the publication:
        Authors:
        Ariana D. Sanchez, Tess C. Branon, Lauren E. Cote, Alexandros Papagiannakis, 
        Xing Liang, Melissa A. Pickett, Kang Shen, Christine Jacobs-Wagner, Alice Y. Ting, Jessica L. Feldman
        Title:
        Proximity labeling at non-centrosomal microtubule-organizing centers reveals VAB-10B and WDR-62 as distinct microtubule regulators
        
        Journal:
        bioRxiv
            
        Year:
        2020
        
        DOI:
        https://doi.org/10.1101/2020.08.29.272369.
    
            
    This class includes a number of filters and functions to segment and track microtubules in C elegans embryos:
        Filters:
            log_filter - LoG filter
            serial_otsu_filters - not used 
            log_filter_adaptive - not used 
            adaptive_threshold_smoothed - adaptive thredholding on the smoothed image
            adaptive_threshold - not used
            embryo_and_gut_filter - adaptive filter used to get the embryo mask and the mask of its intestine
        Functions:
            nd2_snapshots_to_array - converts the nd2 file into numpy arrays
            check_segmentation_parameters - helps determine the comet segmentation parameters
            segment_microtubules - runs segementation on a single time-point
            run_segmentation - runs segmentation on all time-points
            show_segmentation - use to check the segmentation
            microtubule_tracking - tracks the microtubules by linking adhacent comets in subsequent time-points
            show_microtubule_trajectories - shows all the microtubule trajectories
            get_microtubule_speed_and_angles - estimates the speed and the angle of the trajectories
            max_projection - returns the maximum projection of the stream acquisition images
            get_central_line - estimates the central line
            get_gut_mask - constructs the gut mask inside the embryo
            get_embryo_mask - constructs the mask of the embryo
        The following data are saved in the specified path:
            The particle segmentation parameters (dictionary)
            The fluorescence images of the microtubules with the segmented spots (.jpeg)
            A dataframe with the coordinates of the segmented microtubules (pandas dataframe)
            A dataframe with the linked comet positions and the trajectory IDs (pandas dataframe)
            A dataframe with the angle and speed of the microtubule trajectories (pandas dataframe)
            The central line fitted to the manually selected points (list of (x,y) tuples)
            The gut mask (binary numpy array)
            The embryo mask (binary numpy array)
    """
    
    
    def __init__(self, image_path, snapshots_path, save_path, experiment):
        """
        Initializes the class.
        
        The calss initialization imports all the images from the nd2 file and 
        the respective metadata to determine a number of parameters.
        
        Input
            image_path: string - the path of the nd2 documents with the microtubule and phase contrast images
            save_path: string - a path where all the results will be saved
            experiment: sting - an experiment ID which is used to 
            
        Returns
            self.metadata: dictionary which includes the image metadata
            self.frames: a list with a number of frames
            self.interval: the time interval between frames in seconds (rounded up to the first decimal eg 1.5 sec for 1.4942901 sec)
            self.scale: the scale of the images in  μm/px rounded to the second decimal (eg 0.11 μm/px for 0.10811640 μm/px)
            self.n_frames: positive integer - the number of frames
            se;f.channels: a list which includes the channels (e.g. ['TEST 488', 'TEST 405'])
            self.image_arrays: a dictionary which includes all the numpy arrays for the images in all channels.
                The first key of the dictionary corresponds to the channel and the second key to the frame.
                To get the frist frame of the 'TEST 488' channel type: image_arrays['TEST 488'][0]  - EXAMPLE
        """
        
    
        self.image_path = image_path
        self.save_path = save_path
        self.experiment = experiment

        if os.path.isdir(self.save_path) == False:
            raise ValueError('The specified save path is not an valid directory')
        if os.path.exists(self.image_path) == False:
            raise ValueError('The specified image path is not an valid directory')
            
   
        # export the images as numpy arrays for every single channel and timepoint and get the .nd2 metadata
    
        def nd2_to_array(image_path):
            """
            Developer: Alexandros Papagiannakis, CJW group, Stanford University 2020.
            
            This function is used to convert .nd2 images to numpy arrays.
            Input:
                image_path - string: the path of the .nd2 file
            Returns:
                [0] the .nd2 metadata and images (from the ND2Reader)
                [1] a dictionary which contains the images as numpy arrays organized by channels (1st key) and timepoints (2md key)
            """
        
            # The path of the .nd2 file 
            images = ND2Reader(image_path) 
            print('metadata:',images.metadata)
            print('dimensions:',images.sizes)
            
            # get the channels and frames from the .nd2 metadata
            channels = images.metadata['channels']
#            frames = images.metadata['frames']
            
            
            # iterate over the channels and timepoints        
            # create a dictionary were the image numpy arrays will be stored for each channel and timepoint (dictionary within dictionary)
            image_arrays = {}
            image_arrays[channels[0]] = {}
            # iterate over the channels to create a dictionary within each dictionary
            
            f = 0
            
            for img in images:
                
                print('channel:', channels[0],', timepoint:', f)
                        
                #print(fov) # Frame containing one image
#                plt.imshow(np.array(img))
#                plt.show()
                image_arrays[channels[0]][f] = np.array(img)
                f+=1

            return images, image_arrays
        
        
        
        images = nd2_to_array(self.image_path)
        
        self.metadata = images[0].metadata
        self.frames = self.metadata['frames']
        self.interval = round(self.metadata['experiment']['loops'][0]['sampling_interval']/1000,1)
        self.scale = round(self.metadata['pixel_microns'],2)
        self.n_frames = self.metadata['num_frames']
        self.channels = self.metadata['channels']
        self.image_arrays = images[1]
        self.snapshots_path = snapshots_path
        
        
    
    
    def nd2_snapshots_to_array(self, image_path):
        """
        This function is used to convert .nd2 images to numpy arrays.
        Input:
            image_path - string: the path of the .nd2 file
        Returns:
            [0] the .nd2 metadata and images (from the ND2Reader)
            [1] a dictionary which contains the images as numpy arrays organized by channels (1st key) and timepoints (2md key)
        """
    
        # The path of the .nd2 file 
        images = ND2Reader(image_path) 
        print('metadata:',images.metadata)
        print('dimensions:',images.sizes)
        
        # get the channels and frames from the .nd2 metadata
        channels = images.metadata['channels']
        
        # iterate over the channels and timepoints        
        # create a dictionary were the image numpy arrays will be stored for each channel and timepoint (dictionary within dictionary)
        image_arrays = {}
        # iterate over the channels to create a dictionary within each dictionary
        f = 0
    
        for img in images:
            
            print('channel:', channels[f])
                    
            #print(fov) # Frame containing one image
            plt.imshow(np.array(img))
            plt.show()
            image_arrays[channels[f]] = np.array(img)
            f+=1

        return images, image_arrays
    
    
    
    
    def log_filter(self, image, smoothing_factor, laplace_factor, percentile):
        """
        This is a Lacplace of Gaussian filter
        1. The image is smoothed using a Gaussian filter
        
        2. Then a Laplace filter is applied to find fast changes in fluorescence intensity (measures the derivative of intensity changes)
        
        3. A hard threshold is applied on the image after the LoG filter is applied to select for the brigthest spots. 
        This hard threshold is provided as a percentile of the brightest pixels. 
        
        Input:
            image: numpy array of the fluorescence image
            smoothing_factor: the factor of the Gaussian smoothing. This is the standard deviation of fitted gaussian kernel.
            laplace_factor: Size of the discreet Laplace operator.
            percentile: The percentile of the brightest pixels of the LoG image to be considered. Hard threshold.
        Returns:
            [0] The LoG image array
            [1] The thresholded LoG image array - mask
            [2] The augmented image array - inreased contrast around the microtubule edges
        """
        image_smoothed = scipy.ndimage.gaussian_filter(image, sigma=smoothing_factor)
#        image_smoothed = scipy.ndimage.gaussian_filter(image, sigma=smoothing_factor)
        img_laplace = laplace(image_smoothed, laplace_factor)
        
#        dimmest_pixel = min(img_laplace.ravel())
#        brightest_pixel = max(img_laplace.ravel())
        
        #     print(dimmest_pixel)
        #     print(brightest_pixel)
        
        #     plt.hist(img_laplace.ravel(), bins=np.arange(dimmest_pixel,brightest_pixel,0.00001))
        #     plt.show()
        
        
#             plt.figure(figsize=(10,10))
#             plt.imshow(img_laplace*image_smoothed, cmap='gray')
#             plt.clim(-0.000001,0.00008)
#             plt.show()
#        
        
        sorted_data = img_laplace.ravel()
        sorted_data = np.sort(sorted_data)
        
        percentile_intensity = sorted_data[int(len(sorted_data)*(percentile/100))]
        
        thresholded_image = img_laplace > percentile_intensity
#             plt.figure(figsize=(10,10))
#             plt.imshow(thresholded_image, cmap='gray')
#             plt.clim(0,1)
#             plt.show()
       
#        gaussian_image = scipy.ndimage.gaussian_filter(image, sigma=smoothing_factor)
        
        image_augmented = image + thresholded_image * image
#             plt.figure(figsize=(10,10))
#             plt.imshow(image_augmented, cmap='gray')
##             plt.clim(0,1)
#             plt.show()
        
        return img_laplace, thresholded_image, image_augmented



    def serial_otsu_filters(self, image, rounds, smoothing_factor):
        """
        This code applies serial otsu filters on a smoothed image.
        
        Input:
            image: numpy array of the fluorescence image
            smoothing_factor: the factor of the Gaussian smoothing. This is the standard deviation of fitted gaussian kernel.
            rounds: positive integer - the number of otsu filtering rounds.
        Returns:
            [0] The thresholded image after the consecutive otsu filters - mask
        """
        for rnd in range(rounds):
            if rnd == 0:
                
    #             plt.figure(figsize=(5,5))
    #             plt.hist(image.ravel(), bins=50, color='blue')
    #             plt.xlabel('pixel intensity')
    #             plt.ylabel('number of pixels')
    #             plt.title('Raw image')
    #             plt.show()
                
#                image_smoothed = gaussian(image, smoothing_factor)
                otsu = threshold_otsu(image)
    #             print('otsu threshold, round',(rnd+1),'-',otsu)
#                image_otsu_bkg = image <= otsu
                image_otsu = image > otsu
                
                
            elif rnd > 0:
                img = image_otsu * image
                
    #             plt.figure(figsize=(5,5))
    #             plt.hist(img.ravel(), bins=50, color='blue')
    #             plt.xlabel('pixel intensity')
    #             plt.ylabel('number of pixels')
    #             plt.title('Otsu round: '+str(rnd))
    #             plt.show()
                
                otsu = threshold_otsu(img[np.nonzero(img)])
    #             print('otsu threshold, round',(rnd+1),'-',otsu)
#                image_otsu_bkg = image <= otsu
                image_otsu = image > otsu
            

     
    #     plt.figure(figsize=(10,10))
    #     plt.imshow(image_otsu, cmap='gray')
    #     plt.clim(0,1)
    #     plt.show()
        
        return image_otsu



    def log_filter_adaptive(self, image, smoothing_factor, laplace_factor, adaptive_block_size, adaptive_offset):
        """
        This is an adaptive filter applied after a Lacplace of Gaussian (LoG) filter
        1. The image is smoothed using a Gaussian filter
        
        2. Then a Laplace filter is applied to find fast changes in fluorescence intensity (measures the derivative of intensity changes)
        
        3. An adaptive threshold is applied to the LoG image
        
        Input:
            image: numpy array of the fluorescence image
            smoothing_factor: the factor of the Gaussian smoothing. This is the standard deviation of fitted gaussian kernel.
            laplace_factor: Size of the discreet Laplace operator.
            adaptive_block_size: postivie odd integer - Size of pixel neighborhood which is used to calculate the threshold value.
            adaptive_offset: Constant subtracted from weighted mean of neighborhood to calculate the local threshold value. 
        Returns:
            [0] The LoG image array
            [1] The thresholded LoG and adaptively thresholded image - mask
        """
        image_smoothed = scipy.ndimage.gaussian_filter(image, sigma=smoothing_factor)
#        image_smoothed = scipy.ndimage.gaussian_filter(image, sigma=smoothing_factor)
        img_laplace = laplace(image_smoothed, laplace_factor)
    
    #     plt.hist(img_laplace.ravel(), bins=100)
    #     plt.show()
    
    
    #     plt.figure(figsize=(10,10))
    #     plt.imshow(img_laplace, cmap='gray')
    #     plt.clim(-0.000001,0.00008)
    #     plt.show()
        
        
        adaptive_thres = threshold_local(img_laplace, block_size=adaptive_block_size, offset=adaptive_offset)
        image_local = img_laplace > adaptive_thres
    
    #     plt.figure(figsize=(10,10))
    #     plt.imshow(image_local, cmap='gray')
    #     plt.clim(0,1)
    #     plt.show()
        

        
        return img_laplace, image_local



    def adaptive_threshold_smoothed(self, image, smoothing_factor, adaptive_block_size, adaptive_offset):
        """
        This is an adaptive filter applied after a Gaussian filter
        1. The image is smoothed using a Gaussian filter
        
        2. An adaptive threshold is applied to the LoG image
        
        Input:
            image: numpy array of the fluorescence image
            smoothing_factor: the factor of the Gaussian smoothing. This is the standard deviation of fitted gaussian kernel.
            adaptive_block_size: postivie odd integer - Size of pixel neighborhood which is used to calculate the threshold value.
            adaptive_offset: Constant subtracted from weighted mean of neighborhood to calculate the local threshold value. 
        Returns:
            [0] The smoothed image array
            [1] The thresholded LoG and adaptively thresholded image - mask
        """
        image_smoothed = scipy.ndimage.gaussian_filter(image, sigma=smoothing_factor)
#        image_smoothed = scipy.ndimage.gaussian_filter(image, sigma=smoothing_factor)
    
    #     plt.hist(image_smoothed.ravel(), bins=100)
    #     plt.show()
    
    #     plt.figure(figsize=(10,10))
    #     plt.imshow(image_smoothed, cmap='gray')
    #     plt.show()
        
        
        adaptive_thres = threshold_local(image_smoothed, block_size=adaptive_block_size, offset=adaptive_offset)
        image_local = image_smoothed > adaptive_thres
    
#         plt.figure(figsize=(10,10))
#         plt.imshow(image_local, cmap='gray')
#         plt.clim(0,1)
#         plt.show()
        
        
        
        return image_smoothed, image_local
    
    
    
    def adaptive_threshold(self, image, adaptive_block_size, adaptive_offset):
        """
        This is an adaptive filter applied after a Gaussian filter
        
        1. An adaptive threshold is applied to the LoG image
        
        Input:
            image: numpy array of the fluorescence image
            adaptive_block_size: postivie odd integer - Size of pixel neighborhood which is used to calculate the threshold value.
            adaptive_offset: Constant subtracted from weighted mean of neighborhood to calculate the local threshold value. 
        Returns:
            [1] The thresholded LoG and adaptively thresholded image - mask
        """

        
        
        adaptive_thres = threshold_local(image, block_size=adaptive_block_size, offset=adaptive_offset)
        image_local = image > adaptive_thres
    
#         plt.figure(figsize=(10,10))
#         plt.imshow(image_local, cmap='gray')
##         plt.clim(0,1)
#         plt.show()
#        
        
        return image_local
    
    
    
    def check_segmentation_parameters(self, frame, channel):
        """
        This function is used to check and save the spot segmentation parameters.
        The user is guided to select all the parameters for the the filters applied to distinguish the particles in the original image.
        
        At the end of the algorthm, when the user has selected a set of "good" parameters, these are saved in a dictionary in the specified folder (save_path).
        
        The following parameters need to be selected by the user:
            LoG filter paramters (gaussian_param, laplace_param, log_threshold)
            Adaptive threshold parameters (smoothing_factor, adaptive_block_size, adaptive_offset)
            The number of erosion iterations for the contrasted and adaptively thresholded mask (erosion_iterations)
            The maximum size of the spot in pixels (spot_size)
            The minimum aspect ratio (minor/major axis) of the spot
            LoG filter parameters for clustered spots (gaussian_param_zoom, laplace_param_zoom, log_threshold_zoom)
            
            Example of returned dictonary:
                {'log': (1.0, 40, 99.5),
                 'adaptive_threshold': (3.0, 5, 1e-09),
                 'erosion': 2,
                 'spot_area': 15,
                 'aspect_ratio': 1.0,
                 'log_zoom': (1.0, 20, 94.5)}
                
            In this example a set of reccomended parameters is included.
            The dictioanry is also save in the specified path (selg.save_path)
        """
        
        
        image = self.image_arrays[channel][frame]
        
        print('Selecting the parameters for the Laplace of Gaussian filter...')
        
        i = 0
        while i == 0:
            try:
                gaussian_param = float(input('Choose the gaussian smoothing factor (recommended: 1):'))
            except ValueError:
                print('please choose a number')
                gaussian_param = float(input('Choose the gaussian smoothing factor (recommended: 1):'))
            
            try:
                laplace_param = int(input('Choose the laplace filter factor (recommended: 40):'))
            except ValueError:
                print('please choose a number')
                laplace_param = int(input('Choose the laplace filter factor (recommended: 40):'))
            
            try:
                log_threshold = float(input('Choose the hard threshold of the LoG filter - pixel percentage (recommended: 99.5):'))
            except ValueError:
                print('please choose a number')
                log_threshold = float(input('Choose the gaussian smoothing factor - pixel percentage (recommended: 99.5):'))
            
            log_parameters = (gaussian_param, laplace_param, log_threshold)

            image_laplace = self.log_filter(image, *log_parameters)
            image_output = image_laplace[1]
            
            
            
            plt.figure(figsize=(10,10))
            plt.imshow(image_output, cmap='gray')
            plt.show()
            
            j = 0
            while j == 0:
                decision = str(input('if the parameters are good choose "g", esle type "b":'))
                decision = decision.lower()
                
                if decision == 'g':
                    i += 1
                    j += 1
                elif decision == 'b':
                    j += 1
                else:
                    print('wrong input, please try again...')
              
        
        print('Selecting parameters for the adaptive threshold...')
        
        i = 0
        while i == 0:
            try:
                gaussian_param = float(input('Choose the gaussian smoothing factor (recommended: 1):'))
            except ValueError:
                print('please choose a number')
                gaussian_param = float(input('Choose the gaussian smoothing factor (recommended: 1):'))
            
            try:
                adaptive_block_size = int(input('Choose the block size of the adaptive filter - odd number (recommended: 5):'))
            except ValueError:
                print('please choose an odd number')
                adaptive_block_size = int(input('Choose the block size of the adaptive filter - odd number (recommended: 5):'))
            
            try:
                adaptive_offset = float(input('Choose the adaptive offset (recommended: -40):'))
            except ValueError:
                print('please choose a number')
                adaptive_offset = float(input('Choose the adaptive offset (recommended: -40):'))
            
            adaptive_threshold_parameters = (gaussian_param, adaptive_block_size, adaptive_offset)
            image_local = self.adaptive_threshold_smoothed(image, *adaptive_threshold_parameters)
            
            image_local = image_local[1]
            
            plt.figure(figsize=(10,10))
            plt.imshow(image_local, cmap='gray')
            plt.show()
            
            j = 0
            while j == 0:
                decision = str(input('if the parameters are good choose "g", esle type "b":'))
                decision = decision.lower()
                
                if decision == 'g':
                    i += 1
                    j += 1
                elif decision == 'b':
                    j += 1
                else:
                    print('wrong input, please try again...')
              
        
        print('Selecting the number of erosions...')
        i = 0
        while i == 0:
            try:
                erosion_iterations = int(input('Choose the erosion iterations (recommended: 0):'))
            except ValueError:
                print('please choose an integer number')
                erosion_iterations = int(input('Choose the erosion iterations (recommended: 0):'))
            

            if erosion_iterations > 0:
                image_local_er = scipy.ndimage.morphology.binary_erosion(image_local*image_output, iterations = erosion_iterations)
                eroded_masks = image_local_er
            elif erosion_iterations == 0:
                eroded_masks = image_local*image_output
                
            plt.figure(figsize=(10,10))
            plt.imshow(eroded_masks, cmap='gray')
            plt.show()
            
            j = 0
            while j == 0:
                decision = str(input('if the parameters are good choose "g", esle type "b":'))
                decision = decision.lower()
                
                if decision == 'g':
                    i += 1
                    j += 1
                elif decision == 'b':
                    j += 1
                else:
                    print('wrong input, please try again...')
                    
        
        
        print('Selecting the size and aspect ratio of the spots...')
        i = 0
        while i == 0:
            try:
                min_spot_size = int(input('Choose the minimum spot size in pixels (recommended: 2):'))
            except ValueError:
                print('please choose an integer number')
                min_spot_size = int(input('Choose the minimum spot size in pixels (recommended: 2):'))
                  
            try:
                max_spot_size = int(input('Choose the maximum spot size in pixels (recommended: 15):'))
            except ValueError:
                print('please choose an integer number')
                max_spot_size = int(input('Choose the maximum spot size in pixels (recommended: 15):'))
                  
            try:
                spot_aspect_ratio = float(input('Choose the choose the minimum aspect ratio of the spot (recommended: 1):'))
            except ValueError:
                print('please choose an integer number')
                spot_aspect_ratio = float(input('Choose the choose the minimum aspect ratio of the spot: (recommended: 1):'))
            
            
            
            fig, ax = plt.subplots(figsize=(40, 24))
            ax.imshow(image, cmap='viridis')
            
            microtubule_labels = label(eroded_masks)
            
            bad_regions = []
            good_regions = []
            
            for region in regionprops(microtubule_labels):
                centroid = region.centroid
                minor_axis = region.minor_axis_length
                major_axis = region.major_axis_length
                area = region.area
                
    
                # set a particle size threshold
                if area <= max_spot_size and area >= min_spot_size:
                    rect = matplotlib.patches.Rectangle((centroid[1]-major_axis/2,centroid[0]-major_axis/2), major_axis, major_axis, fill=False, edgecolor='coral', linewidth=2.5)
                    ax.add_patch(rect)
                    
                    good_regions.append(region)
               
                elif area > max_spot_size:
                    if minor_axis/major_axis < spot_aspect_ratio:
                        
                        bad_regions.append(region)
                        
                        rect = matplotlib.patches.Rectangle((centroid[1]-major_axis/2,centroid[0]-major_axis/2), major_axis, major_axis, fill=False, edgecolor='pink', linewidth=2.5)
                        ax.add_patch(rect)
            
            plt.show()
            print('The orange marked spots fall within the specified spot size and aspect ratio criteria. The pink marked spots need to be processed furhter...')
            
            j = 0
            while j == 0:
                decision = str(input('if the parameters are good choose "g", esle type "b":'))
                decision = decision.lower()
                
                if decision == 'g':
                    i += 1
                    j += 1
                elif decision == 'b':
                    j += 1
                else:
                    print('wrong input, please try again...')
            
            
            
            print('Post-processing of clustered spots...')
            
        
        print('Selecting the parameters for the Laplace of Gaussian filter to separate the clustered spots...')
        
        i = 0
        while i == 0:
            try:
                gaussian_param_zoom = float(input('Choose the gaussian smoothing factor (recommended: 1):'))
            except ValueError:
                print('please choose a number')
                gaussian_param_zoom = float(input('Choose the gaussian smoothing factor (recommended: 1):'))
            
            try:
                laplace_param_zoom = int(input('Choose the laplace filter factor (recommended: 20):'))
            except ValueError:
                print('please choose a number')
                laplace_param_zoom = int(input('Choose the laplace filter factor (recommended: 20):'))
            
            try:
                log_threshold_zoom = float(input('Choose the hard threshold of the LoG filter - pixel percentage (recommended: 94.5):'))
            except ValueError:
                print('please choose a number')
                log_threshold_zoom = float(input('Choose the gaussian smoothing factor - pixel percentage (recommended: 94.5):'))
            
            log_parameters_zoom = (gaussian_param_zoom, laplace_param_zoom, log_threshold_zoom)
            
            
            for region in bad_regions:
                centroid = region.centroid
                minor_axis = region.minor_axis_length
                major_axis = region.major_axis_length
                area = region.area
                
                
                image_cropped = image[int(centroid[0]-major_axis/2):int(centroid[0]+major_axis/2), int(centroid[1]-major_axis/2):int(centroid[1]+major_axis/2)]
                image_cropped_thres = self.log_filter(image_cropped, *log_parameters_zoom)[1]
                                
                fig, ax = plt.subplots(figsize=(5, 5))
                ax.imshow(image_cropped, cmap='viridis')
                
                sub_microtubule_labels = label(image_cropped_thres)
                    
                for sub_region in regionprops(sub_microtubule_labels):
                    sub_centroid = sub_region.centroid
                    sub_major_axis = sub_region.major_axis_length
                    sub_minor_axis = sub_region.minor_axis_length
                    sub_area = sub_region.area
                    if sub_area > 0 and sub_major_axis > 0:
                        if sub_area >= min_spot_size-2 and sub_area <= max_spot_size and sub_minor_axis/sub_major_axis < spot_aspect_ratio:
                            
                            rect = matplotlib.patches.Rectangle((sub_centroid[1]-sub_major_axis/2,sub_centroid[0]-sub_major_axis/2), sub_major_axis, sub_major_axis, fill=False, edgecolor='coral', linewidth=2.5)
                            ax.add_patch(rect)
                
                plt.show()
                print('Prerss enter to continue to the next region')
                input()
                
                
                        
            j = 0
            while j == 0:
                decision = str(input('if the post processing parameters are good choose "g", esle type "b":'))
                decision = decision.lower()
                
                if decision == 'g':
                    i += 1
                    j += 1
                elif decision == 'b':
                    j += 1
                else:
                    print('wrong input, please try again...')
        
        
        parameters_dict = {}
        parameters_dict['log'] = log_parameters
        parameters_dict['adaptive_threshold'] = adaptive_threshold_parameters
        parameters_dict['erosion'] = erosion_iterations
        parameters_dict['max_spot_area'] = max_spot_size
        parameters_dict['min_spot_area'] = min_spot_size
        parameters_dict['aspect_ratio'] = spot_aspect_ratio
        parameters_dict['log_zoom'] = log_parameters_zoom
        
        with open(self.save_path+'/'+self.experiment+'_segmentation_parameters', 'wb') as handle:
            pickle.dump(parameters_dict, handle)
        
        return parameters_dict
        

    
    def segment_microtubules(self, frame, channel, log_parameters, adaptive_threshold_parameters, erosion_iterations, min_spot_size, max_spot_size, spot_aspect_ratio, log_parameters_zoom):
        """
        This code can be used to segment microtubules in a single frame.
        
        Step 1:
            It first applies an LoG filter to smooth the image and increase the contrast. 
            The LoG image is thresholded using a hard threshold - percentile of the brightest pixels
        Step 2: 
            An adaptive threshold is applied on the LoG masked image to further separate the microtubules.
        Step 3:
            The adaptively thresholded masks are eroded to separate multiple microtubules that are segmented as one.
        Step 4:
            If a segmented microtubule is above a specified size (pixel area) or has an aspect ratio (min_axis/max_axis) below a specified threshold, 
            the algorithms zooms in and applies a stricter LoG fiter on the original pixels.
        Step 5:
            The data are organized in a pandas dataframe
        
        Input:
            frame: integer - specifies the frame in the time-lapse (statring from 0)
            log_parameters: tuple of parameters for the LoG filter (smoothing_factor, laplace_factor, percentile)
            adaptive_threshold_parameters: tuple of the parameters for the adaptive thresholding (smoothing_factor, adaptive_block_size, adaptive_offset)
            erosion_iterations: positive integer or zero - the number of erosion rounds
            spot_size: positive integer - the expected maximum size of the spots. If the detected spot is above this size threshold it is post-processed.
            aspect_ratio: positive float - the expected maximum aspect ration of the spots. If it is set to 1 all spots are post-processed since the min_axis/max_axis ratio will always be < 1
            log_parameters_zoom: tuple of parameters for the LoG filter (smoothing_factor, laplace_factor, percentile) for post-processing of clustered microtubules or other spots. 
            
        
        Returns:
            A dataframe that contains the columns:
                microtubule_df['centroid'] = microtubule_centroid 
                    WARNING: the centroid coordinates are inverted because the numpy.array images have inverted indexes (firsrt index is the y-direction, second index is the x-direction)
                microtubule_df['minor_axis'] = microtubule_minor_axis in pixels
                microtubule_df['major_axis'] = microtubule_major_axis in pixels
                microtubule_df['area'] = microtubule_area in pixels
                microtubule_df['experiment'] = experiment
                microtubule_df['frame'] = frame
            
        """
        
        image = self.image_arrays[channel][frame]
        
        # Apply an LoG filter
        image_laplace = self.log_filter(image, *log_parameters)
        
        image_output = image_laplace[1]
#        image_output_2 = image_laplace[1] * image_laplace[0]
        
        image_local = self.adaptive_threshold_smoothed(image, *adaptive_threshold_parameters)
        
        image_local = image_local[1]
        
        if erosion_iterations > 0:
            image_local_er = scipy.ndimage.morphology.binary_erosion(image_local*image_output, iterations = erosion_iterations)
            eroded_masks = image_local_er
        elif erosion_iterations == 0:
            eroded_masks = image_local*image_output
        
        
#        image_eroded = eroded_masks * image
    
#         plt.figure(figsize=(10,10))
#         plt.imshow(image_eroded, cmap='gray')
#     #     plt.clim(0,1)
#         plt.show()
        
        microtubule_labels = label(eroded_masks)
        
        microtubule_centroid = []
        microtubule_major_axis = []
        microtubule_minor_axis = []
        microtubule_area = []
    
#        i = 0
        for region in regionprops(microtubule_labels):
            
            centroid = region.centroid
            minor_axis = region.minor_axis_length
            major_axis = region.major_axis_length
            area = region.area
            
            # set a particle size threshold
            if region.area <= max_spot_size and region.area >= min_spot_size:
                
                microtubule_centroid.append(centroid)
                microtubule_minor_axis.append(minor_axis)
                microtubule_major_axis.append(major_axis)
                microtubule_area.append(area)
                 
            elif region.area > max_spot_size:
                
                if minor_axis/major_axis < spot_aspect_ratio:
                    
                    image_cropped = image[int(centroid[0]-major_axis/2):int(centroid[0]+major_axis/2), int(centroid[1]-major_axis/2):int(centroid[1]+major_axis/2)]
                    image_cropped_thres = self.log_filter(image_cropped, *log_parameters_zoom)[1]
                    
                    sub_microtubule_labels = label(image_cropped_thres)
                    
                    for sub_region in regionprops(sub_microtubule_labels):
                        
                        sub_centroid = sub_region.centroid
                        sub_major_axis = sub_region.major_axis_length
                        sub_minor_axis = sub_region.minor_axis_length
                        sub_area = sub_region.area
                        
                        if sub_area > 0 and sub_major_axis > 0:
                            if sub_area >= min_spot_size-2 and sub_area <= max_spot_size and sub_minor_axis/sub_major_axis < spot_aspect_ratio:
                                
                                sub_centroid_corrected = (sub_centroid[0]+int(centroid[0]-major_axis/2), sub_centroid[1]+int(centroid[1]-major_axis/2))
                                
                                microtubule_centroid.append(sub_centroid_corrected)
                                microtubule_minor_axis.append(sub_minor_axis)
                                microtubule_major_axis.append(sub_major_axis)
                                microtubule_area.append(sub_area)
    
                elif minor_axis/major_axis >= spot_aspect_ratio:
                    print('label with centroid:',centroid,'neglected')
                    
    
    
        microtubule_df = pd.DataFrame()
        microtubule_df['centroid'] = microtubule_centroid
        microtubule_df['minor_axis'] = microtubule_minor_axis
        microtubule_df['major_axis'] = microtubule_major_axis
        microtubule_df['area'] = microtubule_area
        microtubule_df['experiment'] = self.experiment
        microtubule_df['frame'] = frame
        
        
        # fig, ax = plt.subplots(figsize=(40, 24))
        # ax.imshow(image, cmap='viridis')
    
        
        # for index, row in microtubule_df.iterrows():
        #     rect = matplotlib.patches.Rectangle((row['centroid'][1]-row['major_axis']/2,row['centroid'][0]-row['major_axis']/2), row['major_axis'], row['major_axis'], fill=False, edgecolor='coral', linewidth=2.5)
        #     ax.add_patch(rect)
    
        #     ax.set_axis_off()
        #     plt.tight_layout()
        
        # print(self.save_path)
        
        # if os.path.isdir(self.save_path) == True:
        #     plt.savefig(self.save_path+'/'+str(frame)+'_segmentation.jpeg')
        
        # plt.show()
        
        return microtubule_df
    
    
    
    def run_segmentation(self, channel, log_parameters, adaptive_threshold_parameters, erosion_iterations, min_spot_size, max_spot_size, spot_aspect_ratio, log_parameters_zoom):
        """
        This function is used to run the segment_microtubules function on all the frames of the time-lapse movie.
        It iterates over the frames and the image arrays retrieved from the nd2 metadata and the specified channel and applies the segmentation function on each frame.
        
        Input:
            log_parameters: tuple of parameters for the LoG filter (smoothing_factor, laplace_factor, percentile)
            adaptive_threshold_parameters: tuple of the parameters for the adaptive thresholding (smoothing_factor, adaptive_block_size, adaptive_offset)
            erosion_iterations: positive integer or zero - the number of erosion rounds
            spot_size: positive integer - the expected maximum size of the spots. If the detected spot is above this size threshold it is post-processed.
            aspect_ratio: positive float - the expected maximum aspect ration of the spots. If it is set to 1 all spots are post-processed since the min_axis/max_axis ratio will always be < 1
            log_parameters_zoom: tuple of parameters for the LoG filter (smoothing_factor, laplace_factor, percentile) for post-processing of clustered microtubules or other spots. 
        Returns:
            A dataframe that contains the columns:
                microtubule_df['centroid'] = microtubule_centroid (center of mass)
                    WARNING: the centroid coordinates are inverted because the numpy.array images have inverted indexes (firsrt index is the y-direction, second index is the x-direction)
                microtubule_df['minor_axis'] = microtubule_minor_axis in pixels
                microtubule_df['major_axis'] = microtubule_major_axis in pixels
                microtubule_df['area'] = microtubule_area in pixels
                microtubule_df['experiment'] = experiment
                microtubule_df['frame'] = frame
                microtubule_df['x'] = the x coordinate of the spot center of mass (centroid)
                microtubule_df['y'] = the y coordinate of the spot center of mass (centroid)
                microtubule_df['mean_fluorescence'] = the mean spot fluroescence estimated in a 5x5 pixel array around the center of mass (centroid)
                microtubule_df['median_fluorescence'] = the median spot fluroescence estimated in a 5x5 pixel array around the center of mass (centroid)
                microtubule_df['std_fluorescence'] = the standard deviation of the spot fluroescence estimated in a 5x5 pixel array around the center of mass (centroid)
                microtubule_df['time_sec'] = the actual time that corresponds to each spot
           This dataframe is saved in the specified destination (save_path)
        
        """
        
        for frame in self.frames:
            
            dataframe = self.segment_microtubules(frame, channel, log_parameters, adaptive_threshold_parameters, erosion_iterations, min_spot_size, max_spot_size, spot_aspect_ratio, log_parameters_zoom)
            
            if frame == 0:
                
                microtubule_df = dataframe
            
            elif frame > 0:
                
                microtubule_df = pd.concat([microtubule_df, dataframe])
        
        # Applying an 'x' and a 'y' column in the dataframe which includes the centroid coordinates
        microtubule_df['x'] = microtubule_df['centroid'].apply(lambda x: x[1])
        microtubule_df['y'] = microtubule_df['centroid'].apply(lambda x: x[0])
        
        
        
        # Adding a column for the mean fluorescence of the spot, the median fluorescence and the standard deviation of the fluorescence.
        # All these are measured in an array of 5x5 pixels around the centroid of the spot.
        
        mean_fluorescence = []
        median_fluorescence = []
        std_fluorescence = []
        
        for index, row in microtubule_df.iterrows():
            frame = row['frame']
            x = row['x']
            y = row['y']
            x_int = int(x)
            y_int = int(y)
            
            image = self.image_arrays['TEST 488'][frame].copy()
            
            image_cropped = image[(y_int-2):(y_int+3), (x_int-2):(x_int+3)]
            pixels = image_cropped.ravel()
        #     print(pixels)
        #     plt.imshow(image_cropped)
        #     plt.clim(200,400)
        #     plt.show()
            mean_fluorescence.append(np.mean(pixels))
            median_fluorescence.append(np.median(pixels))
            std_fluorescence.append(np.std(pixels))
        
        microtubule_df['mean_fluorescence'] = mean_fluorescence
        microtubule_df['median_fluorescence'] = median_fluorescence
        microtubule_df['std_fluorescence'] = std_fluorescence
        microtubule_df['time_sec'] = microtubule_df['frame'] * self.interval
        
        
        # reseting the index of the dataframe
        microtubule_df = microtubule_df.reset_index()
        
        
                
        with open(self.save_path+'/'+self.experiment+'_segmented_microtubules_df', 'wb') as handle:
            microtubule_df.to_pickle(path = handle, compression='zip', protocol = pickle.HIGHEST_PROTOCOL)
        
        return microtubule_df
    
    
    
    def show_segmentation(self, channel, frames, save):
        """
        Shows the segmentation of the particles across the specified frame range and saves the images to make a movie
        
        Input:
            channel: string - the channel where the microtubule comets are displayed (e.g. 'TEST 488')
            frames: tuple of integers - the first and the last frame that will be displayed (e.g. (0,50) will display frames 0 to 49)
            save: stirng - the path of the filter where the images will be saved
        Returns:
            Displays the segmented commets on the respective fluorescence image
            If the save string is a valid directory, the images are saved in the specified directory
                If not, error is returned.
        """
        
        
        if os.path.isdir(save) == False:
            print('The specified save path is not an valid directory. The images will not be saved...')
        
        microtubule_df = pd.read_pickle(self.save_path+'/'+self.experiment+'_segmented_microtubules_df', compression='zip')
        
        for fr in range(frames[0], frames[1]):
            
            frame_df = microtubule_df[microtubule_df['frame']==fr]
            
        
            image = self.image_arrays[channel][fr]
            
            fig, ax = plt.subplots(figsize=(40, 24))
            ax.imshow(image, cmap='viridis')
        
            
            for index, row in frame_df.iterrows():
                rect = matplotlib.patches.Rectangle((row['centroid'][1]-row['major_axis']/2,row['centroid'][0]-row['major_axis']/2), row['major_axis'], row['major_axis'], fill=False, edgecolor='coral', linewidth=2.5)
                ax.add_patch(rect)
        
                ax.set_axis_off()
                plt.tight_layout()
            
            if os.path.isdir(save) == True:
                plt.savefig(save+'/'+str(fr)+'_segmentation.jpeg')
                
            plt.show()
    

 
    def microtubule_tracking(self, max_radius, memory):
        """
        This code is used to detect the closest detected spot in the next frame for every spot in the previous frame.
        It measures the distances with all adjacent spots in the next frame, within a circle with a radius = max_radius, and picks the one which is the closest (euclidean distance is considered).
        If a position within the max radius is not linked to the comet in the next frame, the later time points will be used within a specified memory of frames. 
        
        Input:
            max_radius: integer or float - the maximum radius (pixels) whithin which the algorithm is searching for proximal microtubules in the next frame. 
                        Otherwise, this is the maximum distance that a microtubule can cover within a frame.
            memory: integer - how many frames ahead will the algorithm look for linked comet positions
        Returns:
            [0] The original segmentation dataframe is returned with two columns added. 
                The 'linked_index' column includes the index of each linked comet to the current comet.
                The 'microtubule_id' column includes a unique ID for each set of linked comets or microtubule/
                This dataframe is also save in the specified path (self.save_path)
            [1] A dictionary which links the microtubule indexes (list of integers - value) to the trajectory id (string - key)
        """
        
    #     ###
    #     Assumes that the microtubules do not split
    #     ###
        

        dataframe_microtubules = pd.read_pickle(self.save_path+'/'+self.experiment+'_segmented_microtubules_df', compression='zip')

        linkage_dictionary = {}
        trajectories_dictionary = {}
        starting_id_dictionary = {}
        
        for index, row in dataframe_microtubules.iterrows():
            
            if index not in linkage_dictionary.keys():
                frame = row['frame']
                current_x = row['x']
                current_y = row['y']
                
                next_frame = frame + 1
                
                if index not in linkage_dictionary.values():
                    starting_id = 'micro_'+str(index)
                    trajectories_dictionary[starting_id] = [index]
                    starting_id_dictionary[index] = starting_id
                else:
                    starting_id = starting_id_dictionary[index]
    
                
                while next_frame <= frame + memory and next_frame < self.n_frames:
                    
                    next_dataframe = dataframe_microtubules[dataframe_microtubules['frame']==next_frame]
                    next_dataframe['x_diff'] = next_dataframe['x'] - current_x
                    next_dataframe['y_diff'] = next_dataframe['y'] - current_y
                    next_dataframe['distance'] = np.sqrt(next_dataframe['x_diff']**2 + next_dataframe['y_diff']**2)
                    
                    minimum_distance_df = next_dataframe[next_dataframe['distance']==next_dataframe['distance'].min()]
                    
                    if minimum_distance_df['distance'].values[0] < max_radius:
                        
                        linked_index = minimum_distance_df.index.values[0]
                        linkage_dictionary[index] = linked_index
                        starting_id_dictionary[linked_index] = starting_id
                        trajectories_dictionary[starting_id].append(linked_index)
                        
                        print('microtubule index:', index, ', in frame:', frame, ', linked to index:', linked_index, ', in frame:', next_frame)
                        next_frame = frame + memory + 1
                        
                    elif minimum_distance_df['distance'].values[0] >= max_radius:
                        
                        next_frame += 1
        
        
        dataframe_microtubules['linked_index'] = dataframe_microtubules.index.map(linkage_dictionary)
        dataframe_microtubules['microtubule_id'] = dataframe_microtubules.index.map(starting_id_dictionary)
        
        with open(self.save_path+'/'+self.experiment+'_tracked_microtubules_df', 'wb') as handle:
            dataframe_microtubules.to_pickle(path = handle, compression='zip', protocol = pickle.HIGHEST_PROTOCOL)
        
        
        return dataframe_microtubules, trajectories_dictionary
    
    
    
    def show_microtubule_trajectories(self, threshold):
        """
        Shows the microtubule trajectories on the first image of the time-lapse or stream acquisition.
        The points that belong to the same trajectory are linked with a white line.
        The plasma map is used to show the frame of each comet position.
        
        Input:
            threshold: integer - The minimum length (frames) of the trajectories shown.
        """
        
        microtubule_df = pd.read_pickle(self.save_path+'/'+self.experiment+'_tracked_microtubules_df', compression='zip')
        
        plt.figure(figsize=(20,20))
        plt.imshow(self.image_arrays['TEST 488'][0].copy())  
        for traj in microtubule_df['microtubule_id'].unique():
            
            trajectory_df = microtubule_df[microtubule_df['microtubule_id']==traj]
            
            if trajectory_df.shape[0] >= 15:
            
                x_coords = trajectory_df['x']
                y_coords = trajectory_df['y']
                frames = trajectory_df['frame']
            
                plt.plot(x_coords, y_coords, color='white', linewidth=0.5)
                plt.scatter(x_coords, y_coords, c=frames, cmap='plasma', vmin=0, vmax=frames.max())
                #     plt.scatter([x_coords[0], x_coords[-1]], [y_coords[0], y_coords[-1]], c=[frames[0], frames[-1]], cmap='plasma', vmin=0, vmax=dataframe_microtubules['frame'].max())
                
                #plt.clim(200,400)
        plt.colorbar()
        plt.show()



    def get_microtubule_speed_and_angles(self, time_lag):
        """
        This algorithm is used to estimate the angle and speed of the microtubule trajectories.
        It uses the class functions to get the central line and the gut mask in the embryo.
        The microtubule trajectory angles are corrected for the anlge of the central line.
        The minimum distance of the comets from the central line is also estimatied.
        
        Input:
            time_lag: integer - a smoothing window used for the angle and speed estimation.
                This is used to improve the statistics and minimize the contribution 
                of the spot center detection error between subsequent frames 
                to the trajectory speed and angle estimation.
            WARNING:
                The user needs to run the get_central_line, get_gut_mask and
                                          get_embryo mask functions before running this one
        Returns:
            An updated microtubule dataframe with the following columns added:
                dataframe['x_shifted'] # the shifted x positions by the time_lag
                dataframe['y_shifted'] # the shifted y positions by the time_lag
                dataframe['frame_shifted'] # the shifted frame positions by the time_lag
                dataframe['dx'] # horizontal displacement
                dataframe['dy'] # vertical displacement
                dataframe['dframe'] # frame displacement
                dataframe['slope'] # slope (dy/dx)
                dataframe['angle'] # angle of the displacement
                dataframe['central_line_distance_slope'] # a tuple with the minimum distance from the central line and the slope at that position of the central line, as well as the relative coordinates of the central line
                dataframe['central_distance'] # the minimum distance of the comet from the central line in pixels
                dataframe['central_slope'] # the slope of the central line at the point with the minimum distance from the comet
                dataframe['central_x'] # the the relative central line x coordinates
                dataframe['central_y'] # the the relative central line y coordinates
                dataframe['1d_position'] # the 1D position on the central line arch length
                dataframe['central_angle'] # the angle of the central line at the point with the minimum distance from the comet (estimated from the slope)
                dataframe['angle_corrected'] # angle corrected to the central line angle
                dataframe['angle_corrected_360'] # angles above 360 degrees are corrected by subtracting 360
                dataframe['distance_um'] # 2D displacement (sqrt(dx^2 + dy^2))
                dataframe['speed'] # the instantaneous speed of the displacement in μm/sec
                dataframe['central_distance_um'] # the minimum distance of the comet position from the central line in μm
                dataframe['gut_position'] # the position in the gut (True) or outside of the embryo gut (False)
                dataframe['embryo_position'] # the position in the embryo of interest (True) or inside some other embryo in the field of view (False)
            
            This dataframe is saved in the specified directory as "angle_speed_df"
        
        """
        
        microtubule_df = pd.read_pickle(self.save_path+'/'+self.experiment+'_tracked_microtubules_df', compression='zip')
        

        with open(self.save_path+'/'+self.experiment+'_central_line_fitted_coordinates', 'rb') as handle:
            central_line_coordinates = pickle.load(handle)
        
        with open(self.save_path+'/'+self.experiment+'_gut_mask', 'rb') as handle:
            gut_mask = pickle.load(handle)
        
        with open(self.save_path+'/'+self.experiment+'_embryo_mask', 'rb') as handle:
            embryo_mask = pickle.load(handle)
        
        self.central_line_coords = central_line_coordinates
        self.gut_mask = gut_mask
        self.embryo_mask = embryo_mask
            
            
            
#        central_line = self.get_central_line(gut_channel, adaptive_threshold_parameters, depth_threshold)
#        
#        gut_mask = self.get_gut_mask(gut_channel, adaptive_threshold_parameters, close_iter, dilation, object_threshold)
#        
#        embryo_mask = self.get_embryo_mask(gut_channel, adaptive_threshold_parameters, close_iter, dilation, object_threshold)
        
        
        def get_angle_from_slope(dx, dy):
            """
            This function is used to estimate the angle of each microtubule displacement from two adjacent time-points.
            The angle is estmated using the slope of the microtubule displacement and the orientation (negative or positive dx and dy)
            
            Specifically, the slope is converted to an angle.
            Using the dx and dy displacements it is decided to which of the quartiles in a circle the angle belongs and the appropriate adjustments are made.
            
            Input:
                displacements: dx, dy - The displacement between two linked spots in the x and y direction
            Returns:
                angle: float - the angle of the displacement in degrees
            """
            
            if dx != 0:
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
        
        def get_distance_from_line(line_coords, point_coords):
            """
            Estimates the minimum distance of a point from the line and the slope of the line at the point with the minimum distance.
            
            Input:
                line_coords: tuple of two numpy arrays, one for the x and one for the y coordinates
                point_coords: tuple of two floats, the first corresponding to the x_coordinate, and the second to the y_coordinate of the comet
            
            Returns:
                [0] the minimum distance of a comet from the central line with a precision of 0.1 pixels
                [1] the slope of the central line at the point with the minimum distance from the comet
            """
            
            dataframe = pd.DataFrame()
            
            dataframe['line_x'] = line_coords[0]
            dataframe['line_y'] = line_coords[1]
            dataframe['dx_sqr'] = (dataframe['line_x'] - point_coords[0])**2
            dataframe['dy_sqr'] = (dataframe['line_y'] - point_coords[1])**2
            dataframe['distance'] = np.sqrt(dataframe['dx_sqr']+dataframe['dy_sqr'])
            
            # get the arch length of the central line
            line_dx = dataframe['line_x'] - (dataframe['line_x']).shift(1)
            line_dy = dataframe['line_y'] - (dataframe['line_y']).shift(1)
            line_sqr = np.sqrt(line_dx**2 + line_dy**2)
            line_sqr[0] = 0
            arch_length = np.cumsum(line_sqr)
            
            dataframe['arch_length'] = arch_length

            
            min_dataframe = dataframe[dataframe['distance']==dataframe['distance'].min()]
            min_index = min_dataframe.index[0]
            min_x = min_dataframe.line_x.values[0]
            min_y = min_dataframe.line_y.values[0]
            min_arch_length = min_dataframe.arch_length.values[0]
            central_line_distance = dataframe['distance'].min()
            
            if min_index >= 10 and min_index < dataframe.shape[0]-10:
                slope_dataframe = dataframe[(min_index-10):(min_index+11)]
            elif min_index < 10 and min_index < dataframe.shape[0]-10:
                slope_dataframe = dataframe[0:(min_index+11)]
            elif min_index >= 10 and min_index >= dataframe.shape[0]-10:
                slope_dataframe = dataframe[(min_index-10):]
            
            central_line_slope = np.polyfit(slope_dataframe['line_x'], slope_dataframe['line_y'], 1)[0]
            
            
            return (central_line_distance, central_line_slope, min_x, min_y, min_arch_length)
                
        def is_point_in_mask(mask, x, y):
            """
            Checks if a comet is in the masked embryo gut.
            """
            return mask[int(y)][int(x)]
            
        def large_angle_correction(angle):
            """
            Those angles that are larger than 360 have 360 degrees subtracted.
            For instance, and angle of 395o will be adjusted to 35o
            """
            if angle >= 360:
                
                return angle - 360
            
            elif angle < 360:
                
                return angle
        
        
#        time_lag = 4
    
        i = 0
#        plt.figure(figsize=(10,10))
        for traj in microtubule_df['microtubule_id'].unique():
            print('analyzing trajectory:',traj)
            trajectory_df = microtubule_df[microtubule_df['microtubule_id']==traj]
            
            trajectory_df['x_shifted'] = trajectory_df['x'].shift(-time_lag)
            trajectory_df['y_shifted'] = trajectory_df['y'].shift(-time_lag)
            trajectory_df['frame_shifted'] = trajectory_df['frame'].shift(-time_lag)
            trajectory_df['dx'] = trajectory_df['x_shifted'] - trajectory_df['x']
            trajectory_df['dy'] = trajectory_df['y_shifted'] - trajectory_df['y']
            trajectory_df['dframe'] = trajectory_df['frame_shifted'] - trajectory_df['frame']
            trajectory_df['slope'] = trajectory_df['dy']/trajectory_df['dx']
            trajectory_df['angle'] = trajectory_df.apply(lambda row: get_angle_from_slope(row['dx'], row['dy']), axis=1)
            trajectory_df['central_line_distance_slope'] = trajectory_df.apply(lambda row: get_distance_from_line(self.central_line_coords, (row['x'], row['y'])), axis=1)
            trajectory_df[['central_distance', 'central_slope', 'central_x', 'central_y', 'd_position']] = pd.DataFrame(trajectory_df['central_line_distance_slope'].tolist(), index=trajectory_df.index)               
            trajectory_df['central_angle'] = np.rad2deg(np.arctan(trajectory_df['central_slope']))
            trajectory_df['angle_corrected'] = trajectory_df['angle'] - trajectory_df['central_angle']
            trajectory_df['angle_corrected_360'] = trajectory_df.apply(lambda row: large_angle_correction(row['angle_corrected']), axis=1)
            trajectory_df['distance_um'] =  np.sqrt(trajectory_df['dx']**2 + trajectory_df['dy']**2) * self.scale
            trajectory_df['speed'] = trajectory_df['distance_um']/(trajectory_df['dframe']*self.interval)
            trajectory_df['central_distance_um'] = trajectory_df['central_distance'] * self.scale
            trajectory_df['gut_position'] = trajectory_df.apply(lambda row: is_point_in_mask(self.gut_mask ,row['x'],row['y']), axis=1)
            trajectory_df['embryo_position'] = trajectory_df.apply(lambda row: is_point_in_mask(self.embryo_mask ,row['x'],row['y']), axis=1)
            
            if i == 0:
                angle_df = trajectory_df
                i += 1
            elif i > 0:
                angle_df = pd.concat([angle_df, trajectory_df])
        
        
        embryo_df = angle_df[angle_df['embryo_position']==True]
        gut_df = embryo_df[embryo_df['gut_position']==True]
        non_gut_df = embryo_df[embryo_df['gut_position']==False]
        
        plt.figure(figsize=(6,6))
        plt.plot(gut_df['x'], gut_df['y'], 'o', color='red', label='embryo gut', markersize=2)
        plt.plot(non_gut_df['x'], non_gut_df['y'], 'o', color='blue', label='embryo', markersize=2)
        plt.xlabel('x coordinates (px)')
        plt.ylabel('y coordinates (px)')
        plt.legend()
        
        with open(self.save_path+'/'+self.experiment+'_angle_speed_df', 'wb') as handle:
            angle_df.to_pickle(path = handle, compression='zip', protocol = pickle.HIGHEST_PROTOCOL)
        
            
        return angle_df

    
    
    def max_projection(self):
        """
        Returns the max projecton of the stream acquisition images.
        """
        
        stacked = np.empty((self.n_frames, self.metadata['height'], self.metadata['width']))
        for n in range(self.n_frames):
            stacked[n] = self.image_arrays[self.channels[0]][n]
        
        im_max= np.max(stacked, axis=0)
        
        return im_max
    
    
    
    def embryo_and_gut_filter(self, image, adaptive_threshold_parameters, close_iter, dilation, object_threshold):
        """
        This an adaptive filter applied to segment the entire embruo and its gut. 
        It does not only apply an adaptive filter (adaptive_threshold_smoothed filter), but it also bridges the gaps in the masks,
        fills the holes, removes small objects and dilates if selected.
        
        Input:
            image: numpy array - the image where the filter is applied.
            adaptive_threshold_parameters: tuple of integers corresponding to the parameters of the adaptive_threshold_smoothed filter
                (smoothing_factor, adaptive_block_size, adaptive_offset)
            close_iter: integer - the iterations of the scipy.ndimage.binary_closing process
            object_threshold: integer - the minimum number of pixels of the gut mask (objects with lower areas are removed)
                skimage.morphology.remove_small_objects
                recommended #5000
            dilation: integer - the number of dilation iterations in the scipy.ndimage.binary_dilation process
        
        Returns:
            The masked image - binary numpy array 
            
        """
        
        # the adaptive threshold parameters
        smoothing_factor, adaptive_block_size, adaptive_offset = adaptive_threshold_parameters
        
        # Threshold the image
        # 2, 71, -20
        adaptive_image = self.adaptive_threshold_smoothed(image, smoothing_factor, adaptive_block_size, adaptive_offset)
        
        # Bridge the gaps
        adaptive_image_closed = scipy.ndimage.binary_closing(adaptive_image[1], iterations=close_iter)
        #10
        
        # Fill the holes
        adaptive_image_closed_full = scipy.ndimage.morphology.binary_fill_holes(adaptive_image_closed)
        
        # Remove small objects
        adaptive_image_closed_full_dilated = morphology.remove_small_objects(adaptive_image_closed_full, min_size=object_threshold, connectivity=1, in_place=False)
        #100
        
        # Dilate the mask
        if dilation > 0:
            adaptive_image_closed_full_dilated = scipy.ndimage.binary_dilation(adaptive_image_closed_full_dilated, iterations=dilation)
            #1
        
        return adaptive_image_closed_full_dilated
        
    
        
    def get_embryo_mask(self, embryo_channel, adaptive_threshold_parameters, close_iter, dilation, object_threshold):
        """
        This algorithm is used to segment the embryo based on the fluorescence channel where the microtubule comets were imaged.
        First the micrrotubule channel image is thresholded using an adaptive threshold
        Then the gaps in the gut mask are bridged.
        Then the holes in the mask are filled.
        Then the mask is dilated.
        Then the small patches around the mask are removed.
        
        Input:
            gut_channel: string - the channel where the gut marker was imaged (e.g. 'TEST 405')
            adaptive_threshold_parameters: tuple of integers corresponding to the parameters of the adaptive_threshold_smoothed filter
                (smoothing_factor, adaptive_block_size, adaptive_offset)
                recommended # 2, 501, -1
            close_iter: integer - the iterations of the scipy.ndimage.binary_closing process
                recommended #8-15
            object_threshold: integer - the minimum number of pixels of the gut mask (objects with lower areas are removed)
                skimage.morphology.remove_small_objects
                recommended #5000
             dilation: integer - the number of dilation iterations in the scipy.ndimage.binary_dilation process
                 recommended #0
        Returns
            final_gut_mask
                The masked gut image - numpy array
                
                The same mask is stored in the class as a global variable 
                self.embryo_mask
                The masked image is also saved in the specified folder (self.save_path)
        """
        
        snapshot_images = self.nd2_snapshots_to_array(self.snapshots_path)
        
        try:
            image = snapshot_images[1][embryo_channel]
        except KeyError:
            print('The specified channel could not be located. The max projection from the stream acquisition will be used to segment the embryo.')
            image = self.max_projection()

        masked_image = self.embryo_and_gut_filter(image, adaptive_threshold_parameters, close_iter, dilation, object_threshold)
        
        embryo_labels = label(masked_image)
        
        i = 0
        size = 0
        for lbl in range(1, embryo_labels.max()+1):
            embryo_mask = embryo_labels==lbl

            if np.count_nonzero(embryo_mask) > object_threshold:
                size = np.count_nonzero(embryo_mask)
                
                if i == 0:
                    final_embryo_mask = embryo_mask
                    i += 1
                elif i > 0:
                    if size > np.count_nonzero(final_embryo_mask):
                        final_embryo_mask = embryo_mask
        
        try:
            final_embryo_mask = scipy.ndimage.binary_closing(final_embryo_mask, iterations=close_iter*2)
            final_embryo_mask = scipy.ndimage.morphology.binary_fill_holes(final_embryo_mask)
            
            plt.imshow(final_embryo_mask)
            plt.show()
            
            self.embryo_mask = final_embryo_mask
            
            with open(self.save_path+'/'+self.experiment+'_embryo_mask', 'wb') as handle:
                pickle.dump(final_embryo_mask, handle)
            
            return final_embryo_mask
        
        except UnboundLocalError:
            print('This set of parameters does not yield any mask.')
            
        
    
    
    def get_gut_mask(self, gut_channel, adaptive_threshold_parameters, close_iter, dilation, object_threshold):
        """
        This algorithm is used to segment the embryo gut based on the fluorescence channel where the gut marker was imaged?
        First the gut channel image is thresholded using an adaptive threshold
        Then the gaps in the gut mask are bridged.
        Then the holes in the mask are filled.
        Then the small patches around the mask are removed.
        Then the mask is dilated.
            The processes so far are part of the embryo_and_gut_filter function
        Then it selects the mask that is above the specified threshold and it does not span more than 250 pixels
        in both the x and y dimensions. This is useful to separate the borderline of the entire embruo from the borders of the gut.
            
        Input:
            gut_channel: string - the channel where the gut marker was imaged (e.g. 'TEST 405')
            
            For the embryo_and_gut_filter function:
                adaptive_threshold_parameters: tuple of integers corresponding to the parameters of the adaptive_threshold_smoothed filter
                    (smoothing_factor, adaptive_block_size, adaptive_offset)
                    recommended # 2, 71, -20
                close_iter: integer - the iterations of the scipy.ndimage.binary_closing process
                    recommended #9
                object_threshold: integer - the minimum number of pixels of the gut mask (objects with lower areas are removed)
                    skimage.morphology.remove_small_objects
                    recommended #5000
                 dilation: integer - the number of dilation iterations in the scipy.ndimage.binary_dilation process
                     recommended #0
        Returns
            final_gut_mask
                The masked gut image - numpy array
            
                The same mask is stored in the class as a global variable 
                self.gut_mask
                The masked image is also saved in the specified folder (self.save_path)
        """
        
        snapshot_images = self.nd2_snapshots_to_array(self.snapshots_path)
        
        try:
            image = snapshot_images[1][gut_channel]
        except KeyError:
            print('The specified channel could not be located. The max projection from the stream acquisition will be used to segment the gut.')
            image = self.max_projection()
        
        masked_image = self.embryo_and_gut_filter(image, adaptive_threshold_parameters, close_iter, dilation, object_threshold)
            
        gut_labels = label(masked_image)
        
        for lbl in range(1, gut_labels.max()+1):
            gut_mask = gut_labels==lbl
            if np.count_nonzero(gut_mask) > object_threshold:
                print(np.count_nonzero(gut_mask))
                gut_coordinates = np.argwhere(gut_mask == True)
                y_gut, x_gut = map(list, zip(*gut_coordinates))
            
                if max(y_gut)-min(y_gut) < 250 and max(x_gut)-min(x_gut) < 250:
                    plt.imshow(gut_mask)
                    plt.show()
                    final_gut_mask = gut_mask
        
        try:
            self.gut_mask = final_gut_mask
            
            with open(self.save_path+'/'+self.experiment+'_gut_mask', 'wb') as handle:
                pickle.dump(final_gut_mask, handle)
                
            return final_gut_mask
            
        except UnboundLocalError:
            print('This set of parameters does not yield any mask.')
            
        
        
    
    def get_central_line(self, channel, order=1, smoothing=100, lut = (0,2500), rotate=False):
        """
        This code finds the central line on the basis of the gut marker,
        fits a first degree polynomial to the centrai lline pixels,
        and estimates the slope to correct the angles of the microtubule trajectories. 
        
        Input:
            channel: string - the name of the channel that will be used to degine the central line
            rotate: select True if you want to draw the central line on the rotate image (clockwise 90o rotation)
            order: the order of the fitted smoothing spline (integer larger than 0)
            lut: tuple of integers - the min and the max pixel values of the LUTs used to set the contrast and brightness of the image.
            smoothing: the smoothing factor of the fitted smoothing spline (integer with 0 corresponding to interpolation - no smoothing)
        Returns:
            [0] The (x coordinares, y coordinares) of the central line - tuple of 2 numpy arrays
            The same variables are stored in the class as global variables:
                self.central_line_coords = (final_x, final_y)
            The central line manually selected and fitted coordinates are also saved in the specified folder (self.save_path)
        """
        
        snapshot_images = self.nd2_snapshots_to_array(self.snapshots_path)
        
        try:
            image = snapshot_images[1][channel]
        except KeyError:
            print('The specified channel could not be located.  The max projection from the stream acquisition will be used to draw the central line.')
            image = self.max_projection()
        
        if rotate == True:
            plt.imshow(np.rot90(image), cmap='gray', vmin=lut[0], vmax=lut[1])
            coordinates = plt.ginput(n=0, timeout=0)
            
            x_coords, y_coords = map(list, zip(*coordinates))
            x_coords = np.array(x_coords)
            y_coords = np.array(y_coords)
            
            unique_x, index_x = np.unique(x_coords, axis=0, return_index=True)
            unique_y = y_coords[index_x]
            
            permutation = unique_x.argsort()
            x_coords_sorted = unique_x[permutation]
            y_coords_sorted = unique_y[permutation]
            
            spline = UnivariateSpline(x_coords_sorted, y_coords_sorted, k=order, s=smoothing)
            
            smoothed_x = np.arange(min(x_coords_sorted),max(x_coords_sorted), 0.1)
            smoothed_y = spline(smoothed_x)
            
            final_x = -smoothed_y+self.metadata['height']
#            final_x = -smoothed_y+512
            final_y = smoothed_x
            
            
        elif rotate == False:
            plt.imshow(image, cmap='gray', vmin=lut[0], vmax=lut[1])
            coordinates = plt.ginput(n=0, timeout=0)
            
            x_coords, y_coords = map(list, zip(*coordinates))
            x_coords = np.array(x_coords)
            y_coords = np.array(y_coords)
            
            unique_x, index_x = np.unique(x_coords, axis=0, return_index=True)
            unique_y = y_coords[index_x]
            
            permutation = unique_x.argsort()
            x_coords_sorted = unique_x[permutation]
            y_coords_sorted = unique_y[permutation]
            
            spline = UnivariateSpline(x_coords_sorted, y_coords_sorted, k=order, s=smoothing)
            
            smoothed_x = np.arange(min(x_coords_sorted),max(x_coords_sorted), 0.1)
            smoothed_y = spline(smoothed_x)
        
            final_x = smoothed_x
            final_y = smoothed_y
        
        plt.imshow(image, cmap='gray', vmin=lut[0], vmax=lut[1])
        plt.plot(final_x, final_y, color='red')
        plt.show()
        
        central_line_coords = (final_x, final_y)
        self.central_line_coords = central_line_coords
        
        
        
        with open(self.save_path+'/'+self.experiment+'_central_line_fitted_coordinates', 'wb') as handle:
            pickle.dump(central_line_coords, handle)
        
        with open(self.save_path+'/'+self.experiment+'_central_line_manual_coordinates', 'wb') as handle:
            pickle.dump(coordinates, handle)


        return central_line_coords
