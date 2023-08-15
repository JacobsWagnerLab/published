#Python 3.7.0 (default, Jun 28 2018, 07:39:16)
#Type "copyright", "credits" or "license" for more information.
#
#IPython 7.8.0 -- An enhanced Interactive Python.

"""
Created on Mon Aug  5 14:32:50 2019

@author: Alexandros Papagiannakis, Christine Jacobs-Wagner lab, Sarafan ChEM-H, Stanford University 2019
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from skimage import io, filters
from skimage.measure import label, regionprops
import scipy
from skimage.filters import threshold_local
from scipy import ndimage
import pickle
import os

class particle_tracking_microfluidics(object):
    
    def __init__(self, images_path, save_path, experiment, position, channel, edge_width):
        """
        Developer: Alexandros Papagiannakis, Stanford University, May 2021
        
        This class contains all the functions used for particle tracking in the mother machine
        """
        
        def get_fluorescence_image(images_path, xy_position, channel, position, frame, edge_width):
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
        
        def get_channel_fluorescence_images(images_path, xy_position, channel, position, frame_range, edge_width):
            fluorescence_images_dict = {}
            for frame in range(frame_range[0], frame_range[1]):
                fluorescence_images_dict[frame]=get_fluorescence_image(images_path, xy_position, channel, position, frame, edge_width)
            return fluorescence_images_dict

        def background_correction(fluor_image, edge_width):
            """
            Used in the get_fluorescence images function
            """
            return fluor_image - np.mean(np.concatenate((fluor_image[:, 0:edge_width],fluor_image[:, -edge_width:]),axis=1), axis=1)[:,None]
        
        with open(save_path + '/' +experiment+'_pos'+str(position)+ '_curated_watershed_dict', 'rb') as handle:
            curated_watershed_dict = pickle.load(handle)
        self.cell_labels = curated_watershed_dict
        frame_range = (min(self.cell_labels.keys()), max(self.cell_labels.keys())+1)
        
        xy_position = experiment[experiment.find('xy'):]
        self.xy_position = xy_position
        self.particle_images = get_channel_fluorescence_images(images_path, xy_position, channel, position, frame_range, edge_width) # background corrected images
        self.save_path = save_path # the watershed labels are saved and where the particle positions will be saved
        self.experiment = experiment # the experiment ID including the xy_position
        self.xy_position = xy_position # the XY -position (string)
        self.position = position # the channel position (integer)
        self.position_string = 'pos'+str(position)
        self.sensor = self.particle_images[0].shape
        
        
    #--------- USED FUNCTIONS ----------#
    # The 2D Gaussian with rotation fit functions are used to the estimate 
    # the particle center with subpixel resolution
    def gaussian(self, height, center_x, center_y, width_x, width_y, rotation):
        """
        Returns a gaussian function with the given parameters
        
        Reference
        ---------
        Originally implemented by Andrew Giessel
        Code available on GitHub:
            https://gist.github.com/andrewgiessel/6122739
        
        Notes
        -----
        Modified by Alexandros Papagiannakis to correct the rotation of the 2D Gaussian around the central pixel
        """
        width_x = float(width_x)
        width_y = float(width_y)
        
        rotation = np.deg2rad(rotation)
    #    center_x = center_x * np.cos(rotation) - center_y * np.sin(rotation)
    #    center_y = center_x * np.sin(rotation) + center_y * np.cos(rotation)
        
        def rotgauss(x,y):
            x = x-center_x # modification
            y = y-center_y # modification
            xp = x * np.cos(rotation) - y * np.sin(rotation)
            yp = x * np.sin(rotation) + y * np.cos(rotation)
            xp = xp + center_x
            yp = yp + center_y
            g = height*np.exp(
                -(((center_x-xp)/width_x)**2+
                  ((center_y-yp)/width_y)**2)/2.)
            return g
        return rotgauss
    
    def fitgaussian(self, data):
        """Returns (height, x, y, width_x, width_y)
        the gaussian parameters of a 2D distribution found by a fit"""
        
        def moments(data):
            """Returns (height, x, y, width_x, width_y)
            the gaussian parameters of a 2D distribution by calculating its
            moments """
            total = data.sum()
            X, Y = np.indices(data.shape)
            x = (X*data).sum()/total
            y = (Y*data).sum()/total

            col = data[:, int(y)]

            width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())

            row = data[int(x), :]

            width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
            height = data.max()
            return height, x, y, width_x, width_y, 0.0
        
        params = moments(data)
        errorfunction = lambda p: np.ravel(self.gaussian(*p)(*np.indices(data.shape)) - data)
        p, success = scipy.optimize.leastsq(errorfunction, params)
        return p
    

    
    #--------- IMAGE FILTERS ---------#
    # Customized image filters used to sharpen or smooth, as well as threshold images.
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

    
    #---------- SEGMENTATION AND TRACKING ALGORITHMS -----------#
    # These functions are used to segment and track the particles

    def particle_segmentation(self, back_sub_image, log_adaptive_parameters, min_spot_size, max_spot_size, min_spot_aspect_ratio, post_processing_threshold):
        """
        This function is used to segment particle in the fluorescence image.
        It first applies an log_adaptive filer (see line 883) to filter and segment the image.
        If the detected spots are larger than the max_spot_size parameter, or have a very low aspect ratio (minor axis / major axis) then the algorithm zooms into
        the spot and applies new segmentation parameters. This happens since very large or elongated spots are probably attributed to many particles recognized as one.
        By applying stricter parameters it is possible to separate the clustered particles. 
        
        Parameters
        ----------
        back_sub_image: 2D numpy ndarray: the particle fluorescence image to be segmented (background corrected)
        log_adaptive_parameters: list - the parameters for segmentation (log_adaptive filter parameters)
        min_spot_size: integer - the expected min size of the particle in pixels
        max_spot_size: integer - the expected max size of the particle in pixels
        min_spot_aspect_ratio: float - the expected minimum aspect ration (minor/major axis) of thge particle
        post_processing_threshold: float or integer - The hard threshold for segmented the clustered particles as a percentage of brightest pixels
    
        Returns
        -------
        A pandas dataframe:
            columns:
            'centroid' microtubule_centroid -> the centroid coordinates of the particle (x,y tuple)
            'minor_axis' microtubule_minor_axis -> the minor axis of the segmemted particle spot in px
            'major_axis' microtubule_major_axis -> the major axis of the segmented particle spot in px
            'aspect_ratio' particle_aspect_ratio_list -> the minor/major axis ratio
            'area' microtubule_area -> the area of the segmented particle spot in px
            'experiment' self.experiment -> the experiment (inherited from the class)
        """
        def post_processing(back_sub_image, particle_centroid_yx, particle_major_axis, particle_minor_axis, post_processing_threshold, particle_lists):
            """
            This function is used to reprocess the particle ROI in otder to separate clustered particles.
            It zooms into the segmented spot and applies a new LoG filter to separate multiple particles that
            were segmented as one (very proximal to each other).
            This function is applied if the segmented spot is larger than the maximum spot area selected by the user
            or if the aspect ratio (minor/major axis) of the spot is lower than expected.
            
            If the pos-processed regoins of interest result in smaller or rounder segmented spots within the specified ranges of 
            particle areas and aspect ratios, then the particle coordinates and dimenions lists are updated including these
            post-processed particles.
            
            Parameters
            ----------
            back_sub_image: a 2D numpy ndarray of the backgrund subtracted image
            particle_centroid_yx: tuple - the (y,x) coordinates of the segmented spot
            particle_major_axis: float - the major axis of the segmented spot
            particle_minor_axis: float - the minor axis of the segmented spot
            post_processing_threshold: float or integer - The hard threshold for segmented the clustered particles as a percentage of brightest pixels
            particle_lists: a list containing the particle_centroid, particle_minor_axis, particle_major_axis, particle_aspect_ratio and particle_area lists.
                [particle_centroid_list, particle_minor_axis_list, particle_major_axis_list, particle_aspect_ratio_list, particle_area_list]
            Returns
            -------
            The updated particle_centroid [0], particle_minor_axis [1], particle_major_axis [2], particle_aspect_ratio [3] and particle_area lists [4].
            """
            particle_centroid_list, particle_minor_axis_list, particle_major_axis_list, particle_aspect_ratio_list, particle_area_list = particle_lists
            
            if type(post_processing_threshold) == int or type(post_processing_threshold) == float:
                if int(particle_centroid_yx[0]-particle_major_axis/2) >= 0 and int(particle_centroid_yx[0]+particle_major_axis/2) < self.sensor[1] and int(particle_centroid_yx[1]-particle_major_axis/2) >= 0 and int(particle_centroid_yx[1]+particle_major_axis/2) < self.sensor[0]:
                    print(particle_centroid_yx[::-1], particle_area, particle_major_axis, particle_minor_axis,'post processing...')
                    # crop the image around the spot position - zooming into the spot of interest and apply a stricter log filter
                    cropped_image = back_sub_image[int(particle_centroid_yx[0]-particle_major_axis/2):int(particle_centroid_yx[0]+particle_major_axis/2), int(particle_centroid_yx[1]-particle_major_axis/2):int(particle_centroid_yx[1]+particle_major_axis/2)]
                    fluorescence_threshold = np.sort(cropped_image.ravel())[-int(len(cropped_image.ravel())*((100-post_processing_threshold)/100))]
                    
                    cropped_filtered_image = cropped_image > fluorescence_threshold
                    # emumerate the segmented spots in the cropped image
                    sub_particle_labels = label(cropped_filtered_image)
                    # get the features of the segmented spots
                    for sub_particle_region in regionprops(sub_particle_labels):
                        sub_particle_centroid_yx = sub_particle_region.centroid
                        sub_particle_centroid_yx_corrected = (sub_particle_centroid_yx[0] + int(particle_centroid_yx[0]-particle_major_axis/2), sub_particle_centroid_yx[1] + int(particle_centroid_yx[1]-particle_major_axis/2))
                        sub_particle_minor_axis = sub_particle_region.minor_axis_length
                        sub_particle_major_axis = sub_particle_region.major_axis_length
                        sub_particle_area = sub_particle_region.area
                        
                        if sub_particle_area >= min_spot_size  and sub_particle_minor_axis > 0 and sub_particle_major_axis > 0:
                            sub_particle_aspect_ratio = sub_particle_minor_axis / sub_particle_major_axis

                            particle_centroid_list.append(sub_particle_centroid_yx_corrected[::-1])
                            particle_minor_axis_list.append(sub_particle_minor_axis)
                            particle_major_axis_list.append(sub_particle_major_axis)
                            particle_aspect_ratio_list.append(sub_particle_aspect_ratio)
                            particle_area_list.append(sub_particle_area)
            # Only if the post-processing threshold is a % float or integer of brightest pixels it is then used to split clustered spots within the region
            elif post_processing_threshold == None:
                print(particle_centroid_yx[::-1], particle_area, particle_major_axis, particle_minor_axis, 'the post-processing parameters were not defined. post-processing aborted...')
            
            return particle_centroid_list, particle_minor_axis_list, particle_major_axis_list, particle_aspect_ratio_list, particle_area_list

        print('Filtering image...')
        filtered_image = self.log_adaptive_filter(back_sub_image, log_adaptive_parameters)[2]
        # ennumerate the separated masks in the image using the labels functioon from the skimage.measure library
        particle_labels = label(filtered_image)
        # These lists will store the particle data if the particles satisfy all the conditions
        particle_centroid_list = []
        particle_minor_axis_list = []
        particle_major_axis_list = []
        particle_aspect_ratio_list = []
        particle_area_list = []
    
        for particle_label in regionprops(particle_labels):
            particle_centroid_yx = particle_label.centroid
            particle_minor_axis = particle_label.minor_axis_length
            particle_major_axis = particle_label.major_axis_length
            particle_area = particle_label.area
            # This condition is important to avoid numerical erros (very small particles of 1 or 2 pixels in size) and avoids division by zero
            if (particle_area > 0 and particle_minor_axis > 0 and particle_major_axis > 0):
                particle_aspect_ratio = particle_minor_axis/particle_major_axis
                if particle_area >= min_spot_size and particle_area <= max_spot_size and particle_aspect_ratio >= min_spot_aspect_ratio:
                    # reversing the centroid ro correspond to the xy and not the yx coordinates
                    particle_centroid_list.append(particle_centroid_yx[::-1])
                    particle_minor_axis_list.append(particle_minor_axis)
                    particle_major_axis_list.append(particle_major_axis)
                    particle_aspect_ratio_list.append(particle_aspect_ratio)
                    particle_area_list.append(particle_area)
                elif particle_area >= min_spot_size and particle_area <= max_spot_size and particle_aspect_ratio < min_spot_aspect_ratio:
                    print(particle_area, particle_minor_axis, particle_major_axis)
                    particle_lists = [particle_centroid_list, particle_minor_axis_list, particle_major_axis_list, particle_aspect_ratio_list, particle_area_list]
                    particle_centroid_list, particle_minor_axis_list, particle_major_axis_list, particle_aspect_ratio_list, particle_area_list = post_processing(back_sub_image, particle_centroid_yx, particle_major_axis, particle_minor_axis, post_processing_threshold, particle_lists)
                elif (particle_area > max_spot_size):
                    print(particle_area, particle_minor_axis, particle_major_axis)
                    particle_lists = [particle_centroid_list, particle_minor_axis_list, particle_major_axis_list, particle_aspect_ratio_list, particle_area_list]
                    particle_centroid_list, particle_minor_axis_list, particle_major_axis_list, particle_aspect_ratio_list, particle_area_list = post_processing(back_sub_image, particle_centroid_yx, particle_major_axis, particle_minor_axis, post_processing_threshold, particle_lists)
            elif (particle_area <= 0 or particle_minor_axis <= 0 or particle_major_axis <= 0):
                print(particle_centroid_yx[::-1], particle_area, particle_major_axis, particle_minor_axis, 'Error: division by zero, spot aborted')
        # create a pandas dataframe with the segmented particle coordinates and statistics per frame
        particle_segmentation_df = pd.DataFrame()
        particle_segmentation_df['centroid'] = particle_centroid_list
        particle_segmentation_df['minor_axis'] = particle_minor_axis_list
        particle_segmentation_df['major_axis'] = particle_major_axis_list 
        particle_segmentation_df['area'] = particle_area_list
        particle_segmentation_df['aspect_ratio'] = particle_aspect_ratio_list
        particle_segmentation_df['experiment'] = self.experiment
        
        return particle_segmentation_df
    
    
    def check_particle_segmentation_parameters(self, frame, post_processing=False):
        """
        This code is used to check the segmentation parametetrs.
        The image uneven background is subtracted. 
        Then the user is asked to input the particle segmentation parameters. Default values are recommended.
        All the segmentation parameters are used in the self.particle_segmentation function (line 933)
        
        The segmentation is shown in the particle fluorescence channel.
        
        Parameters
        ----------
        frame: non-negative integer - the frame in the fast time-lapse to be checked
        post_processing: binary - if True the user is asked to provide post-processing parameters

        Returns
        -------
        [0] log_adaptive_parameters
        [1] the maximum particle area
        [2] the minimum particle aspect ratio
        [3] log_adaptive post_processing parameters (if post_processing is False an empty list is returned.
        Providing an empty list in the self.particle_segmentation function will abort post-processing)
        """
        # good default parameters
#        ([4.0, 1000, 98.5, 2.0, 7, -2.0, 0], 3, 150, 0.3, [])
        
        print('The particles image in frame:', frame)
        back_sub_image = self.particle_images[frame]

        i = 0
        while i == 0:
            try:
                gaussian_param = float(input('Choose the gaussian smoothing factor (recommended: 4):'))
            except ValueError:
                print('please choose a number')
                gaussian_param = float(input('Choose the gaussian smoothing factor (recommended: 4):'))
            
            try:
                laplace_param = int(input('Choose the laplace filter factor (recommended: 1000):'))
            except ValueError:
                print('please choose a number')
                laplace_param = int(input('Choose the laplace filter factor (recommended: 1000):'))
            
            try:
                log_threshold = float(input('Choose the hard threshold of the LoG filter (recommended: 97):'))
            except ValueError:
                print('please choose a number')
                log_threshold = float(input('Choose the gaussian smoothing factor (recommended: 97):'))
            
            try:
                adaptive_smoothing = float(input('Choose the gaussian smoothing factor before the adaptive thresholding (recommended: 2):'))
            except ValueError:
                print('please choose a number')
                adaptive_smoothing = float(input('Choose the gaussian smoothing factor before the adaptive thresholding (recommended: 2):'))
            
            try:
                block_size_in = int(input('Choose the block size for the adaptive smoothing (recommended: 9, odd number):'))
            except ValueError:
                print('please choose an odd integer positive number')
                block_size_in = int(input('Choose the block size for the adaptive smoothing (recommended: 9, odd number):'))
            
            try:
                offset_in = float(input('Choose the offset for the adaptive smoothing (recommended: -2):'))
            except ValueError:
                print('please choose a number')
                offset_in = float(input('Choose the offset for the adaptive smoothing (recommended: -2):'))
            
            try:
                erosion_in = int(input('Choose the number of erosion rounds for the adaptively thresholded mask (reccomdned: 0):'))
            except ValueError:
                print('please choose a non-negative integer')
                erosion_in = int(input('Choose the number of erosion rounds for the adaptively thresholded mask (reccomdned: 0):'))
            
            try:
                min_particle_size = int(input('Choose the expected minimum size of the particle (recommended: 3):'))
            except ValueError:
                print('please choose a number')
                min_particle_size = int(input('Choose the expected minimum size of the particle (recommended: 3):'))
            
            try:
                max_particle_size = int(input('Choose the expected maximum size of the particle (recommended: 150):'))
            except ValueError:
                print('please choose a number')
                max_particle_size = int(input('Choose the expected maximum size of the particle (recommended: 150):'))
            
            try:
                min_particle_aspect_ratio = float(input('Choose the expected minimum aspect ratio of the aprticle (recommended: 0.3):'))
            except ValueError:
                print('please choose a number')
                min_particle_aspect_ratio = float(input('Choose the expected minimum aspect ratio of the aprticle (recommended: 0.3):'))
            
            if post_processing == True:
                try:
                    post_processing_threshold = float(input('Choose the post processing threshold (recommended: 90):'))
                except ValueError:
                    print('please choose a number')
                    post_processing_threshold = float(input('Choose the post processing threshold (recommended: 90):'))
            elif post_processing ==False:
                post_processing_threshold = None
                

            log_adaptive_parameters = [gaussian_param, laplace_param, log_threshold, adaptive_smoothing, block_size_in, offset_in, erosion_in]
            particle_segmentation_df = self.particle_segmentation(back_sub_image, log_adaptive_parameters, min_particle_size, max_particle_size, min_particle_aspect_ratio, post_processing_threshold)
            filtered_image = self.log_adaptive_filter(back_sub_image, log_adaptive_parameters)
            
            print('The LoG image')
            plt.figure(figsize=(20,20))
            plt.imshow(filtered_image[0])
            plt.show()
            
            print('Press enter to continue')
            input()
            
            print('The adaptively thresholded image')
            plt.figure(figsize=(20,20))
            plt.imshow(filtered_image[1])
            plt.show()
            
            print('Press enter to continue')
            input()
            
            print('Combining the LoG hard thresholded and the adaptively thresholded images')
            plt.figure(figsize=(20,20))
            plt.imshow(filtered_image[2])
            plt.show()
            
            print('Press enter to continue')
            input()
            
            # Print a figure with all the segmented particles and get the particle centers
            print('showing segmentation...')
            fig, ax = plt.subplots(figsize=(20, 20))
            # The vmin and vmax values can be adjusted to change the LUTs, as well as the colormap
            ax.imshow(back_sub_image, cmap='viridis', vmin=0, vmax=300)
            
            for index, row in particle_segmentation_df.iterrows():
                rect = matplotlib.patches.Rectangle((row['centroid'][0]-row['major_axis']/2,row['centroid'][1]-row['major_axis']/2), row['major_axis'], row['major_axis'], fill=False, edgecolor='coral', linewidth=1)
                ax.add_patch(rect)
                ax.set_axis_off()
                plt.tight_layout()
            
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
        
        return log_adaptive_parameters, min_particle_size, max_particle_size, min_particle_aspect_ratio, post_processing_threshold

    
    def show_gaussian_particle_fit(self, cropped_particle_image, fitted_subpx, fitted_px,  param):
        """
        This function plots the particle fluorescence 2D gaussian fit and the particle center estimation.
        
        Parameters
        ----------
        cropped_particle_image: 2D numpy array - the cropped particle image (within the bounding box)
        fitted_subpx: 2D numpy array - the blurred particle image
        fitted_px: 2D numpy array - the raw particle image
        param: height, y, x, width_y, width_x, rotation - parameters of the 2D Gaussian fit
        """
        (height, y, x, width_y, width_x, rotation) = param
        # PLOT 1
        fig, ax = plt.subplots(figsize=(10, 5))
        plt.imshow(fitted_subpx)
        ax = plt.gca()
        plt.text(0.95, 0.05, """
                 amp: %.1f
                 x : %.1f
                 y : %.1f
                 width_x : %.1f
                 width_y : %.1f""" %(height, x, y, width_x, width_y),
                 fontsize=16, horizontalalignment='right',
                 verticalalignment='bottom', transform=ax.transAxes)
        plt.show()
        # PLOT 2
        # Fit a citcle at the center of each particle, which corresponds to the peak of the Gaussian
        circle = matplotlib.patches.Circle((x, y), radius=0.1, color='red', fill=False, linewidth=3)
        fig, ax = plt.subplots(figsize=(10, 5))
        plt.imshow(cropped_particle_image)
        plt.contour(fitted_px, cmap=plt.cm.Greys)
        ax.add_patch(circle)
        print('2D gaussian params:',param)
        plt.show()  
    

    def estimate_particle_center(self, bkg_cor_image, particle_center, box_size, gaussian_fit_show=False):
        """
        This function fits the 2D Gaussian to estimate the particle center.
        
        Parameters
        ----------
        bkg_cor_image: 2D numpy array - the background corrected particle fluorescence image
        particle_center: tuple of floats (x,y) corresponding to the center of mass of the particle label
        box_size: odd integer - the size of the box that is used to fit the 2D gaussian for the particle center estimation (5 or 11 suggested)
        gaussian_fit_show: bool - True to show the 2D Gaussian fit
                                  False otherwise
        Returns
        -------
        gaussian_particle_center: tupel of (x,y) floats: the estimate particle center (subpixel)
        brightest_raw_pixels: numpy array - the values of the 60% brightest raw pixels
        brightest_fitted_pixels: numpy array - the values of the 60% brightest smoothed pixels
            returns 'none' if one of the following exceptions is met.
        
        Exception
        ---------
        If the approximated particle center falls outside the bounding box, 
                the spot is not diffratcion limited and the particle is aborted.
                        'none' is returned
        If the particle bounding box falls outside the sensor dimensions, 
            'none' is returned
        If the particle eccentricity is below 0.1 or above 19,
            This spot probably corresponds to multiple clustered particles.
                'none' is returned
        """
        half_side = int((box_size-1)/2)
        cropped_particle_image = bkg_cor_image[(int(particle_center[1])-half_side):(int(particle_center[1])+(half_side+1)), (int(particle_center[0])-half_side):(int(particle_center[0])+(half_side+1))]
        if (cropped_particle_image.shape[0] > 0 and cropped_particle_image.shape[1] > 0):
            # create a meshgrid with 0.1 pixel resolution
            Xin, Yin = np.mgrid[0:cropped_particle_image.shape[1]:0.1,0:cropped_particle_image.shape[0]:0.1]
            # create a meshgrid with single pixel resolution
            Xin2, Yin2 = np.mgrid[0:cropped_particle_image.shape[1]:1,0:cropped_particle_image.shape[0]:1]
            # fit a 2D gaussian with rotation to the cropped image
            try:
                param = self.fitgaussian(cropped_particle_image)
                (height, y, x, width_y, width_x, rotation) = param
                # set an eccentricity threshold and the center coordinates fall within the box
                if (width_y/width_x > 0.1 and width_y/width_x < 10 and x<box_size and y<box_size and x>0 and y>0):
                    gaussian_fit = self.gaussian(*param)
                    # fit a guassian to a lattice of Xin by Yin pixels
                    fitted_subpx = gaussian_fit(Xin, Yin)
                    fitted_px = gaussian_fit(Xin2, Yin2)
                    # Get the 60% brightest smoothed pixels
                    brightest_fitted_pixels = fitted_px[fitted_px>np.percentile(fitted_px,40)]
                    # Get the 60% brightest raw pixels
                    brightest_raw_pixels = cropped_particle_image[cropped_particle_image>np.percentile(cropped_particle_image, 40)]
                    # Correct the gaussian coordinates to the original image dimensions (from the cropped dimensions)
                    gaussian_x_coord = x + int(particle_center[0])-half_side
                    gaussian_y_coord = y + int(particle_center[1])-half_side
                    gaussian_particle_center = (gaussian_x_coord, gaussian_y_coord)
                    if gaussian_fit_show == True:
                        self.show_gaussian_particle_fit(cropped_particle_image, fitted_subpx, fitted_px,  param)
            
                    return gaussian_particle_center, brightest_raw_pixels, brightest_fitted_pixels, param
                else:
                    print('This particle did not pass the eccentricity and position conditions. Particle position aborted:', particle_center)
                    return 'none'
            except IndexError:
                print('The approximate particle center is outside the square bounds. Particle position aborted:', particle_center)
                return 'none'
        else:
            print('This particle position spec ranges out of bounds. Particle position aborted:', particle_center)
            return 'none'

    
    def get_particle_fluorescence(self, metric, operation, brightest_raw_pixels, brightest_fitted_pixels, gaussian_vol):
        """
        This function is used to estimate the particle fluorescence at each position
        
        Paramters
        ---------
        metric: string, the column of the pandas dataframe which is used as a proxy for particle size. 
                Choose 'gaussian volume', 'raw pixels' or 'smoothed pixels'
        operation: function, the mathematical operation applied to the bin_column as a particle size proxy. 
                    Choose 'mean', 'median', 'sum' or 'max'
        brightest_raw_pixels: numpy array - the 60% brightest raw pixels 
        brightest_fitted_pixels: numpy array - the 60% brightest smoothed pixels 
        gaussian_vol: the volume of the fitted 2D Gaussian to the particle pixels
        
        Returns
        -------
        particle_vol: float - the particle fluorescence value corresponding to the chosen metric and operation
        
        Raises
        ------
        ValueError if an invalid metric or operation are included as inputs
        """
        # estimate the particle volume statistic based on a given metric and operation
        if metric == 'gaussian volume':
            particle_vol = gaussian_vol.copy()
            operation = 'none'
        elif metric == 'raw pixels':
            if operation == 'mean':
                particle_vol = np.mean(brightest_raw_pixels)
            elif operation == 'median':
                particle_vol = np.median(brightest_raw_pixels)
            elif operation == 'sum':
                particle_vol = np.sum(brightest_raw_pixels)
            elif operation == 'max':
                particle_vol = np.max(brightest_raw_pixels)
            else:
                raise ValueError("Choose 'mean', 'median', 'sum' or 'max' as the operation.")
        elif metric == 'smoothed pixels':
            if operation == 'mean':
                particle_vol = np.mean(brightest_fitted_pixels)
            elif operation == 'median':
                particle_vol = np.median(brightest_fitted_pixels)
            elif operation == 'sum':
                particle_vol = np.sum(brightest_fitted_pixels)
            elif operation == 'max':
                particle_vol = np.max(brightest_fitted_pixels)
            else:
                print("wrong operation for particle fluorescence estimation")
                raise ValueError("Choose 'mean', 'median', 'sum' or 'max' as the operation.")
        else:
            print("wrong metric for particle fluorescence estimation")
            raise ValueError("Choose 'gaussian volume', 'raw pixels', or 'smoothed pixels' as a metric")
        
        return particle_vol
        

    
    def getting_the_particles(self, log_adaptive_parameters, min_particle_size, max_particle_size, min_particle_aspect_ratio, 
                              post_processing_threshold, box_size, 
                              analysis_range, metric, operation,  
                              gaussian_fit_show = False):
        """
        This function is used to run the segmentation of the particles and track them.
        Gaussian distributions are fitted to the difraction limited particles to estimate the volume and the center of each particle. 
        
        Parameters
        ----------
        analysis_range: tuple of non-negaive integers: the first and the last frame to be analyzed
        log_adaptive_parameters: list of the particle segmentation parameters. See particle_segmentation function.
        min_particle_size: Non-negative integer: the minimum expected particle size in pixels. 
        max_particle_size: Positive integer: the maximum expected particle size in pixels.
        min_particle_aspect_ratio: float: the minimum expected particle aspect ration (minor/major axis length)
        post_processing_threshold: float or integer: the % threshold of brightest pixels to separate clustered spots.
        box_size: odd integer - the size of the box that is used to fit the 2D gaussian for the particle center estimation (5 or 11 suggested)
        analysis_range: tuple - the frame range of the analysis 
        metric: string, the column of the pandas dataframe which is used as a proxy for particle size. Choose 'gaussian volume', 'raw pixels' or 'smoothed pixels'
        operation: function, the mathematical operation applied to the bin_column as a particle size proxy. Choose 'mean', 'median', 'sum' or 'max'
        
        Returns
        -------
        A dataframe with the following fields:
            'experiment' - A string with the experiment ID
            'position' - Non-negative integer showing the xy position of the experiment
            'label' - the cell label
            'frame' - the time frame in the analysis (Non-negative integer)
            'particle_center' - the particle center approximated during particle segmentation (particle mask centroid)
            'max_fluorescence' - the maxmimum fluorescence of the particle estimated from the particle mask
            'particle_brightest_pixels' - the 60% of the particle pixels around the particle centroid
            'smoothed_brightest_pixels' - the 60% brightest smoothed particle pixels (after fitting a 2D Gaussian) around the particle centroid
            'gaussian_center' - the center of the particle estimated from the peak of the 2D gaussian fitted
            'gaussian_amplitude' - the amplitude of the fitted 2D gaussian
            'gaussian_std' - the std of the fitted 2D gaussian in the x and y dimension (tuple)
            'gaussian_volume' - the volume of the fitted 2D gaussian given by this function (2 π A σx σy) - https://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function
            'gaussian_rotation' - the rotation of the fitted 2D gaussian
            'particle_fluorescence' - the proxy of particle size estimated using a mathematical operation (parameter: operation) on a given particle metric (parameter: metric)
        
        This dataframe is also saved in the designated folder:
            self.save_path+'/'+self.experiment+'_'+self.position_string+'_particle_df'
        
        Raises
        -----
        ValueError if the box_size is not an odd integer (for fitting the 2D Gaussian)
                   if the metric input is not valid (for the estimation of the particle fluorescence)
                   if the operation input is not valid (for the estimation of the particle fluorescence)
        """
        if box_size%2 == 0:
            raise ValueError('The box_size parameter is even. Choose an odd integer.')
        if metric not in ['gaussian volume', 'raw pixels', 'smoothed pixels']:
            raise ValueError("Choose 'gaussian volume', 'raw pixels', or 'smoothed pixels' as a metric")
        if operation not in ['mean', 'median', 'sum', 'max']:
             raise ValueError("Choose 'mean', 'median', 'sum' or 'max' as the operation.")
        
        if analysis_range[0] == 0:
            # Initiate a pandas dataframe were all the results will stored
            particle_df = pd.DataFrame(columns=['experiment', 'position', 'label', 'frame', 'particle_center', 'max_fluorescence', 'particle_brightest_pixels', 'smoothed_brightest_pixels', 'gaussian_center', 'gaussian_amplitude', 'gaussian_std', 'gaussian_volume', 'gaussian_rotation', 'particle_fluorescence'])
            # resume the analysis from a previous frame if it exists, otherwise start a new analysis from a non-zero frame
        elif analysis_range[0] > 0:
            if os.path.exists(self.save_path+'/'+self.experiment+'_'+self.position_string+'_particles_df') == True:
                particle_df = pd.read_pickle(self.save_path+'/'+self.experiment+'_'+self.position_string+'_particles_df', compression='zip')
            else:
                particle_df = pd.DataFrame(columns=['experiment', 'position','label', 'frame', 'particle_center', 'max_fluorescence', 'particle_brightest_pixels', 'smoothed_brightest_pixels', 'gaussian_center', 'gaussian_amplitude', 'gaussian_std', 'gaussian_volume', 'gaussian_rotation', 'particle_fluorescence'])
        # Iterate through the specified frames in the stream acquisition
        for fr in range(analysis_range[0], analysis_range[1]):
            print('frame '+str(fr+1)+' out of '+str(len(self.particle_images))+', position: '+self.position_string)
            # get the signal image
            bkg_cor_image = self.particle_images[fr]

            #--------- PARTICLE SEGMENTATION ----------#
            particle_segmentation_df = self.particle_segmentation(bkg_cor_image, log_adaptive_parameters, min_particle_size, max_particle_size, min_particle_aspect_ratio, post_processing_threshold)
            for index, row in particle_segmentation_df.iterrows():
                particle_center = row['centroid']
                # get the particle center and the 2D Gaussian parameters
                particle_center_stats = self.estimate_particle_center(bkg_cor_image, particle_center, box_size, gaussian_fit_show)
                if particle_center_stats != 'none':
                    gaussian_particle_center, brightest_raw_pixels, brightest_fitted_pixels, param = particle_center_stats
                    cell_label = self.cell_labels[fr][int(round(gaussian_particle_center[1],0))][int(round(gaussian_particle_center[0],0))]
                    (height, y, x, width_y, width_x, rotation) = param
                    gaussian_vol = 2*np.pi*height*width_x*width_y
                    # get the particle fluorescence metric
                    particle_vol = self.get_particle_fluorescence(metric, operation, brightest_raw_pixels, brightest_fitted_pixels, gaussian_vol)
                    # organize the data into a new pandas row
                    pre_df = [self.experiment, self.position, cell_label, fr, particle_center, brightest_raw_pixels.max(), brightest_raw_pixels, brightest_fitted_pixels, gaussian_particle_center, height, (width_x, width_y), gaussian_vol, rotation, particle_vol]
                    # append the row to the pandas dataframe
                    particle_df = particle_df.append(pd.DataFrame([pre_df], columns=list(particle_df.columns)),ignore_index=True)

        with open(self.save_path+'/'+self.experiment+'_'+self.position_string+'_particles_df', 'wb') as handle:
            particle_df.to_pickle(path = handle, compression='zip', protocol = pickle.HIGHEST_PROTOCOL)
        
        return particle_df
    
    

    
    
    