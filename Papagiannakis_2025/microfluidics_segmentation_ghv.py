# -*- coding: utf-8 -*-
"""
Created on Wed May 26 14:23:37 2021

@author: Alexandros Papagiannakis, Christine Jacobs-Wagner lab, Sarafan ChEM-H, Stanford University 2021
"""

import scipy.ndimage as ndi
import scipy.signal as scisig
import numpy as np
import matplotlib.pyplot as plt
from skimage import filters, io, segmentation, measure
import os
import pickle
import pandas as pd
import warnings
warnings.filterwarnings('ignore')


class mother_machine_segmentation(object):
    """
    Developer: Alexandros Papagiannakis, Stanford University, May 2021
    
    This class contains all the functions used for cell segmentration and tracking in the mother machine
    """
    def __init__(self, deback_phase_path, deback_fluor_path, phase_interval, experiment, position, save_path):
        """
        The class is iniialized.
        
        Input:
            ___File paths___
            deback_phase_path: string - the path to the folder with the background corrected phase images
            deback_fluor_path: string - the path to the folder with the fluorescent images
            
            ___General parameters___
            experiment: a string describing the experiment (e.g. '07202919_seqA')
            position: a non-negatice integer corresponding to the channel position in the experimenbt (e.g. 1 for the first position)
            phase_interval: the time interval of the time lapse for the phase channel in min
            save_path: the directory where the segmentation results will be saved
            
                
            The _init_ function loads all the phase and fluorescence images and stores them in the class:
                self.phase_images
                self.phase_interval
                self.experiment
                self.position
                self.position_string
                self.save_path
                self.n_frames
        """
        
        def get_images_dict(image_path):
            images_names = os.listdir(image_path)
            images_dict = {}
            i = 0
            for img in images_names:
                images_dict[i] = io.imread(image_path + '/' + img)
                i += 1
            return images_dict
            
#        images_names = os.listdir(deback_phase_path)
#        phase_dict = {}
#        i = 0
#        for img in images_names:
#            phase_dict[i] = io.imread(deback_phase_path + '/' + img)
#            i += 1

        self.phase_images = get_images_dict(deback_phase_path)
        if os.path.isdir(deback_fluor_path):
            self.fluor_images = get_images_dict(deback_fluor_path)
        else:
            self.fluor_images = 'none'
            
        
#        self.phase_images = phase_dict
        self.phase_interval = phase_interval
#        self.fluor_interval = fluor_interval
        self.experiment = experiment
        self.position = position
#        self.n_frames = len(phase_dict)
        self.n_frames = len(self.phase_images)
        self.save_path = save_path
        
                
        if position <= 9:
            position_string = 'pos0'+str(position)
        elif position > 9:
            position_string = 'pos'+str(position)
            
        self.position_string = position_string

    
    def get_median_filtered(self, projection, difference_threshold=7):
        """
        This function is used to locate the positions where the sausage is cut.
        """
        signal = projection.copy()
        difference = np.abs(signal - np.median(signal))
        median_difference = np.median(difference)
        if median_difference == 0:
            s = 0
        else:
            s = difference / float(median_difference)
        mask = s > difference_threshold
        
        return mask
    
    
    def adaptive_filter(self, image, gaussian_sigma, block_size, offset):
        """
        Implementation of the scipy.ndimage adaptive filter
        """
        image_smooth = ndi.gaussian_filter(image, sigma = gaussian_sigma)
        adaptive_threshold = filters.threshold_local(image_smooth, block_size=block_size, offset=offset)
        return image_smooth>adaptive_threshold
    

    def log_filter(self, image, gaussian_sigma, log_kernel, percentile):
        """
        Implementation of the scipy.ndimage Laplace of Gaussian filter
        """
        image_smooth = ndi.gaussian_filter(image, sigma = gaussian_sigma)
        image_laplace = filters.laplace(image_smooth, log_kernel)
        log_threshold = np.percentile(image_laplace.ravel(), percentile)
        return image_smooth>log_threshold
    
    
    def watershed_function(self, image_mask, markers='none', min_distance=4):
        """
        Custom implementation of the scipy.ndimage watershedding
        """
        distance = ndi.distance_transform_edt(image_mask)
        if markers=='none':
            markers = distance * (distance>=min_distance)
            markers, _ = ndi.label(markers)
        labels = segmentation.watershed(-distance, markers, mask=image_mask)
        return labels
    
    
    def get_max_horizontal_projection_range(self, image, window):
        """
        Get the center of the channel where the channel background corrected phase signal is the brightest
        """
        x_max = image.shape[1]
        vertical_mean_projection = []
        for x_win in range(0,x_max-1-window):
            croped_image = image[:, x_win:(x_win+window)]
            vertical_mean_projection.append(croped_image.mean())
        return np.argmax(vertical_mean_projection) + int(window/2)


    def get_vertical_projection(self, image, line_center, half_window, function):
        """
        Get the vertical projection of the channel.
        """
        croped_image = image[:,(line_center-half_window):(line_center+half_window+1)]
        if function == 'mean':
            product = np.mean(croped_image, axis=1)
        elif function == 'median':
            product = np.median(croped_image, axis=1)
        elif function == 'min':
            product = np.min(croped_image, axis=1)
        elif function == 'max':
            product = np.max(croped_image, axis=1)
        return product


    def background_correction(self, image, window=10, orientation='down'):
        """
        Correct the phase background using the top-left for the 'down' orientation,
        or the bottom-right corner of 10x10 pixels in the 'up' orientation
        """
        if orientation == 'down':
            bkr = image[0:window, 0:window].mean()
        elif orientation == 'up':
            bkr = image[(image.shape[0]-window):, (image.shape[1]-window):].mean()
        image_cor = image-bkr
        image_cor[image_cor<0]=0
        return image_cor


    def get_upper_lower_bounds(self, image):
        """
        Get the bounds of the Otsu mask used to crop the cell mask.
        """
        otsu_mask = image>filters.threshold_otsu(image)
        y_coords, x_coords = np.nonzero(otsu_mask)
        return np.min(y_coords), np.max(y_coords)


    def get_sausage_troughs(self, frame, orientation, channel_end, channel_start, 
                            central_line_window, smoothing_factor, projection_width,
                            std_threshold, show_threshold, save_threshold_path='none'):
        """
        Get the troughs of the background corrected phase intensity which is used to split the dividing cells
        """
        phase_image = self.phase_images[frame]
        image_cor = self.background_correction(phase_image)  # correct the image background
        
        upper_y, lower_y = self.get_upper_lower_bounds(image_cor) # get the upprt and lower bounds of the channel
        if orientation == 'down':
            lower_y = lower_y - channel_end # if the channel is deformed towards the edge, remove a number of pixels from the edge of the mask
            upper_y = upper_y + channel_start # add 10 pixels to avoid the cell edges
        elif orientation == 'up':
            upper_y = upper_y + channel_end # if the channel is deformed towards the edge, remove a number of pixels from the edge of the mask
            lower_y = lower_y - channel_start # add 10 pixels to avoid the cell edges
        
        crop_image = image_cor[upper_y:lower_y, :] # crop the image using the upper/lower bounds
        central_line_position = self.get_max_horizontal_projection_range(crop_image, window=central_line_window) # get the central line with the max hirozontal projection going through the center of the channel
        smooth_image = ndi.gaussian_filter(crop_image, sigma=smoothing_factor) # smooth the phase images
        projection = self.get_vertical_projection(smooth_image, line_center=central_line_position, half_window=projection_width, function='mean') # get the mean projection along the channel
        troughs = scisig.find_peaks(1/projection) # get the troughs of the prokection
#            thresholded_projection = projection - (np.mean(projection)-std_threshold*np.std(projection))
        thresholded_projection = projection - (np.mean(projection)-std_threshold*np.std(projection)) # set a threshold using the mean and the std of the projection
        curated_troughs = troughs[0][np.argwhere(thresholded_projection[troughs[0]]<0)] # only the troughs below the degined dynamic threshold are considered
        
        if show_threshold == True:
            self.visualize_segmentation_threshold(frame, thresholded_projection, 
                                     curated_troughs, save_threshold_path)
        
        curated_troughs = curated_troughs + upper_y # correct the trough coordinates for the cropping mask
        
        return curated_troughs


    def get_watershed_labels(self, frame, orientation, channel_end, channel_start, 
                            central_line_window, smoothing_factor, projection_width,
                            std_threshold, adaptive_smooth, adaptive_window, adaptive_offset,
                            erosions, distance_threshold, area_threshold,
                            show_labels, show_threshold, save_labels_path='none', save_threshold_path='none'):
        """
        Apply the watershed function to the masked cells, using the labels from the intensity thresholding 
        and the distance threshold.
        """
        phase_image = self.phase_images[frame]
        image_cor = self.background_correction(phase_image)  # correct the image background
        
        curated_troughs = self.get_sausage_troughs(frame, orientation, channel_end, channel_start, 
                                                   central_line_window, smoothing_factor, projection_width,
                                                   std_threshold, show_threshold, save_threshold_path)
        
        curated_troughs_2d = np.reshape(curated_troughs, (-1,1)) # transform the 1D to a 2D numpy array to index the image rows (along the Y dimension)
        otsu_mask = image_cor>filters.threshold_otsu(image_cor) # apply the Otsu threshold to get the cell mask
        adaptive_mask = self.adaptive_filter(image_cor, adaptive_smooth, adaptive_window, adaptive_offset)   
        markers = adaptive_mask.copy()
        markers[curated_troughs_2d]=False # exclude the rows at the troughs from the cell mask
        markers = ndi.morphology.binary_erosion(markers, iterations=erosions) # erode the cell masks to generate the labels for watershedding
        distance = ndi.distance_transform_edt(adaptive_mask) # get the distance transformation of the cell masks
        distance_mask = distance > distance_threshold # use a distance threshold to generate the labels for watershedding
        markers = distance_mask * markers # combine the distance with the trough separated masks
        markers, _ = ndi.label(markers)
        watershed_labels = self.watershed_function(image_mask=otsu_mask, markers=markers, min_distance='none')
        for cell_label in range(1, watershed_labels.max()+1):
            cell_area = watershed_labels[watershed_labels==cell_label].shape[0]
            if cell_area < area_threshold:
                watershed_labels[watershed_labels==cell_label]=0

        if show_labels == True:
            self.visualize_labels(background_corrected_phase=image_cor, frame=frame, 
                                  watershed_labels=watershed_labels, save_fig_path=save_labels_path)
        
        return watershed_labels


    def visualize_labels(self, background_corrected_phase, frame, watershed_labels, save_fig_path='none'):
        """
        This function can be used to visualize the watershedding labels for each specified frame.
        """
        fig = plt.figure(figsize=(3,10))
        gs = fig.add_gridspec(1, 2, hspace=0, wspace=0)
        (ax1, ax2) = gs.subplots(sharey='row')
        ax1.imshow(background_corrected_phase, cmap='gray')
        plt.text(x = -50, y =10, s = '0.066um/px', size=10, color='white', weight='bold')
        ax1.tick_params(axis='x',          # changes apply to the x-axis
                        which='both',      # both major and minor ticks are affected
                        bottom=False,      # ticks along the bottom edge are off
                        top=False,         # ticks along the top edge are off
                        labelbottom=False)
#            plt.imshow(markers)
        ax2.imshow(watershed_labels, cmap='jet')
#            plt.imshow(watershed_labels)
        plt.text(x = 0, y =10, s = 'time: '+str(frame*self.phase_interval)+'min', size=10, color='white', weight='bold')
        ax2.tick_params(axis='x',          # changes apply to the x-axis
                        which='both',      # both major and minor ticks are affected
                        bottom=False,      # ticks along the bottom edge are off
                        top=False,         # ticks along the top edge are off
                        labelbottom=False)
        if os.path.isdir(save_fig_path) == True:
            plt.savefig(save_fig_path+'/'+self.experiment+'_'+self.position_string+'_'+str(frame)+'_watershed_labels.jpeg')
        plt.show()
        

    def visualize_segmentation_threshold(self, frame, thresholded_projection, 
                                         curated_troughs, save_fig_path='none'):
        """
        This function can be used to visualize the dynamic phase intensity threshold for a given frame.
        """
        plt.plot(thresholded_projection, label='frame:'+str(frame))
        plt.axhline([0], color='black', linewidth=1)

        for trf in curated_troughs:
            plt.axvline(trf, color='red')
        plt.ylabel('inverted phase intensity (a.u.)')
        plt.xlabel('channel length (px)')
        plt.legend(loc='upper right')
        if os.path.isdir(save_fig_path) == True:
            plt.savefig(save_fig_path+'/'+self.experiment+'_'+self.position_string+'_'+str(frame)+'_dynamic_threshold.jpeg')
        plt.show()


    def fluorescence_thresholding(self, cell_labels, fluor_image, fluor_smoothing, fluor_threshold):
        """
        Masks the cell labels using a threshold on the cellular fluorescence (e.g. the ribosomal fluorescence).
        """
        fluor_mask = ndi.gaussian_filter(fluor_image, fluor_smoothing) > fluor_threshold
        
        return cell_labels * fluor_mask


    def check_threshold(self, frames):
        """
        This function can be used to check the phase intensity thresholding parameters for any given frame.
        """
        i = 0
        while i == 0:
            orientation = str(input('Choose "down" or "up" for the orientation of the channel:'))
            orientation = orientation.lower()
            if orientation == 'up' or orientation=='down':
                i += 1
            else:
                print('wrong input, try again.')
        i = 0
        while i == 0:
            try:
                channel_end = int(input('Choose the number of pixels from the end of the channel to be exluded:'))
                i+=1
            except ValueError:
                print('wrong input, choose a number.')
        i = 0
        while i == 0:
            try:
                channel_start = int(input('Choose the number of pixels from the start of the channel to be exluded:'))
                i+=1
            except ValueError:
                print('wrong input, choose a number.')
        i = 0
        while i == 0:
            try:
                central_line_window = int(input('Choose the window of pixels used to scan for the central line position:'))
                i+=1
            except ValueError:
                print('wrong input, choose a number.')
        i = 0
        while i == 0:
            try:
                smoothing_factor = int(input('Choose the smoothing kernel size in pixels. The smoothing is applied on the phase image:'))
                i+=1
            except ValueError:
                print('wrong input, choose a number.')
        i = 0
        while i == 0:
            try:
                projection_width = int(input('Choose the width of the vertical projection in pixels:'))
                i+=1
            except ValueError:
                print('wrong input, choose a number.')
        i = 0
        while i == 0:
            try:
                std_threshold = float(input('Choose the standard deviation threshold below which the troughs are considered:'))
                i+=1
            except ValueError:
                print('wrong input, choose a number.')
        
        for fr in range(frames[0], frames[1]):
            self.get_sausage_troughs(fr, orientation, channel_end, channel_start, 
                                     central_line_window, smoothing_factor, projection_width,
                                     std_threshold, show_threshold=True, save_threshold_path='none')
        
        return orientation, channel_end, channel_start, central_line_window, smoothing_factor, projection_width, std_threshold
    

    def check_watershed(self, frames, threshold_parameters):
        """
        This function can be used to check the watershedding parameters for any given frame.
        """
        orientation, channel_end, channel_start, central_line_window, smoothing_factor, projection_width, std_threshold = threshold_parameters

        i = 0
        while i == 0:
            try:
                adaptive_smooth = int(input('Choose the smoothing kernel size for the adaptive filter in pixels:'))
                i+=1
            except ValueError:
                print('wrong input, choose a number.')    
        i = 0
        while i == 0:
            try:
                adaptive_window = int(input('Choose the adaptive window side size (odd number) in pixels:'))
                if adaptive_window%2 == 0:
                    print('Choose an odd positive integer.')
                elif adaptive_window%2 != 0:
                    i+=1
            except ValueError:
                print('wrong input, choose a number.')
        i = 0
        while i == 0:
            try:
                adaptive_offset = int(input('Choose the adaptive offset (negative integer):'))
                i+=1
            except ValueError:
                print('wrong input, choose a number.')
        i = 0
        while i == 0:
            try:
                erosions = int(input('Choose the number of erosions for the thresholded masks:'))
                i+=1
            except ValueError:
                print('wrong input, choose a number.')
        i = 0
        while i == 0:
            try:
                distance_threshold = int(input('Choose the pixels threshold of the distance transformed masks to divide constricting cells:'))
                i+=1
            except ValueError:
                print('wrong input, choose a number.')
        i = 0
        while i == 0:
            try:
                area_threshold = int(input('Choose the minimum area threshold of the segmented masks in pixels:'))
                i+=1
            except ValueError:
                print('wrong input, choose a number.')
       
        for fr in range(frames[0], frames[1]):
            self.get_watershed_labels(fr, orientation, channel_end, channel_start, 
                                central_line_window, smoothing_factor, projection_width,
                                std_threshold, adaptive_smooth, adaptive_window, adaptive_offset,
                                erosions, distance_threshold, area_threshold, 
                                show_labels=True, show_threshold=False, save_labels_path='none', save_threshold_path='none')
        
        return  adaptive_smooth, adaptive_window,adaptive_offset,erosions,distance_threshold,area_threshold


    def check_fluorescence_threshold(self, frames):
        """
        This function is used to check the fluorescence threshold
        """
        if self.fluor_images != 'none':
            if len(self.fluor_images) == len(self.phase_images):
                i = 0
                while i == 0:
                    try:
                        fluor_smoothing = int(input('Choose the smoothing factor of the fluorescence image:'))
                        i+=1
                    except ValueError:
                        print('wrong input, choose a number.')
                
                i = 0
                while i == 0:
                    try:
                        fluor_threshold = int(input('Choose the fluorescence intensity threshold above which cell pixels are considered:'))
                        i+=1
                    except ValueError:
                        print('wrong input, choose a number.')
                
                fluor_list = []
                for fr in range(frames[0], frames[1]):
                    fluor_mask = ndi.gaussian_filter(self.fluor_images[fr],fluor_smoothing)>fluor_threshold
                    plt.imshow(fluor_mask)
                    plt.show()
                    fluor_list.append(self.fluor_images[fr][np.nonzero(fluor_mask)].mean())
                plt.figure(figsize=(10,5))
                plt.plot(fluor_list, 'o')
                plt.xlabel('Frame')
                plt.ylabel('Mean fluorescence per channel')
                plt.axhline(fluor_threshold)
                plt.show()
                
            elif len(self.fluor_images) < len(self.phase_images):
                print('There is not a single fluorescence image for every phase image.')
                print('The fluorescence threshold was set to zero...')
                fluor_threshold = 0
                fluor_smoothing = 0
        else:
            print('No fluorescence images were loaded. The fluorescence threhsold was set to zero...')
            fluor_threshold = 0
            fluor_smoothing = 0
        
        return fluor_threshold, fluor_smoothing
    

    def check_segmentation(self, frames):
        """
        This function can be used to check all the parameters and the cell segmentation.
        """
        if os.path.isfile(self.save_path+'/'+self.experiment+'_'+self.position_string+'_'+'parameters_dict') == True:
            with open(self.save_path+'/'+self.experiment+'_'+self.position_string+'_'+'parameters_dict', 'rb') as handle:
                parameters_dict = pickle.load(handle)
        else:
            parameters_dict = {}
        print('Parameters dictionary definition: "parameters_dict" for a range of frames:', frames)
        parameters_dict[frames] = {}
        
        i = 0
        while i == 0:  
            threshold_parameters = self.check_threshold(frames)
            d = 0
            while d == 0:
                decision = str(input('Are these parameters good? yes/no:'))
                decision = decision.lower()
                if decision == 'yes':
                    i+=1
                    d+=1
                    parameters_dict[frames]['orientation'] = threshold_parameters[0]
                    parameters_dict[frames]['channel_end'] = threshold_parameters[1]
                    parameters_dict[frames]['channel_start'] = threshold_parameters[2]
                    parameters_dict[frames]['central_line_window'] = threshold_parameters[3]
                    parameters_dict[frames]['smoothing_factor'] = threshold_parameters[4]
                    parameters_dict[frames]['projection_width'] = threshold_parameters[5]
                    parameters_dict[frames]['std_threshold'] = threshold_parameters[6]
    
                elif decision == 'no':
                    print('Repeating the parameter selection.')
                    d+=1
                else:
                    print('wronng input. Choose "yes" or "no".')
        
        b = 0
        while b == 0:
            watershed_parameters = self.check_watershed(frames, threshold_parameters)
            c = 0
            while c == 0:
                decision = str(input('Is the segmentation good? yes/no:'))
                decision = decision.lower()
                if decision == 'yes':
                    parameters_dict[frames]['adaptive_smooth'] = watershed_parameters[0]
                    parameters_dict[frames]['adaptive_window'] = watershed_parameters[1]
                    parameters_dict[frames]['adaptive_offset'] = watershed_parameters[2]
                    parameters_dict[frames]['erosions'] = watershed_parameters[3]
                    parameters_dict[frames]['distance_threshold'] = watershed_parameters[4]
                    parameters_dict[frames]['area_threshold'] = watershed_parameters[5]
                    c += 1
                    b += 1
                elif decision == 'no':
                    e = 0
                    while e == 0:
                        decision_2 = str(input('Check all the parameters from the beginning? yes/no:'))
                        decision_2 = decision_2.lower()
                        if decision_2 == 'yes':
                            print('Re-checking thresholding.')
                            self.check_segmentation(frames)
                        elif decision_2 == 'no':
                            print('Re-checking watershedding.')
                            e += 1
                            c += 1
                        else:
                            print('wronng input. Choose "yes" or "no".')
                else:
                    print('wronng input. Choose "yes" or "no".')    
            
        h = 0
        while h == 0:
            fluorescence_threshold, fluorescence_smoothing = self.check_fluorescence_threshold(frames)
            k = 0
            while k == 0:
                decision = str(input('Is the fluorescence mask good? yes/no:'))
                decision = decision.lower()
                if decision == 'yes':
                    parameters_dict[frames]['fluor_threshold'] = fluorescence_threshold
                    parameters_dict[frames]['fluor_smoothing'] = fluorescence_smoothing
                    k += 1
                    h += 1
                elif decision == 'no':
                    l = 0
                    while l == 0:
                        decision_3 = str(input('Check all the parameters from the beginning? yes/no:'))
                        decision_3 = decision_2.lower()
                        if decision_3 == 'yes':
                            print('Re-checking fluorescence thresholding.')
                            self.check_fluorescence_threshold(frames)
                        elif decision_3 == 'no':
                            print('Re-checking the fluorescence thresholding.')
                            l += 1
                            k += 1
                        else:
                            print('wronng input. Choose "yes" or "no".')
                else:
                    print('wronng input. Choose "yes" or "no".')      
        
                            
        with open(self.save_path+'/'+self.experiment+'_'+self.position_string+'_'+'parameters_dict', 'wb') as handle:
            pickle.dump(parameters_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        return parameters_dict


    def run_parametric_segmentation(self):
        """
        This function uses the saved parameters for the specified frame intervals to store the segmented cell labels.
        """
        if os.path.isfile(self.save_path+'/'+self.experiment+'_'+self.position_string+'_'+'parameters_dict') == True:
            with open(self.save_path+'/'+self.experiment+'_'+self.position_string+'_'+'parameters_dict', 'rb') as handle:
                parameters_dict = pickle.load(handle)
        else:
            print('No parameters dictionary was detected. Specify the parameters for different frame intervals')

        watershed_dict = {}
        for frame_range in parameters_dict:
            seg_params = parameters_dict[frame_range]
            print('analyzing frame range:', frame_range)
            for fr in range(frame_range[0], frame_range[1]):
                print('segmenting cells in frame:', fr)
                watershed_dict[fr] = self.get_watershed_labels(fr, orientation=seg_params['orientation'], channel_end=seg_params['channel_end'], 
                            channel_start=seg_params['channel_start'], central_line_window=seg_params['central_line_window'], 
                            smoothing_factor=seg_params['smoothing_factor'], projection_width=seg_params['projection_width'],
                            std_threshold=seg_params['std_threshold'], adaptive_smooth=seg_params['adaptive_smooth'],
                            adaptive_window=seg_params['adaptive_window'], adaptive_offset=seg_params['adaptive_offset'], erosions=seg_params['erosions'], 
                            distance_threshold=seg_params['distance_threshold'], area_threshold=seg_params['area_threshold'],  
                            show_labels=False, show_threshold=False)
                
                fluor_threshold=seg_params['fluor_threshold']
                fluor_smoothing=seg_params['fluor_smoothing']
                orientation = seg_params['orientation']
                
                if orientation == 'up':
                    print('Inverting the labels to match the orientation of the channel')
                    watershed_dict[fr] = self.invert_labels(watershed_dict[fr])

                if fluor_threshold > 0:
                    watershed_dict[fr] = self.fluorescence_thresholding(watershed_dict[fr], self.fluor_images[fr], fluor_smoothing, fluor_threshold)
                    
        with open(self.save_path+'/'+self.experiment+'_'+self.position_string+'_'+'watershed_labels_dict', 'wb') as handle:
            pickle.dump(watershed_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        return watershed_dict


    def invert_labels(self, watershed_array):
        """
        This function is used to invert the array labels depending on the orientation of the microfluidics channel
        """
        max_label = watershed_array.max()
        sub_label = int(max_label + 1)
        curated_watershed_dict = watershed_array - sub_label
        curated_watershed_dict[curated_watershed_dict==-sub_label]=0
        
        return np.absolute(curated_watershed_dict)
        
        
    def run_segmentation(self, frames, orientation='down', channel_end=50, channel_start=15, 
                         central_line_window=3, smoothing_factor=9, projection_width=8,
                         std_threshold=0.8, adaptive_smooth=5, adaptive_window=5, adaptive_offset=-8,
                         erosions=3, distance_threshold=2, area_threshold=100,
                         show_labels=False, show_threshold=False, fluor_smoothing=2, fluor_threshold=230):
        """
        Use the phase intensity thresholding and watershedding to segment the cells.
        """
        watershed_dict = {}
        for fr in range(frames[0], frames[1]):
            if show_labels == False and show_threshold == False:
                print('segmenting cells in frame:', fr)
            watershed_dict[fr] = self.get_watershed_labels(fr, orientation, channel_end, channel_start, 
                                                        central_line_window, smoothing_factor, projection_width,
                                                        std_threshold,  adaptive_smooth, adaptive_window, adaptive_offset,
                                                        erosions, distance_threshold, 
                                                        area_threshold, show_labels, show_threshold)
            
            if orientation == 'up':
                print('Inverting the labels to match the orientation of the channel')
                watershed_dict[fr] = self.invert_labels(watershed_dict[fr])
                
            if fluor_threshold > 0 and self.fluor_images != 'none': 
                if len(self.fluor_images) == len(self.phase_images):
                    watershed_dict[fr] = self.fluorescence_thresholding(watershed_dict[fr], self.fluor_images[fr], fluor_smoothing, fluor_threshold)
                else:
                    print('There is not a single fluorescence image for each phase image. Fluorescence thresholding aborted...')
            else:
                print('There were no images loaded or the fluorescence threshold was set to zero. Fluorescence thresholding aborted...')
            
        with open(self.save_path+'/'+self.experiment+'_'+self.position_string+'_'+'watershed_labels_dict', 'wb') as handle:
            pickle.dump(watershed_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return watershed_dict    


    def generate_stacks(self, watershed_dict, length=300, height=3, edge=15):
        """
        This function is used to generate virtual stacks for manual curation of cell segmentation in microfluidics.
        """
        step = int(round((length/height),0))
        channel_length, channel_width = watershed_dict[0][:, edge:-edge].shape
        batch_dictionary = {}
        
        for batch in range(0, len(watershed_dict)-1, length):
            arrays_dict = {}
            for i in range(batch, batch+length, step):
                print(i)
                if i < len(watershed_dict):
                    first_array = watershed_dict[i][:, edge:-edge]
                    for fr in range(i+1, i+step):
                        if fr < len(watershed_dict):
                            first_array = np.hstack((first_array, watershed_dict[fr][:, edge:-edge]))
                        arrays_dict[i] = first_array
            
            first_stack = arrays_dict[batch]
            first_stack_dimensions = first_stack.shape
            for indx in arrays_dict:
                if indx > batch:
                    if arrays_dict[indx].shape == first_stack_dimensions:
                        first_stack = np.vstack((first_stack, arrays_dict[indx]))
                    else:
                        array_dimensions = arrays_dict[indx].shape
                        zeros_array = np.zeros((first_stack_dimensions[0], first_stack_dimensions[1]-array_dimensions[1]))
                        corrected_array = np.hstack((arrays_dict[indx], zeros_array))
                        first_stack = np.vstack((first_stack, corrected_array))
                        
        
            batch_dictionary[batch] = first_stack
        
        return batch_dictionary, (channel_length, channel_width)
    
    
    def assign_coordinates(self, watershed_dict, length=300, height=3, edge=15, cmap_='jet'):
        """
        This algorithm is used to manually assign splitting or merging coordinates and returns the curation dataframe.
        """
        batch_dictionary, (channel_length, channel_width) = self.generate_stacks(watershed_dict, length, height, edge)
        step = int(round((length/height),0))
        coordinates_dict = {}
        
        plt.figure(figsize=(500,1000))
        for bch in batch_dictionary:
            plt.imshow(batch_dictionary[bch], cmap=cmap_)
            plt.xticks([])
            plt.yticks([])
            coordinates = plt.ginput(n=0, timeout=0)
            try:
                coordinates_dict[bch] = coordinates
            except:
                print('no coordinates selected')
            plt.close()
        
        i = 0
        for bch in coordinates_dict:
            coords = coordinates_dict[bch]
            try:
                x_coords, y_coords = map(list, zip(*coords))
                x_coords = np.array(x_coords)
                y_coords = np.array(y_coords)
                frames = (bch + (x_coords/channel_width).astype(int)) + (y_coords/channel_length).astype(int)*step # - ((y_coords/channel_length).astype(int)+1)
                positions = np.subtract(y_coords, (y_coords/channel_length).astype(int)*channel_length)
                if i == 0:
                    curation_df = pd.DataFrame()
                    curation_df['frame'] = frames
                    curation_df['position'] = positions
                    i+=1
                elif i > 0:
                    temp_curation_df = pd.DataFrame()
                    temp_curation_df['frame'] = frames
                    temp_curation_df['position'] = positions
                    curation_df = pd.concat([curation_df, temp_curation_df])
                    
            except ValueError:
                print('No coordinates assigned for frames in batch:', bch)
        try:
            return curation_df
        except UnboundLocalError:
            print('No selections to curate.')
            
            
    def optimize_manual_split(self, splitting_coord, phase_image):
        """
        Optimizes the manually splitting coordinates to the point with the minimum 
        projection around the manually assigned coordinates
        """
        splitting_coord = int(round(splitting_coord,0)) # convert the splitting coord to integer to represent pixels
        image_cor = self.background_correction(phase_image)  # correct the image background
        upper_y, lower_y = splitting_coord+10, splitting_coord-10 # select a region around the splitting height
        crop_image = image_cor[lower_y:upper_y, :] # crop the image using the upper/lower bounds
        central_line_position = self.get_max_horizontal_projection_range(crop_image, window=5)
        projection = self.get_vertical_projection(crop_image, line_center=central_line_position, half_window=5, function='mean')
        trough = np.argmin(projection)
        
        return trough + lower_y
    
    
    def correction_watershed_labels(self, watershed_labels, phase_image, frame_df, orientation, merging_threshold=15):
        """
        Corrects the watershed labels using the curation dataframe
        """
        watershed_mask = watershed_labels>0
        merge_positions = []
        split_positions = []
        used_indexes = []
        
        if frame_df.shape[0] == 1:
            pos = frame_df.position.values[0]
            opt_pos = self.optimize_manual_split(pos, phase_image)
            split_positions.append(int(pos))
        elif frame_df.shape[0] > 1:
            for indx, row in frame_df.iterrows():
                if indx not in used_indexes:
                    used_indexes.append(indx)
                    
                    pos = row.position
                    pos_df = frame_df.copy()
                    pos_df = pos_df[pos_df.index!=indx]
                    pos_df['pos_dif'] = (pos_df['position'] - pos).abs()
                    
                    if pos_df.pos_dif.min() <= merging_threshold:
                        pos_2 = frame_df.loc[pos_df.pos_dif.idxmin(), 'position']
                        used_indexes.append(pos_df.pos_dif.idxmin())
                        merge_positions.append((int(round(pos,0)), int(round(pos_2,0))))
                    
                    elif pos_df.pos_dif.min() > merging_threshold:
                        opt_pos = self.optimize_manual_split(pos, phase_image)
                        split_positions.append(int(opt_pos))
            
        if len(split_positions) > 0:
            for pos in split_positions:
                print('Splitting at:', pos, 'px')
                if orientation == 'down':
                    watershed_labels[pos:,:]+=1 
                elif orientation == 'up':
                    watershed_labels[0:pos+1,:]+=1
                watershed_labels[~watershed_mask]=0
        
        if len(merge_positions) > 0:
            for pos in merge_positions:
                label_1 = np.max(watershed_labels[pos[0]])
                label_2 = np.max(watershed_labels[pos[1]])
                print('Merging labels',label_1, 'and',label_2,'at position:', int(np.mean(pos)), 'px')
                watershed_labels[watershed_labels==label_2]=label_1
        
        return watershed_labels
    
    
    def segmentation_curation(self, watershed_dict, length=300, height=3, edge=15, orientation='down', merging_threshold=15, cmap_='jet'):
        """
        Corrects all the segmented positions using the curation dataframe.
        """
        corrected_watershed_dict = {}
        
        curation_df = self.assign_coordinates(watershed_dict, length, height, edge, cmap_)
        for fr in watershed_dict:      
            phase_image = self.phase_images[fr].copy()
            watershed_labels = watershed_dict[fr].copy()
            frame_df = curation_df[curation_df.frame==fr]
            if frame_df.shape[0]  == 0:
#                    print('No curations for frame:', fr)
                corrected_watershed_dict[fr] = watershed_labels
            elif frame_df.shape[0] > 0:     
                print('Correcting frame:', fr)
                corrected_watershed_dict[fr] = self.correction_watershed_labels(watershed_labels, phase_image, frame_df, orientation, merging_threshold)
        return corrected_watershed_dict

            


class mother_machine_tracking(mother_machine_segmentation):
    """
    This class inherits all variables and attributes from the mother_machine_segmentation class.
    It includes all the functions needed to track segmented cells in microfluidics.
    """
    
    def __init__(self, deback_phase_path, deback_fluor_path, phase_interval, experiment, position, save_path):
        super().__init__(deback_phase_path,  deback_fluor_path, phase_interval, experiment, position, save_path)

    
    def remove_out_of_bounds(self, curated_watershed_dict, channel_bounds, frame_range):
        """
        Removes segmented lavels above and below a certain channel height.
        The channel orientation is NOT considered. 
        """
        for fr in range(frame_range[0], frame_range[1]+1):
            fr_labels = curated_watershed_dict[fr]
            if fr_labels[:channel_bounds[0]].shape[0]>0:
                upper_bad_labels = np.delete(np.unique(fr_labels[:channel_bounds[0]]),0)
                fr_labels[np.isin(fr_labels, upper_bad_labels)]=0
            if fr_labels[channel_bounds[1]+1:].shape[0]>0:
                lower_bad_labels =  np.delete(np.unique(fr_labels[channel_bounds[1]+1:]),0)
                fr_labels[np.isin(fr_labels, lower_bad_labels)]=0
            
            curated_watershed_dict[fr] = fr_labels
        
        return curated_watershed_dict
    
    
    def remove_bad_edges(self, curated_watershed_dict, channel_bounds_list, frame_range_list):
        """
        The channel_bounds_list and the frame_range_list should have the same dimensions.
        """
        for mov_seg in range(len(channel_bounds_list)):
            
            frame_range = frame_range_list[mov_seg]
            channel_bounds = channel_bounds_list[mov_seg]
            curated_watershed_dict = self.remove_out_of_bounds(curated_watershed_dict, channel_bounds, frame_range)
        
        return curated_watershed_dict
    
    
    def get_cell_IDs(self, curated_watershed_dict, min_area, frame_range):
        """
        This function generates cell IDs within a specified frame range
        for cell masks that are larger than a specified pixel area.
        """
        
        cell_id = []
        cell_area = []
        cell_centroid_xy = []
        cell_frame = []
        cell_label = []
        
        for fr in range(frame_range[0], frame_range[1]+1):
            for msk_lb in measure.regionprops(curated_watershed_dict[fr]):
                if msk_lb.area >= min_area:
                    cell_id.append(self.experiment+'_'+self.position_string+'_fr'+str(fr)+'_cell'+str(msk_lb.label))
                    cell_frame.append(fr)
                    cell_label.append(msk_lb.label)
                    cell_centroid_xy.append(msk_lb.centroid[::-1])
                    cell_area.append(msk_lb.area)
        
        cell_id_dataframe = pd.DataFrame()
        cell_id_dataframe['cell_id'] = cell_id
        cell_id_dataframe['frame'] = cell_frame
        cell_id_dataframe['label'] = cell_label
        cell_id_dataframe['centroid_xy'] = cell_centroid_xy
        cell_id_dataframe['area_px'] = cell_area
        cell_id_dataframe[['x','y']] = pd.DataFrame(cell_id_dataframe['centroid_xy'].tolist(), index=cell_id_dataframe.index)
        
        return cell_id_dataframe
    
            
    def link_cellpairs(self, curated_watershed_dict, min_area, frame_range, max_distance, area_ratio_range):
        """
        This function is used to link cells in pairs based on their distance and relative area.
        """
        cell_id_df = self.get_cell_IDs(curated_watershed_dict, min_area, frame_range)
        
        cell_linkage_dict = {}
        
        for index, row in cell_id_df.iterrows():
            
            fr = row.frame
            cl = row.cell_id
            ex = row.x
            wi = row.y
            ara = row.area_px
            
            if fr < frame_range[1]:
                fr_df = cell_id_df[cell_id_df.frame==fr+1]
                fr_df['distance'] = np.sqrt((fr_df.x-ex)**2+(fr_df.y-wi)**2)
                fr_df['area_ratio'] = fr_df.area_px / ara
                fr_df = fr_df[(fr_df.distance<=max_distance) & (fr_df.area_ratio.between(area_ratio_range[0], area_ratio_range[1]))]
                if fr_df.shape[0] == 0:
                    cell_linkage_dict[cl] = np.nan
                    print('cell:', cl, 'in frame:', fr, 'not lnked to another cell in the next frame')
                elif fr_df.shape[0] == 1:
                    cell_linkage_dict[cl] = fr_df.cell_id.values[0]
                    print('cell:', cl, 'in frame:', fr, 'linked to cell', cell_linkage_dict[cl], 'in frame', fr+1)
                elif fr_df.shape[0] > 1:
                    min_df = fr_df.loc[fr_df.distance.idxmin()]
                    cell_linkage_dict[cl] = min_df.cell_id
                    print('cell:', cl, 'in frame:', fr, 'linked to cell', cell_linkage_dict[cl], 'in frame', fr+1)
        
        cell_id_df['link_cell_id'] = cell_id_df.cell_id.map(cell_linkage_dict)
        
        return cell_id_df, cell_linkage_dict
    

    def track_cells(self, curated_watershed_dict, min_area, frame_range, max_distance, area_ratio_range):
        """
        Track the linked cells to generate trajectories and assign unique trajectory IDs.
        This function includes a double recursion to link cells into trajectories.
        """
        cell_id_df, cell_linkage_dict = self.link_cellpairs(curated_watershed_dict, min_area, 
                                                            frame_range, max_distance, area_ratio_range)
        
        def get_linked_cells(cell, temp_cell_id_df, cell_list):
            if type(cell)==str and cell!='nan':
#                print(cell)
                linked_cell = temp_cell_id_df[temp_cell_id_df.cell_id==cell].link_cell_id.values[0]
                cell_list.append(linked_cell)
                cell = linked_cell
                get_linked_cells(cell, cell_id_df, cell_list) # recursion 1

        # This function can be used for dial recursion bit the recusrion depth is reached this way
        def get_trajectories(temp_cell_id_df, traj, all_cells):
            cell_list = [temp_cell_id_df.cell_id.values[0]] #global variable for recursion 1
            cell = cell_list[0]
#            print(cell)
            get_linked_cells(cell, temp_cell_id_df, cell_list)
            cell_list.pop()
            trajectory_dict[traj] = cell_list
            all_cells += cell_list
            temp_cell_id_df = temp_cell_id_df[~temp_cell_id_df.cell_id.isin(cell_list)]
            if temp_cell_id_df.shape[0] > 0:
                traj+=1
                get_trajectories(temp_cell_id_df, traj, all_cells) # recursion 2
        
        # global variables
        trajectory_dict = {}
        all_cells = [] # global variable for all recursions
        try:
            temp_cell_id_df = cell_id_df.copy()
            traj=0
            get_trajectories(temp_cell_id_df, traj, all_cells)
        except RecursionError: # the the max depth of recursion is reached continue from where it stopped
            temp_cell_id_df = cell_id_df.copy()
            temp_cell_id_df = temp_cell_id_df[~temp_cell_id_df.cell_id.isin(all_cells)]
            traj = max(trajectory_dict.keys())+1
            get_trajectories(temp_cell_id_df, traj, all_cells)
        
        linkage_trajectory_dict = {v: k for k, values in trajectory_dict.items() for v in values}
        
        cell_id_df['cell_trajectory_int'] = cell_id_df.cell_id.map(linkage_trajectory_dict)
        cell_id_df['cell_trajectory_id'] = self.experiment+'_'+self.position_string+'_traj'+cell_id_df.cell_trajectory_int.astype(str)
        
        return cell_id_df
    
    
    def visualize_tracked_labels(self, background_corrected_phase, frame, watershed_labels, cell_id_df, tail_length=10, save_fig_path='none', cmap_='jet'):
        """
        This function can be used to visualize the watershedding labels for each specified frame.
        """
        fig = plt.figure(figsize=(3,10))
        gs = fig.add_gridspec(1, 2, hspace=0, wspace=0)
        (ax1, ax2) = gs.subplots(sharey='row')
        ax1.imshow(background_corrected_phase, cmap='gray')
        fr_df = cell_id_df[cell_id_df.frame==frame]
        trajectories_list = fr_df.cell_trajectory_int.unique().tolist()
        traj_df = cell_id_df[cell_id_df.cell_trajectory_int.isin(trajectories_list)]
        traj_df = traj_df[traj_df.frame<=frame]
        ax1.scatter(fr_df.x.values, fr_df.y.values, s=80, facecolors='none', edgecolors='r')
        if traj_df.shape[0] > 0:
            for traj in trajectories_list:
                single_traj_df = traj_df[traj_df.cell_trajectory_int==traj]
                if single_traj_df.shape[0] > tail_length:
                    single_traj_df = single_traj_df.tail(tail_length)
                ax1.plot(single_traj_df.x.values, single_traj_df.y.values)
        plt.text(x = -50, y =10, s = '0.066um/px', size=10, color='white', weight='bold')
        ax1.tick_params(axis='x',          # changes apply to the x-axis
                        which='both',      # both major and minor ticks are affected
                        bottom=False,      # ticks along the bottom edge are off
                        top=False,         # ticks along the top edge are off
                        labelbottom=False)
#            plt.imshow(markers)
        ax2.imshow(watershed_labels, cmap=cmap_)
#            plt.imshow(watershed_labels)
        plt.text(x = 0, y =10, s = 'time: '+str(frame*self.phase_interval)+'min', size=10, color='white', weight='bold')
        ax2.tick_params(axis='x',          # changes apply to the x-axis
                        which='both',      # both major and minor ticks are affected
                        bottom=False,      # ticks along the bottom edge are off
                        top=False,         # ticks along the top edge are off
                        labelbottom=False)
        if os.path.isdir(save_fig_path) == True:
            plt.savefig(save_fig_path+'/'+self.experiment+'_'+self.position_string+'_'+str(frame)+'_watershed_labels.jpeg')
        plt.show()


                
                
        
        
        
            
    
    
    


            
            