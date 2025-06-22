# -*- coding: utf-8 -*-
"""
Code for screening good from bad segmentation of mother machine data. 
Apply binary_fill_holes function from scipy.ndimage to attempt recovery of poorly segmented cells. 
Include empty trenches (no masks found) in the training dataset.

@author: fragasso
"""

from scipy.signal import lfiltic, lfilter
import skimage
from scipy import io
import os
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import keyboard
import seaborn as sns
import math
import torch
import time
import cupy as cp
import tifffile
from scipy.interpolate import splprep, splev, splrep, interp1d, UnivariateSpline, BSpline
from scipy.signal import lfilter, savgol_filter
from scipy.stats import spearmanr
from skimage import graph, morphology, measure
import imageio.v3 as iio
from scipy.ndimage import gaussian_filter, gaussian_filter1d, binary_fill_holes, binary_closing, binary_opening, binary_erosion
from scipy import ndimage
from PIL import Image as im
from skimage.morphology import skeletonize
from skimage.morphology import binary_closing
from skimage.filters import threshold_otsu, threshold_local
from skimage.segmentation import find_boundaries
from skimage import measure, morphology
from numpy.lib.stride_tricks import as_strided


''' Inputs:
 - folder_path = path to images to screen
 - save_folder = destination folder where to store the curated images
 - exp_name = name of the experiment
 - tw, offx = trench width and x-offset (px), to crop individual trenches for evaluation 
 - y_threshold = y-axis coordinates above which to remove masks, as there is overlap with the chip features
'''

mask_path = folder_path + '/masks'
masks_list = os.listdir(mask_path)
phase_list = os.listdir(phase_path)

for i in range(0, len(masks_list)):    
    print(masks_list[i])
    mask_temp = iio.imread(mask_path+'/'+masks_list[i])
    dim = mask_temp.shape
    mask_temp_bin = mask_temp > 0
    phase_temp = iio.imread(phase_path+'/'+phase_list[i])

    ### remove top masks
    labeled_array, num_features = ndimage.label(mask_temp_bin)
    center_of_mass = ndimage.center_of_mass(mask_temp_bin, labeled_array, range(1, num_features + 1))
    filtered_labels = [i for i, (y, x) in enumerate(center_of_mass, start=1) if y > y_threshold]
    mask_filtered = np.isin(labeled_array, filtered_labels)
    mask_filtered = mask_filtered[y1:y2,:]
    phase_temp = phase_temp[y1:y2,:]
    
    mask_2 = ndimage.binary_fill_holes(
        mask_filtered).astype(np.int16)

    fig = plt.figure(figsize=(12, 6))
    fig.suptitle('image '+str(i))
    ax1 = fig.add_subplot(2, 1, 1)
    plt.imshow(mask_2, cmap='gray')
    plt.title('Tile phase')

    ax1 = fig.add_subplot(2, 1, 2)
    plt.imshow(phase_temp, cmap='gray')
    plt.title('Tile mask')
    plt.show()

    n_trench = int(dim[1]/tw)
    for k in range(0,n_trench):
        mask_t = mask_2[:,(k*tw+offx):((k+1)*tw+offx)]
        if np.all(mask_t==0):    ## includes some good empty trenches in final training dataset
            trench_id='frame_'+str(i) + '_empty_trench_' + str(k)
            print(trench_id+' is an empty trench')
            continue      
        phase_t = phase_temp[:,(k*tw+offx):((k+1)*tw+offx)]
        labeled_t, _ = ndimage.label(mask_t)
        trench_id='frame_'+str(i) + '_trench_' + str(k)

    
        contours = measure.find_contours(mask_t, level=0.5)
        fig = plt.figure(figsize=(4, 4))
        fig.suptitle(trench_id)
        ax1 = fig.add_subplot(1, 3, 1)
        plt.imshow(labeled_t)
        plt.title('mask')

        ax1 = fig.add_subplot(1, 3, 2)
        plt.imshow(phase_t, cmap='gray')
        plt.title('phase')
        if contours != False:
            for contour in contours:
                ax1.plot(contour[:, 1], contour[:, 0],
                         linewidth=1, color='yellow')
        
        ax1 = fig.add_subplot(1, 3, 3)
        plt.imshow(phase_t, cmap='gray')
        plt.title('phase')
        
        plt.tight_layout()
        plt.show()
        
        while True:
            if keyboard.is_pressed('right'):
                path = save_folder + '/' + exp_name + '_' + pos+ '_'+ trench_id
                if save: 
                    tifffile.imsave(path+'.tif', mask_t)
                    tifffile.imsave(path+'_masks.tif', phase_t)
                print("Trench "+trench_id+" is moved to the good folder")
                break
            elif keyboard.is_pressed('left'):
                print("Trench "+trench_id+" is bad or empty")
                break
       

