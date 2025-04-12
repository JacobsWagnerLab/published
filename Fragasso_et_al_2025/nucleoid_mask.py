# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 14:45:18 2025

@author: fragasso
"""
import numpy as np
from scipy import ndimage
from skimage.filters import threshold_otsu
'''
This function is used to compute the nucleoid mask
'''
def get_fluor_mask(fluor_image,parameters=[2.5,0.95]):
    """
    Parameters are 
    [0] sigma of the gaussian filter
    [1] offset
    """
    log = -ndimage.gaussian_laplace(fluor_image, sigma=parameters[0])
    log_norm = (log- np.min(log)) / (np.max(log) - np.min(log))
    threshold_value = threshold_otsu(log_norm)
    mask = log_norm > threshold_value*parameters[1]
    return mask
