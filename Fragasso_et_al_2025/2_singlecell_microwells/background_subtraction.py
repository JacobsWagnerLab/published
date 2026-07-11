# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 14:47:43 2025

@author: fragasso
"""
import cupy as cp
import numpy as np
import cupy as cp
import cupyx.scipy.ndimage as ndi
import matplotlib.pyplot as plt
from analysis_functions_library import get_median_back_img

'''
This function applies background subtraction to inpput img, using mask to remove cell mask regions (after dilation). Code adapted from (Papagiannakis et al.,eLife, 2025).
It's modified to convert numpy arrays to cupy arrays to run on GPU.
'''


def get_median_back_img(img, img_cell_free, box_size, show_):
    img_dim = img.shape
    n_row = int(img_dim[0]/box_size)   # round to smallest integer
    n_col = int(img_dim[1]/box_size)
    # dim_new = n_row*box_size-sigma_, n_col*box_size-sigma_
    median_back_img = cp.zeros(img_dim)
    expand_window = False
    for i in range(0,n_row):
        for j in range(0,n_col):
            box_temp = img_cell_free[i*box_size:(i+1)*box_size,j*box_size:(j+1)*box_size]
            if np.any(box_temp>0)==False:
                if show_:
                    print('Using larger rolling window, too high cell density')
                expand_window = True
                return median_back_img, expand_window
            else:    
                median_back_img[i*box_size:(i+1)*box_size,j*box_size:(j+1)*box_size]=cp.median(box_temp[box_temp>0])
    if show_:  print('Found correct box_size = ' + str(box_size))
    return median_back_img, expand_window

def local_bkg_sub_cp(img, mask, pos, exp, frame='',dilations=15, box_size = 128, sigma_=60, show=False):
    binary_mask = cp.asarray(mask>0)
    img_cp = cp.asarray(img)
    img_cell_free = img_cp.copy()
    mask_dil = ndi.binary_dilation(binary_mask,iterations=dilations,brute_force=True)
    img_cell_free[cp.where(mask_dil==1)]=0  ## set to 0 values of fluorescence image within dilated masks
    img_dim = img.shape
    n_row = int(img_dim[0]/box_size)   # round to smallest integer
    n_col = int(img_dim[1]/box_size)
    dim_new = n_row*box_size-sigma_, n_col*box_size-sigma_
    expand_window = True
    bs = box_size
    while expand_window: 
        n_row = int(img_dim[0]/bs)   # round to smallest integer
        n_col = int(img_dim[1]/bs)
        dim_new = n_row*bs-sigma_, n_col*bs-sigma_
        median_back_temp, expand_window  = get_median_back_img(img, img_cell_free, bs,  show)
        bs = bs*2
        
    if expand_window:
        print('Cannot background subtract, too high cell density, try reducing mask dilation')
        return [], [], []
    
    median_back_img = median_back_temp
    
    back_img = ndi.gaussian_filter(median_back_img, sigma = sigma_)
    back_sub = img_cp - back_img
    
    if show:
        fig=plt.figure(figsize=(12, 6))
        fig.suptitle('Background subtraction at frame '+frame+ ', at pos ' + pos +', ' + exp)
        
        plt.subplot(3, 2, 1)
        plt.imshow(img_cp[0:dim_new[0],0:dim_new[1]].get())
        plt.title('Original Image')
        
        plt.subplot(3, 2, 2)
        plt.imshow(mask_dil[0:dim_new[0],0:dim_new[1]].get())
        plt.title('Dilated masks')
        
        plt.subplot(3, 2, 3)
        plt.imshow(img_cell_free[0:dim_new[0],0:dim_new[1]].get())
        plt.title('Cell free image')
        
        plt.subplot(3, 2, 4)
        plt.imshow(median_back_img[0:dim_new[0],0:dim_new[1]].get())
        plt.title('Local median Bkg')
        
        plt.subplot(3, 2, 5)
        plt.imshow(back_img[0:dim_new[0],0:dim_new[1]].get())
        plt.title('Local smooth Bkg')
        
        plt.subplot(3, 2, 6)
        plt.imshow(back_sub[0:dim_new[0],0:dim_new[1]].get())
        plt.title('Bkg sub img')
        plt.tight_layout()
        plt.show()
    
    return back_sub[0:dim_new[0],0:dim_new[1]].get(), back_img[0:dim_new[0],0:dim_new[1]].get(), median_back_img[0:dim_new[0],0:dim_new[1]].get()

