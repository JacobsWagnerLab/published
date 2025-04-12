# -*- coding: utf-8 -*-
"""
Functions used in the code for 

Created on Tue Jun  6 18:19:10 2023

@author: fragasso
"""

import os
import pandas as pd
import numpy as np
import pickle
from pims import ND2Reader_SDK
import seaborn as sns
import matplotlib.pyplot as plt
import skimage
import math 
import cupy as cp
import cupyx.scipy.ndimage as ndi
from matplotlib.colors import LinearSegmentedColormap, Normalize
from cupyx.profiler import benchmark
from scipy.spatial.distance import pdist, squareform
from numpy.lib.stride_tricks import as_strided
from skimage.morphology import remove_small_objects, binary_erosion
from skimage.transform import resize
from skimage import filters, morphology, measure
from skimage.registration import phase_cross_correlation
from skimage.measure import label, regionprops
from skimage.feature import peak_local_max
from skimage.segmentation import watershed, expand_labels
from skimage.filters import threshold_otsu, threshold_local
from scipy import ndimage
from scipy.ndimage import zoom, label, center_of_mass, sum as ndi_sum
from scipy.signal import find_peaks, lfiltic, lfilter
from scipy.optimize import least_squares, curve_fit
from scipy.ndimage import gaussian_filter, gaussian_filter1d, binary_fill_holes, binary_opening, binary_dilation, binary_erosion
from scipy.interpolate import splprep, splev
from scipy.signal import lfilter, savgol_filter
import imageio.v3 as iio
import scipy.misc 
from PIL import Image as im

'''
Helper functions for creating list of xy positions, or formatting given number n into cell_id cell00000n
'''

def create_pos_list(start, end, double_digit = True):
    pos_list = []
    if double_digit:
        pos_list = [f'{i:02d}' for i in range(start, end + 1)]
    else:
        pos_list =  [str(i) for i in range(start, end + 1)]
    pos_list = ['xy'+pos for pos in pos_list]
    return pos_list


def format_number_as_string(number):
    formatted_number = f"cell{number:07d}"
    return formatted_number


'''
This function applies background subtraction to inpput img, using mask to remove cell mask regions (after dilation). Code adapted from (Papagiannakis et al.,bioRxiv, 2024).
It's modified to convert numpy arrays to cupy arrays to run on GPU.
'''

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


'''
This function inputs the cell_stack founf in SuperSegger cell file and returns quantities: 
A (cell area), rcm (center of mass), r_off (global coordinate of top left pixel), angle (angle between major axis and y=0), 
x and y coordinates of bounding box, edge_flag (1 if cell is too close to the image edge)
'''
def get_cell_coordinates(cell_stack_temp):   # given cell stack from supersegger, extract cell coordinates
    cell_coord = cell_stack_temp['coord'][0, 0][0, 0]
    rcm = cell_coord[8][0]
    A = cell_coord[0][0]
    r_off = cell_stack_temp['r_offset'][0, 0][0]
    angle = cell_coord[0][0][0]
    edge_flag = cell_stack_temp['edgeFlag'][0, 0][0, 0]
    BB = np.zeros(4, dtype=np.uint16)
    for k in range(4):
        BB[k] = cell_stack_temp['BB'][0][0][0][k]
    x1 = BB[0]-1
    y1 = BB[1]-1
    x2 = BB[0]+BB[2]
    y2 = BB[1]+BB[3]
    
    return A, rcm, r_off, angle, x1, x2, y1, y2, edge_flag


'''Returns cropped image of the mask and pads coordinates'''
def get_pads(cell_mask,pad_off=100):    
    b=cp.where(cell_mask>0)
    pads = cp.array([cp.min(b[1])-pad_off,cp.max(b[1])+pad_off,
                     cp.min(b[0])-pad_off,cp.max(b[0])+pad_off])
    return cell_mask[pads[2]:pads[3],pads[0]:pads[1]], pads


'''Medial axis functions'''

'''This function converts a 2D array of pixel intensities to a 3-column array with y, x, coordinates and pixel intensities'''
def convert_to_coordinates(intensity_array):
    rows, columns = intensity_array.shape
    pixels_3d = np.zeros((rows * columns, 3))
    # Fill the 3-column array with coordinates and pixel intensities
    pixels_3d[:, 0] = np.repeat(np.arange(rows), columns)
    pixels_3d[:, 1] = np.tile(np.arange(columns), rows)
    pixels_3d[:, 2] = intensity_array.flatten()
    return pixels_3d

'''This function calculates euclidean distances between pixel coordinates, and the cumulative distance, namely the length, along the set of points. Specify augmentation factor used if sub-pixel calculations were used, otherwise aug = 1.'''

def calculate_euclidean_distance(coords, pixel_size, aug = 1):
    # Calculate the difference in x and y coordinates between consecutive pixels
    dx = np.diff(coords[:, 0]) / aug
    dy = np.diff(coords[:, 1]) / aug
    # Calculate the Euclidean distance between consecutive pixels
    distances = pixel_size * np.sqrt(dx**2 + dy**2) 
    absolute_distances = np.insert(np.cumsum(distances), 0, 0)
    return distances[1:], absolute_distances[1:]


'''This function adds padding and rotates the image by specified angle. Uses cupy arrays'''
def rotateImage_cp(img, angle, rcm, r_off, aug, rotate=True):
    pivot = np.array([round(rcm[0]-r_off[0])*aug,round(rcm[1]-r_off[1])*aug])     # get relative center of mass (not global, but within the cropped region) 
    padX = [img.shape[1] - pivot[0], pivot[0]]
    padY = [img.shape[0] - pivot[1], pivot[1]]
    if any(s<0 for s in padX) or any(s<0 for s in padY):
        print('..Something went wrong with padding and rotation, medial_axis will not be computed')
        return [],[]
    imgP = cp.pad(img, [padX, padY], 'constant')
    if rotate:
        imgR = ndi.rotate(imgP, angle, reshape=False, order=0)
        return imgR, imgP
    else:
        return imgP

'''This function rotates adds padding the image by specified angle. Uses numpy arrays'''
def rotateImage(img, angle, rcm , r_off, aug):
    pivot = np.array([round(rcm[0]-r_off[0])*aug,round(rcm[1]-r_off[1])*aug])     # get relative center of mass (not global, but within the cropped region) 
    padX = [img.shape[1] - pivot[0], pivot[0]]
    padY = [img.shape[0] - pivot[1], pivot[1]]
    imgP = np.pad(img, [padX, padY], 'constant')
    imgR = ndimage.rotate(imgP, angle, reshape=False,order=0)
    return imgR[(padY[0]):(-padY[1]), (padX[0]) : (-padX[1])],imgR

'''
This function computes medial axis given the cropped cell mask:
    - cropped_mask: 2D binary array containing the cropped cell mask 
    - r_off: global coordinate of the top left pixel of the cropped image
    - rcm: center of mass of cell mask in global coordinates
    - angle: angle of the major axis relative to y=0
    - aug: augmentation factor for sub-pixel computation
    - extra_p: excess pixels that are extrapolated during fitting of the medial axis

Returns: 
    - medial_axis_f: Computed medial axis as numpy array of x and y coordinates (relative to ther cropped augmented mask)
    - mask_aug_cp: The cropped augmented mask (default is 20X upsampling) as a cupy 2D array
    - bad_cell_flag: boolean, if 1 means the computation failed and cell is considered a bad cell, e.g. when cell is too small
This function uses cupy arrays to work on GPU. To compute on CPU: convert cupy to numpy arrays
'''

def get_medial_axis_cp(cropped_mask, r_off, rcm, angle, show=True, aug=20, extra_p=50):
    cropped_mask_cp=cp.asarray(cropped_mask)
    mask_aug_cp = tile_array_cp(cropped_mask_cp, aug, aug)                              # upsample by aug factor
    cell_mask, cell_mask_pad = rotateImage_cp(mask_aug_cp, angle-10, rcm , r_off, aug)  # rotate image. Provides both rotated image and non-rotated, only padded image
    if len(cell_mask)==0:
        bad_cell_flag = True
        return [], [], bad_cell_flag
    cell_mask_filt = ndi.binary_opening(cell_mask, iterations=90,brute_force=True)  # smoothen edges. Needed for skeletonize step  
    if cp.size(cp.array(cp.where(cell_mask_filt>0)).T)<50:
        if show:
            print('Cell has a weird shape or it is too small')
        medial_axis_try_f=[] 
        bad_cell_flag = True
        return medial_axis_try_f, mask_aug_cp, bad_cell_flag
    else:
        bad_cell_flag = False
        skeleton = morphology.skeletonize(cell_mask_filt.get())                                # need to convert back to numpy for skeleton computation
        skeleton = cp.asarray(skeleton)
        skeleton_list = []
        for idx in cp.column_stack(cp.where(skeleton)):
            skeleton_list.append(idx)
        sk_x = cp.zeros(len(skeleton_list))
        sk_y = cp.zeros(len(skeleton_list))
        for i, idx in enumerate(skeleton_list):
            sk_x[i] = idx[1]
            sk_y[i] = idx[0]
    
        x_min = cp.min(cp.where(cell_mask > 0)[1])
        x_max = cp.max(cp.where(cell_mask > 0)[1])
    
        idx_medial = cp.where(sk_x == cp.round(cp.median(sk_x)))
        medial_center = cp.array([sk_x[idx_medial], sk_y[idx_medial]]).flatten()
    
        '''
        Find segments of the skeleton where to fit linear regression to find the poles
        Segments are defined as the first 3/8 of the skeleton for the left pole, and the last 3/8 for the right pole 
        '''
        thr_r = abs(cp.max(sk_x) + medial_center[0]) / 2
        thr_rr = abs(thr_r + medial_center[0]) / 2
        thr_l = abs(cp.min(sk_x) + medial_center[0]) / 2
        thr_ll = abs(thr_l + medial_center[0]) / 2
        idx_r = cp.where(sk_x > thr_rr)[0]
        idx_l = cp.where(sk_x < thr_ll)[0]
        fit_vec_r_x, fit_vec_r_y = sk_x[idx_r], sk_y[idx_r]
        fit_vec_l_x, fit_vec_l_y = sk_x[idx_l], sk_y[idx_l]
    
        '''
        Fit linear regressions to segments. 
        If segments have at least 20 points (if cell is large enough), then method 1 is used, 
        where cell poles are extrapolated from linear regressions to these segments and a third-order polynomial fit is performed using the found poles and skeleton
        '''
        if len(fit_vec_r_x)>20 and len(fit_vec_r_y)>20:
            if show:
                print('Use method 1 for medial axis')
            coefs_r = cp.polyfit(fit_vec_r_x, fit_vec_r_y, deg=1)
            xx_r = cp.linspace(medial_center[0], x_max, 1000)
            ffit_r = cp.polyval(coefs_r, xx_r)
            r_img = cp.zeros_like(cell_mask)
            r_img[ffit_r.astype(int), xx_r.astype(int)] = 1
            r_img_2 = r_img * cell_mask
            idx_r = cp.column_stack(cp.where(r_img_2 > 0))
            idx_pole_r = cp.where(idx_r[:, 1] == cp.max(idx_r[:, 1]))[0]
            pole_r = idx_r[idx_pole_r].flatten()
        
            coefs_l = cp.polyfit(fit_vec_l_x, fit_vec_l_y, deg=1)
            xx_l = cp.linspace(x_min, medial_center[0], 100)
            ffit_l = cp.polyval(coefs_l, xx_l)
            l_img = cp.zeros_like(cell_mask)
            l_img[ffit_l.astype(int), xx_l.astype(int)] = 1
            l_img_2 = l_img * cell_mask
            idx_l = cp.column_stack(cp.where(l_img_2 > 0))
            idx_pole_l = cp.where(idx_l[:, 1] == cp.min(idx_l[:, 1]))[0]
            pole_l = idx_l[idx_pole_l].flatten()
            
            medial_axis_fit_x = cp.append(cp.append(pole_l[1], sk_x), pole_r[1])
            medial_axis_fit_y = cp.append(cp.append(pole_l[0], sk_y), pole_r[0])
            weights = cp.zeros_like(medial_axis_fit_y) + 0.1
            weights[0] = weights[-1] = 10
            weights[1] = weights[-2] = 1
            weights[len(weights) // 2] = 100
    
            coefs_m = cp.polyfit(medial_axis_fit_x, medial_axis_fit_y, deg=3, w=weights)
            xx_m = cp.linspace(medial_axis_fit_x[0]-extra_p, medial_axis_fit_x[-1]+extra_p, 1000)
            ffit_m = cp.polyval(coefs_m, xx_m)
            medial_axis = cp.column_stack((xx_m, ffit_m))
            medial_axis_or =cp.array(rotate_coordinates(medial_axis,cell_mask.shape,(angle-10))).astype(int)

        else:    
            '''
            If segments have at < 20 points (if cell is small), then method 2 is used, 
            where the cell mask is rotated to be fully horinzontal (major axis parallel to y=0), poles are found by taking the median of y-coordinates of pixels with either minimum x-coordinates (left pole) or maxium x-coordinates (right pole).
            Fit is performed using a second-order polynomial function through skeleton points and the found poles.
            '''
            if show:
                print('Use method 2 for medial axis')
            cell_mask, cell_mask_pad = rotateImage_cp(mask_aug_cp, angle, rcm , r_off, aug)  # rotate image. Provides both rotated image and non-rotated, only padded image
            cell_mask_filt = ndi.binary_opening(cell_mask, iterations=90,brute_force=True)
            nonzero_coords= cp.argwhere(cell_mask_filt)
            if len(nonzero_coords)==0:
                if show:
                    print('Cell has a weird shape or it is too small')
                medial_axis_try_f=[] 
                bad_cell_flag = True
                return medial_axis_try_f, mask_aug_cp, bad_cell_flag     
            nonzero_coords= cp.argwhere(cell_mask_filt)
            coord_y = nonzero_coords[:,0]
            coord_x = nonzero_coords[:,1]
            min_x = cp.min(coord_x)
            max_x = cp.max(coord_x)
            
            idx_l = cp.array(cp.where(coord_x==min_x)).T
            idx_r = cp.array(cp.where(coord_x==max_x)).T
            
            pole_l_y=cp.median(coord_y[idx_l])
            pole_r_y=cp.median(coord_y[idx_r])
            
            pole_l = cp.array([pole_l_y,min_x])
            pole_r = cp.array([pole_r_y,max_x])
            
            medial_axis_fit_x = cp.append(cp.append(pole_l[1], sk_x), pole_r[1])
            medial_axis_fit_y = cp.append(cp.append(pole_l[0], sk_y), pole_r[0])
            weights = cp.zeros_like(medial_axis_fit_y) + 0.1
        
            coefs_m = cp.polyfit(medial_axis_fit_x, medial_axis_fit_y, deg=2)
            xx_m = cp.linspace(medial_axis_fit_x[0]-extra_p, medial_axis_fit_x[-1]+extra_p, 1000)
            ffit_m = cp.polyval(coefs_m, xx_m)
            medial_axis = cp.column_stack((xx_m, ffit_m))
            medial_axis_or =cp.array(rotate_coordinates(medial_axis,cell_mask.shape,(angle))).astype(int)        # rotate medial axies back to original image
    
    
        cell_mask_cropped, pads=get_pads(cell_mask_pad)
        medial_axis_x = medial_axis_or[:,1]-pads[2]
        medial_axis_y = medial_axis_or[:,0]-pads[0]
        medial_axis = cp.array([medial_axis_x,medial_axis_y]).T
        
        cell_mask_idx = cp.array(cp.where(mask_aug_cp>0)).T
        medial_axis_f=medial_axis[(medial_axis[:, None] == cell_mask_idx).all(-1).any(-1)]
        medial_axis_f=unique_rows_preserve_order(medial_axis_f.get())
    
        
        if show:
            plt.imshow(mask_aug_cp.get())
            plt.plot(medial_axis_f[:, 1], medial_axis_f[:, 0])
            plt.show()
    
    return medial_axis_f, mask_aug_cp, bad_cell_flag



'''This function maps pixels from 2D arrays (fluor_ch) to medial axis coordinates.
Returns dataframes:
   - medial_df_sorted that contains for each (x_px, y_px) pixel coordinates of the original cell mask (croped_mask), the corrsponding coordinates of the its projection on (minimum distance from) the medial axis (m_proj_x, m_proj_y)
   - mean_px_df that contains mean intensity projection of pixel intensities on the medial axis from the whole cell mask
   - medial_df_thr that is equal medial_df_sorted, except it filters out pixels that are more than 5-pixels away from the medial axis
   - mean_px_df_thr that is equal to mean_px_df, except it filters out pixels that are more than 5-pixels away from the medial axis
'''

def get_1D_proj_medial(fluor_ch, cropped_mask, medial_axis, aug=20, px_dist_thr = 5):   ## provide 1D projection on medial axis of cell pixels
    
    fluor_1_coord = convert_to_coordinates(fluor_ch)        # convert 2D image into array of x,y coord, and pix intensities
    cell_mask_coord = convert_to_coordinates(cropped_mask)      # convert 2D image into array of x,y coord, and pix intensities
    m_x = medial_axis[:, 1]
    m_y = medial_axis[:, 0]
    f1_coord_cropped = fluor_1_coord.copy()
    # select only pixels within cell mask
    f1_coord_cropped = f1_coord_cropped[np.where(cell_mask_coord[:, 2] > 0)]
    # expand coordinates to match the sub-pixel medial axis coordinates
    f1_coord_cropped[:, 0] *= aug
    # add offset to bring the coordinate to the pixel center in the expanded space
    f1_coord_cropped[:, 0] += round(aug/2)
    f1_coord_cropped[:, 1] *= aug
    f1_coord_cropped[:, 1] += round(aug/2)
    f1_x = f1_coord_cropped[:, 1]
    f1_y = f1_coord_cropped[:, 0]
    medial_df = pd.DataFrame()
    # list of medial coordinates where pixels of f1_coord are projected
    medial_proj_coord = []
    px_int = []                # pixel intensity
    px_dist = []               # distance from medial axis
    for i in range(0, len(f1_coord_cropped)):
        x_px = f1_x[i]
        y_px = f1_y[i]
        # for each pixel in f1_coord, calculate all distances to medial axis pixel coordinates
        distances = np.sqrt((x_px - m_x)**2 + (y_px - m_y)**2)
        # find minimum of such distances, and return its indices along medial axis array
        medial_proj_coord.append(np.argmin(distances))
        # append the pixel value of the pixel
        px_int.append(f1_coord_cropped[i, 2])
        # append pixel distance of the pixel from the closest medial axis pixel
        px_dist.append(np.min(distances))

    m_proj_coord = [tuple(row) for row in medial_axis[medial_proj_coord]]
    medial_df['x_px'] = f1_x
    medial_df['y_px'] = f1_y
    # coordinates on medial axis where fluor coordinates i are projected
    medial_df['m_proj_x'] = medial_axis[medial_proj_coord][:, 1]
    # coordinates on medial axis where fluor coordinates i are projected
    medial_df['m_proj_y'] = medial_axis[medial_proj_coord][:, 0]
    medial_df['m_proj_coord'] = m_proj_coord
    medial_df['px_dist'] = px_dist
    medial_df['px_int'] = px_int

    medial_axis_tup = [tuple(row) for row in medial_axis]
    medial_ordered = []
    i = 0
    idx = np.zeros(medial_df.shape[0])
    for coord in medial_axis_tup:
        if coord in medial_df['m_proj_coord'].tolist():
            match = np.where(coord == medial_df['m_proj_coord'])
            idx[i:(i+len(match[0]))] = match[0]
            medial_ordered.append([coord])
            i += len(match[0])

    medial_df_sorted = medial_df.loc[idx].reset_index()
    medial_df_sorted['OriginalIndex'] = range(len(medial_df_sorted))
    mean_df = medial_df_sorted.groupby('m_proj_coord').mean().reset_index()
    mean_df = mean_df.sort_values('OriginalIndex').reset_index(drop=True)
    mean_df.set_index(mean_df['m_proj_coord'],inplace=True)
    mean_px_df = mean_df['px_int']
    
    medial_df_thr = medial_df_sorted[medial_df_sorted['px_dist']<5*aug]
    mean_df_thr = medial_df_thr.groupby('m_proj_coord').mean().reset_index()
    mean_df_thr = mean_df_thr.sort_values('OriginalIndex').reset_index(drop=True)
    mean_df_thr.set_index(mean_df_thr['m_proj_coord'],inplace=True)
    mean_px_df_thr = mean_df_thr['px_int']
    return medial_df_sorted, mean_px_df, medial_df_thr, mean_px_df_thr
