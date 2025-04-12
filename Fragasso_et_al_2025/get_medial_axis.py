# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 16:12:59 2025

@author: fragasso
"""


import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
from scipy import ndimage
from skimage import morphology



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

'''This function adds padding and rotates the image by specified angle. Uses numpy arrays'''
def rotateImage(img, angle, rcm , r_off, aug):
    pivot = np.array([round(rcm[0]-r_off[0])*aug,round(rcm[1]-r_off[1])*aug])     # get relative center of mass (not global, but within the cropped region) 
    padX = [img.shape[1] - pivot[0], pivot[0]]
    padY = [img.shape[0] - pivot[1], pivot[1]]
        
    imgP = np.pad(img, [padX, padY], 'constant')
    imgR = ndimage.rotate(imgP, angle, reshape=False,order=0)
    return imgR[(padY[0]):(-padY[1]), (padX[0]) : (-padX[1])],imgR


'''This function rotates adds padding the image by specified angle. Uses numpy arrays'''
def rotateImage(img, angle, rcm , r_off, aug):
    pivot = np.array([round(rcm[0]-r_off[0])*aug,round(rcm[1]-r_off[1])*aug])     # get relative center of mass (not global, but within the cropped region) 
    padX = [img.shape[1] - pivot[0], pivot[0]]
    padY = [img.shape[0] - pivot[1], pivot[1]]
    imgP = np.pad(img, [padX, padY], 'constant')
    imgR = ndimage.rotate(imgP, angle, reshape=False,order=0)
    return imgR[(padY[0]):(-padY[1]), (padX[0]) : (-padX[1])],imgR

'''Eliminate repeated rows in 2D array while keeping the original order of the array'''
def unique_rows_preserve_order(arr):
    _, idx = np.unique(arr, axis=0, return_index=True)
    unique_arr = arr[np.sort(idx)]
    return unique_arr

'''Upsamples cupy array a by b0 and b1'''
def tile_array_cp(a, b0, b1):
    r, c = a.shape                                    # number of rows/columns
    x = cp.repeat(a, b0, axis=0)                      # repeat rows
    x = cp.repeat(x, b1, axis=1)                      # repeat columns
    return x

'''Returns cropped image of the mask and pads coordinates'''
def get_pads(cell_mask,pad_off=100):    
    b=cp.where(cell_mask>0)
    pads = cp.array([cp.min(b[1])-pad_off,cp.max(b[1])+pad_off,
                     cp.min(b[0])-pad_off,cp.max(b[0])+pad_off])
    return cell_mask[pads[2]:pads[3],pads[0]:pads[1]], pads

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

