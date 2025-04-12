# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 16:21:30 2025

@author: fragasso
"""

import numpy as np
import pandas as pd

'''This function converts a 2D array of pixel intensities to a 3-column array with y, x, coordinates and pixel intensities'''
def convert_to_coordinates(intensity_array):
    rows, columns = intensity_array.shape
    pixels_3d = np.zeros((rows * columns, 3))
    # Fill the 3-column array with coordinates and pixel intensities
    pixels_3d[:, 0] = np.repeat(np.arange(rows), columns)
    pixels_3d[:, 1] = np.tile(np.arange(columns), rows)
    pixels_3d[:, 2] = intensity_array.flatten()
    return pixels_3d

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