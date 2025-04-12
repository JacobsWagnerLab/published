# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 16:40:55 2025

@author: fragasso
"""

import numpy as np
import pandas as pd




"""
Calculate Euclidean distances between consecutive coordinates
"""

def calculate_euclidean_distance(coords, pixel_size, aug, medial_axis_=[]):
    # Calculate the difference in x and y coordinates between consecutive pixels
    if len(medial_axis_)>0:   # if medial axis is provided, add start coordinate to get the correct absolute coordinates of the mean proj xx vector
        coords = np.vstack((medial_axis_[0],coords))
    dx = np.diff(coords[:, 0]) / aug
    dy = np.diff(coords[:, 1]) / aug

    # Calculate the Euclidean distance between consecutive pixels
    distances = pixel_size * np.sqrt(dx**2 + dy**2) 
    absolute_distances = np.insert(np.cumsum(distances), 0, 0)
    return distances[1:], absolute_distances[1:]

"""
Calculate cell length from the medial projection DataFrame for a single cell 
(or a single frame).
"""
def calculate_cell_length(cell_medial_proj_df, px_size, aug):
    # Stack the medial projection coordinates into an (N, 2) array.
    coords = np.vstack([
        cell_medial_proj_df['m_proj_x'].values,
        cell_medial_proj_df['m_proj_y'].values
    ]).T
    # Compute cumulative Euclidean distance.
    _, cumulative = calculate_euclidean_distance(coords, pixel_size=px_size, aug=aug)
    return cumulative[-1]


"""
Calculate Euclidean distances between consecutive coordinates and group them into
contiguous blocks of valid steps. A valid step is defined as one where the distance
does not exceed gap_factor times the median step distance.

This is used to ensure that in the case of multiple nucleoids, actual nucleoid length is not computed over the gaps in between the nucleoids
"""

def calculate_euclidean_distance_nuc(coords, pixel_size, aug, medial_axis_=[], gap_factor=1.5):

    # If a medial axis is provided, prepend its first coordinate.
    if len(medial_axis_) > 0:
        coords = np.vstack((medial_axis_[0], coords))
    
    # Calculate differences (scaled)
    dx = np.diff(coords[:, 0]) / aug
    dy = np.diff(coords[:, 1]) / aug
    
    # Compute Euclidean distances between consecutive points.
    distances = pixel_size * np.sqrt(dx**2 + dy**2)
    
    # Define the threshold for a valid step.
    median_distance = np.median(distances)
    threshold = gap_factor * median_distance
    valid = distances <= threshold
    
    # Split distances into contiguous blocks of valid steps.
    valid_blocks = []
    cumulative_blocks = []
    
    current_block = []
    current_cum = []
    cum = 0  # cumulative distance for the current block
    
    for i, dist in enumerate(distances):
        if valid[i]:
            cum += dist
            current_block.append(dist)
            current_cum.append(cum)
        else:
            # End of current valid block due to a gap.
            if current_block:
                valid_blocks.append(np.array(current_block))
                cumulative_blocks.append(np.array(current_cum))
                # Reset for the next block.
                current_block = []
                current_cum = []
                cum = 0
            # The gap distance is skipped.
    
    # Add the last block if it exists.
    if current_block:
        valid_blocks.append(np.array(current_block))
        cumulative_blocks.append(np.array(current_cum))
    
    return valid_blocks, cumulative_blocks
 """
 Calculate nucleoid length from the medial projection DataFrame after filtering
 for nucleoid pixels.
"""
def calculate_nucleoid_length(cell_medial_proj_df, nuc_coords_set, px_size, aug):
   
    df = cell_medial_proj_df.copy()
    # Mark pixels that are in the nucleoid mask.
    df['found'] = df.apply(lambda row: 1 if (row['x_px'], row['y_px']) in nuc_coords_set else 0, axis=1)
    nuc_df = df[df['found'] == 1]
    if nuc_df.empty:
        return 0.0
    nuc_coords = np.vstack([
        nuc_df['m_proj_x'].values,
        nuc_df['m_proj_y'].values
    ]).T
    _, cumulative = calculate_euclidean_distance_nuc(nuc_coords, px_size, aug, medial_axis_=[], gap_factor=3)
    return cumulative[-1]
