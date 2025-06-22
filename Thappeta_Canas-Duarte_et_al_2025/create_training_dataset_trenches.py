# -*- coding: utf-8 -*-
"""

Code for generating training dataset for Omnipose retraining of mother machine data.
Imports curated images and creates 2 replicates for each image, by flipping the image along the vertical axis.
Creates training data by forming tile images, with 5 trenches per image.
Saves the entire training dataset in one folder.


Created on Tue Aug 29 16:55:04 2023

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
import imageio.v3 as iio
import random
from scipy import ndimage
from PIL import Image as im


''' Inputs:
 - folder_path = path to images to screen
 - save_folder = destination folder where to store the curated images
'''
images_list = os.listdir(folder_path)
images_list = sorted(images_list, key=lambda x: os.path.getmtime(
    os.path.join(folder_path, x)))

masks_list = [s for s in images_list if 'masks' in s]
phase_list = [s for s in images_list if 'masks' not in s]

number_list = []

for num in range(1, 10001):
    formatted_num = f"{num:05d}"  # Format with leading zeros for 5 digits
    number_list.append(formatted_num)

k = 0
for i in range(0, len(images_list)):
    img_name = images_list[i][:-4]
    img_1 = iio.imread(folder_path+'/'+images_list[i])
    img_2 = np.fliplr(img_1)

    fig = plt.figure(figsize=(12, 6))

    fig.suptitle(img_name)
    ax1 = fig.add_subplot(1, 4, 1)
    plt.imshow(img_1, cmap='gray')
    plt.title('img')

    ax1 = fig.add_subplot(1, 4, 2)
    plt.imshow(img_2, cmap='gray')
    plt.title('flipped')

    fig.suptitle(img_name)
    ax1 = fig.add_subplot(1, 4, 1)
    plt.imshow(img_1, cmap='gray')
    plt.title('img')

    ax1 = fig.add_subplot(1, 4, 2)
    plt.imshow(img_2, cmap='gray')
    plt.title('flipped')

    plt.show()

    if 'masks' not in img_name:
        k += 1
        file_name = 'train_empty_' + number_list[k]+'_0' + '_masks'
        tifffile.imwrite(save_path+'/'+file_name+'.tif', img_1)
        file_name = 'train_empty_' + number_list[k]+'_1' + '_masks'
        tifffile.imwrite(save_path+'/'+file_name+'.tif', img_2)

    elif 'masks' in img_name:
        file_name = 'train_empty_' + number_list[k]+'_0'
        tifffile.imwrite(save_path+'/'+file_name+'.tif', img_1)
        file_name = 'train_empty_' + number_list[k]+'_1'
        tifffile.imwrite(save_path+'/'+file_name+'.tif', img_2)


'''
Create tile images that are 240x295 big, with 5 trenches per image
'''

new_folder_path = "path to saved images after first augmentation"
new_save_folder_path = "path to output training dataset folder"

images_list_aug = os.listdir(new_folder_path)

masks_list_aug = [s for s in images_list_aug if 'masks' in s]
phase_list_aug = [s for s in images_list_aug if 'masks' not in s]

size = len(masks_list_aug)
unique_vector = np.arange(size)
np.random.shuffle(unique_vector)

masks_list_shuffled = [masks_list_aug[i] for i in unique_vector]
phase_list_shuffled = [phase_list_aug[i] for i in unique_vector]

num_rows = 1
num_cols = 5

# Create a blank tile image of the desired size

img_w = 59
img_h = 240

tile_width = num_cols * img_w
tile_height = num_rows * img_h


i1 = 0
k = 0
while i1 < size:

    tile_image_mask = np.zeros((tile_height, tile_width))
    tile_image_phase = np.zeros((tile_height, tile_width))

    for row in range(num_rows):
        for col in range(num_cols):
            y_start = row * img_w
            y_end = y_start + img_h
            x_start = col * img_w
            x_end = x_start + img_w
            curr_mask = iio.imread(new_folder_path+'/'+masks_list_shuffled[i1])
            curr_phase = iio.imread(
                new_folder_path+'/'+phase_list_shuffled[i1])
            tile_image_mask[y_start:y_end, x_start:x_end] = curr_mask
            tile_image_phase[y_start:y_end, x_start:x_end] = curr_phase
            i1 += 1

    mask_name = 'train_' + number_list[k]+'_masks'
    tifffile.imwrite(new_save_folder_path+'/'+mask_name+'.tif', tile_image_mask)

    phase_name = 'train_' + number_list[k]
    tifffile.imwrite(new_save_folder_path+'/' +
                    phase_name+'.tif', tile_image_phase)

    k += 1

    fig = plt.figure(figsize=(12, 6))
    fig.suptitle(phase_name)
    ax1 = fig.add_subplot(1, 2, 1)
    plt.imshow(tile_image_phase, cmap='gray')
    plt.title('Tile phase')

    ax1 = fig.add_subplot(1, 2, 2)
    plt.imshow(tile_image_mask)
    plt.title('Tile mask')
    plt.show()
