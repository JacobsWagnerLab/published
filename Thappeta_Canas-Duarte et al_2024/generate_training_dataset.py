"""
Code for generating training dataset for Omnipose retraining.
Imports curated images and creates 3 replicates for each image, by rotating the original image by 90,180,270 degrees.
Saves the entire training dataset in one folder.

@author: Alessio Fragasso
March 12, 2024
"""

import skimage
from scipy import io
import os
import numpy as np
import matplotlib.pyplot as plt
import keyboard
import seaborn as sns
import torch
import time
import tifffile
import imageio.v3 as iio
import random
from scipy import ndimage
from PIL import Image as im
from matplotlib.colors import ListedColormap

''' Inputs:
 - folder_path = path to images to screen
 - save_folder = destination folder where to store the curated images
'''

images_list = os.listdir(folder_path)
images_list = sorted(images_list, key=lambda x: os.path.getmtime(os.path.join(folder_path, x)))

masks_list = [s for s in images_list if 'masks' in s]
phase_list = [s for s in images_list if 'masks' not in s]
number_list = []

for num in range(1, 10001):
    formatted_num = f"{num:05d}"  # Format with leading zeros for 5 digits
    number_list.append(formatted_num)

k = 0
for i in range(0, len(images_list)):
    img_name = images_list[i][:-4]
    rot_0 = iio.imread(folder_path+'/'+images_list[i])
    rot_1 = np.rot90(rot_0, k=1)
    rot_2 = np.rot90(rot_0, k=2)
    rot_3 = np.rot90(rot_0, k=3)
    
    if 'mask' in img_name:
        unique_integers = np.unique(rot_0)  # Replace with your data or class labels
        random_colors = np.random.rand(np.max(unique_integers), 3)
        color_list = [[1, 1, 1]]  # White color for zero
        color_list.extend(random_colors)
        random_colormap = ListedColormap(color_list)
        
        print(random_colormap)
    
        fig = plt.figure(figsize=(12, 6))
    
        fig.suptitle(img_name)
        ax1 = fig.add_subplot(1, 4, 1)
        plt.imshow(rot_0, cmap=random_colormap, interpolation='none')
        plt.title('0 deg')
    
        ax1 = fig.add_subplot(1, 4, 2)
        plt.imshow(rot_1, cmap=random_colormap, interpolation='none')
        plt.title('90 deg')
    
        ax1 = fig.add_subplot(1, 4, 3)
        plt.imshow(rot_2, cmap=random_colormap, interpolation='none')
        plt.title('180 deg')
    
        ax1 = fig.add_subplot(1, 4, 4)
        plt.imshow(rot_3, cmap=random_colormap, interpolation='none')
        plt.title('270 deg')
    
        plt.show()
    else:
        fig = plt.figure(figsize=(12, 6))
    
        fig.suptitle(img_name)
        ax1 = fig.add_subplot(1, 4, 1)
        plt.imshow(rot_0, cmap='gray')
        plt.title('0 deg')
    
        ax1 = fig.add_subplot(1, 4, 2)
        plt.imshow(rot_1, cmap='gray')
        plt.title('90 deg')
    
        ax1 = fig.add_subplot(1, 4, 3)
        plt.imshow(rot_2, cmap='gray')
        plt.title('180 deg')
    
        ax1 = fig.add_subplot(1, 4, 4)
        plt.imshow(rot_3, cmap='gray')
        plt.title('270 deg')
    
        plt.show()
        

    if 'masks' in img_name:
        file_name = 'train_' + number_list[k]+'_0' + '_masks'
        tifffile.imsave(save_path+'/'+file_name+'.tif', rot_0)
        file_name = 'train_' + number_list[k]+'_1' + '_masks'
        tifffile.imsave(save_path+'/'+file_name+'.tif', rot_1)
        file_name = 'train_' + number_list[k]+'_2' + '_masks'
        tifffile.imsave(save_path+'/'+file_name+'.tif', rot_2)
        file_name = 'train_' + number_list[k]+'_3' + '_masks'
        tifffile.imsave(save_path+'/'+file_name+'.tif', rot_3)
    else:
        k += 1
        file_name = 'train_' + number_list[k]+'_0'
        tifffile.imsave(save_path+'/'+file_name+'.tif', rot_0)
        file_name = 'train_' + number_list[k]+'_1'
        tifffile.imsave(save_path+'/'+file_name+'.tif', rot_1)
        file_name = 'train_' + number_list[k]+'_2'
        tifffile.imsave(save_path+'/'+file_name+'.tif', rot_2)
        file_name = 'train_' + number_list[k]+'_3'
        tifffile.imsave(save_path+'/'+file_name+'.tif', rot_3)


