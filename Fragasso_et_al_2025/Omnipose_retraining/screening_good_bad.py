"""
Code for screening good from bad segmentation. 
Apply binary_fill_holes function from scipy.ndimage to attempt recovery of poorly segmented cells. 
Skips if masks are not found in the current image.

@author: Alessio Fragasso
March 12, 2024
"""

from scipy import io
import os
import numpy as np
import matplotlib.pyplot as plt
import keyboard
import seaborn as sns
from matplotlib.colors import ListedColormap
import tifffile
from skimage import graph, morphology, measure
import imageio.v3 as iio
from scipy.ndimage import gaussian_filter, gaussian_filter1d, binary_fill_holes, binary_closing, binary_opening, binary_erosion
from scipy import ndimage

''' Inputs:
 - folder_path = path to images to screen
 - save_folder = destination folder where to store the curated images
 - exp_name = name of the experiment
 - y1, y2 = dimensions of final images
'''
pos_list=os.listdir(folder_path) 
for pos in pos_list:
    if pos in bad_pos:
        continue
    exp_path = folder_path + '/' +pos
    mask_path = exp_path + '/masks'
    phase_path = exp_path + '/phase'
    masks_list = os.listdir(mask_path)
    phase_list = os.listdir(phase_path)

    for i in range(0, len(masks_list)):    # scans through  time points at fixed position
        print(masks_list[i])
        mask_temp = iio.imread(mask_path+'/'+masks_list[i])
        print(np.max(mask_temp))
        print(mask_temp.shape)
        unique_integers = np.unique(mask_temp)  # Replace with your data or class labels
        random_colors = np.random.rand(len(unique_integers) - 1, 3)
        color_list = [[1, 1, 1]]  # White color for zero
        color_list.extend(random_colors)
        random_colormap = ListedColormap(color_list)
        
        dim = mask_temp.shape
        mask_temp_bin = mask_temp > 0
        phase_temp = iio.imread(phase_path+'/'+phase_list[i])
    
        nx_box = int(dim[1]/y1)
        ny_box = int(dim[0]/y2)
        for j in range(nx_box):
            for k in range(ny_box):
                mask_cropped = mask_temp[y1*j:(j+1)*y1,y2*k:(k+1)*y2]
                phase_cropped = phase_temp[y1*j:(j+1)*y1,y2*k:(k+1)*y2]
                
                if np.all(mask_cropped==0):
                    print('empty image')
                    continue   
                
                cell_id = exp_name + '_' + pos+ '_' + str(i)+'_'+str(j)+str(k)
    
                contours = measure.find_contours(mask_cropped, level=0.5)
                
                fig = plt.figure(figsize=(12, 6))
                fig.suptitle('cell '+cell_id)
                ax1 = fig.add_subplot(1, 3, 1)
                plt.imshow(mask_cropped, cmap=random_colormap, interpolation='none')
                plt.title('Mask')            
            
                ax1 = fig.add_subplot(1, 3, 2)
                plt.imshow(phase_cropped, cmap='gray')
                plt.title('Phase')
                if contours != False:
                    for contour in contours:
                        ax1.plot(contour[:, 1], contour[:, 0],
                                 linewidth=1, color='yellow')
              
                
                ax1 = fig.add_subplot(1, 3, 3)
                plt.imshow(phase_cropped, cmap='gray')
                plt.title('Phase')
                
                plt.tight_layout()
                plt.show()

                ''' 
                Manual screening. Press 'Right' if image is good. Press 'Left' if image is bad and should be discarded. Press 'Up' if you want to apply a fill_holes function an reassess.
                '''
                while True:
                    if keyboard.is_pressed('right'):
                        path = save_folder + '/' + cell_id
                        tifffile.imsave(path+'.tif', phase_cropped)
                        tifffile.imsave(path+'_masks.tif', mask_cropped)
                        print("Cell "+cell_id+" moved to 'good' folder.")
                        break
                   
                    elif keyboard.is_pressed('left'):
                        print("Cell "+cell_id+" is bad")
                        break

                    elif keyboard.is_pressed('up'):   ##apply fill holes
                        print('Try fill holes')
                        mask_filled = ndimage.binary_fill_holes(
                            mask_cropped).astype(np.int16)
                        mask_labelled, _ = ndimage.label(mask_filled)
                        
                        contours = measure.find_contours(mask_labelled, level=0.5)
                        fig = plt.figure(figsize=(12, 6))
                        fig.suptitle('image '+str(i))
                        ax1 = fig.add_subplot(1, 3, 1)
                        plt.imshow(mask_labelled, cmap=random_colormap, interpolation='none')
                        plt.title('Mask')                    
                    
                        ax1 = fig.add_subplot(1, 3, 2)
                        plt.imshow(phase_cropped, cmap='gray')
                        plt.title('Phase')
                        if contours != False:
                            for contour in contours:
                                ax1.plot(contour[:, 1], contour[:, 0],
                                          linewidth=1, color='yellow')
                      
                        
                        ax1 = fig.add_subplot(1, 3, 3)
                        plt.imshow(phase_cropped, cmap='gray')
                        plt.title('Phase')
                        
                        plt.tight_layout()
                        plt.show()
                  
    
    