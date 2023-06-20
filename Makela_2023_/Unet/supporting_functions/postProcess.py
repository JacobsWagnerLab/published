import numpy as np
from skimage import io, transform, filters, segmentation, util, morphology, measure, feature
from scipy import ndimage as ndi
from matplotlib import pyplot as plt
import time
from supporting_functions import fileLoaders as fl

def par_watershed(im_path, h, max_ws, max_area, max_hole):
    
    """
    Script to post process mask images after segmentation. Fills holes in the
    cells, uses watershed to split joined cells and removes too small mask areas.
    Watershed function has additionally parameter to set a maximum length for
    watershed line (longer lines are not applied). 
    
    Note! It overwrites the original mask images
    
    Watershed based on this: 
    https://scikit-image.org/docs/dev/auto_examples/segmentation/plot_watershed.html
 
    Jarno Makela, 3/23/2021
    """
    # load image
    label_im = io.imread(im_path)
    
    # remove holes in the cells and convert to binary image
    bin_im = morphology.remove_small_holes(label_im > 0, max_hole)

    # distance transform
    distance = ndi.distance_transform_edt(bin_im)
    
    # define connectivity kernel for 2d image
    # connectivity size impacts some artifacts.
    conn = np.ones((5,5))

    # use morphological reconstruction, similar to the matlab example
    distance2 = morphology.reconstruction(distance-h,distance)
    local_maxi = feature.peak_local_max(distance2, indices=False, footprint=conn)
    markers = measure.label(local_maxi)
    watershed_labels = segmentation.watershed(-distance, markers, mask=bin_im, watershed_line=True)
    
    # find watershed lines (ridges) and give them unique labels
    watershed_labels[watershed_labels==0]=10e6
    ridges = measure.label((bin_im*watershed_labels==10e6))
    unique_areas = np.unique(ridges)
    
    # loop through watershed lines and apply them if area is smaller than ws_area
    for x in unique_areas[unique_areas > 0]:
        area = sum(sum(ridges==x))
        # remove regions small enough
        if area <= max_ws:
            bin_im[ridges==x]=0
 
    # remove small objects after watershedding
    bin_im = morphology.remove_small_objects(bin_im, max_area)
    labels = measure.label(bin_im,connectivity=1)
    
    # save labels
    io.imsave(im_path, labels.astype('uint16'), compress=6, check_contrast=False)
