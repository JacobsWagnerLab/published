from __future__ import division
import numpy as np
import torch
from skimage import io, transform, filters, exposure, util
from PIL import Image

###############################
#prepare custom transformations based on scimage to use in PyTorch. 
#Follow up: use Pytorch functions; requires PIL library
###############################

class RandomStretch(object):
    
    """Rescale the image in a sample to a given size.
    Args:
        scale_x, scale_y = the scale of transfomration in x and y dimensions.
        A new image of the same size will be produced. Always scales up, never down.
        
        Tranining transformation
    """

    def __init__(self):
        pass
        #assert isinstance(scale_x, (int, float, tuple)) and isinstance(scale_y, (int, float, tuple))
        #self.scale_x = scale_x
        #self.scale_y = scale_y

    def __call__(self,sample):
        
        scale_x = np.random.randint(65,240)/100
        scale_y = np.random.randint(65,240)/100
        
        image, mask = sample['im'], sample['mask']

        h, w = sample['im'].shape[:2]
        #find size of crop window
        
        h2, w2 = round(h * 1/scale_y), round(w*1/scale_x)
        #cut a crop region, then use transformation to automatically stretch it to fit        
        if h > h2:
            top = np.random.randint(0, h - h2)
        else:
            top = 0
        
        if w > w2:
            left = np.random.randint(0, w - w2)
        else:
            left = 0
            
        image = image[top: top + h2,
                      left: left + w2]
        
        mask = mask[top: top + h2,
                      left: left + w2]
        
        #scimage
        img =  transform.resize(image, (h, w))
        mask = transform.resize(mask, (h, w))
        
        if len(np.unique(mask)) > 1:
            mask_thresh = filters.threshold_otsu(mask)
            mask = mask > mask_thresh
        
        return {'im': img, 'mask': mask}

class RandomFlip90(object):
    """rotate image horizonstaly, adjust the size to match initial size
    
    Tranining transformation"""

    def __init__(self):
        pass
    
    def __call__(self,sample):
        
        image, mask = sample['im'], sample['mask']

        #50% chance to rotate
        rollDice = np.random.randint(0,2)
        
        if rollDice == 0:       
            image = np.rot90(image)
            mask = np.rot90(mask)

        return {'im' : image, 'mask' : mask}

class RandomCrop(object):
    """Crop randomly the image in a sample.

    Args:
        output_size (tuple or int): Desired output size. If int, square crop
            is made.
            
    Tranining transformation
    """

    def __init__(self,output_size):
        pass
        assert isinstance(output_size, (int, tuple))
        if isinstance(output_size, int):
            self.output_size = (output_size, output_size)
        else:
            assert len(output_size) == 2
            self.output_size = output_size
        
    def __call__(self, sample):
        #find out if the crop region is not bigger than the image
        
        image, mask = sample['im'], sample['mask']

        h, w = image.shape[:2]
        
        new_h, new_w = self.output_size
        
        if new_h < h:
            top = np.random.randint(0, h - new_h)
        else:
            top = 0
            
        if new_w < w:
            left = np.random.randint(0, w - new_w)
        else:
            left = 0
            
        image = image[top: top + new_h,
                      left: left + new_w]
        
        mask = mask[top: top + new_h,
                      left: left + new_w]


        return {'im' : image, 'mask' : mask}
    
class RandomRotate(object):
    """Rotation transformation
    
    Tranining transformation"""
    
    def __init__(self):
        pass
        #assert isinstance(angle, (int, float, tuple))
        #self.angle = angle
    
    def __call__(self, sample):

        angle = np.random.randint(-20,20)
        
        image, mask = sample['im'], sample['mask']
        
        #scimage
        image = transform.rotate(image,angle)
        mask = transform.rotate(mask,angle)
        
        #swithed off to see if still going to work
        #if len(np.unique(mask)) > 1:
        #    mask_thresh = filters.threshold_otsu(mask)
        #    mask = mask > mask_thresh
            
        return {'im' : image, 'mask' : mask}
 
class RemoveBackground(object):
    """change brighness of image
    
    Tranining transformation"""
    
    def __init__(self, kernel):
        self.kernel = kernel
    
    def __call__(self, sample):
        
        image = sample['im']
        
        #scimage
        image =  image - filters.gaussian(image, self.kernel)
        
        return {'im' : image, 'mask' : sample['mask']}
    
class CustomToTensor(object):
    """Rotation transformation
    
    Tranining transformation"""
    #shold be random rotation?
    
    def __init__(self):
        pass
        #assert isinstance(angle, (int, float, tuple))
        #self.angle = angle
    
    def __call__(self, sample):

        image, mask = sample['im'], sample['mask']
        dtype = torch.FloatTensor

        #scimage
        image = torch.from_numpy(image).unsqueeze(0).type(dtype)
        mask = torch.from_numpy(mask).unsqueeze(0).type(dtype)

        return {'im' : image, 'mask' : mask}
    
class CustomNormalize(object):
    """custom normalization of a tensor
    
    Tranining transformation"""
    
    def __init__(self):
        pass
    
    def __call__(self,sample):
        sample['im'] = (sample['im'] - torch.mean(sample['im'])) / torch.std(sample['im'])
        
        return sample
    
class CustomScale(object):
    """scale image by factor
    
    Tranining transformation"""
    
    def __init__(self, scale_factor):
        self.scale_factor = scale_factor
    
    def __call__(self, sample):
        image, mask = sample['im'], sample['mask']
        
        #scimage
        image =  transform.rescale(image, self.scale_factor, anti_aliasing=False)
        mask =  transform.rescale(mask, self.scale_factor, anti_aliasing=False)
        
        return {'im' : image, 'mask' : mask} 

class CustomBrightness(object):
    """change brightness of the image
    
    Tranining transformation"""
    
    def __init__(self):
        pass
    def __call__(self, sample):
        
        #initialise random change in the brightness
        gamma_factor = np.random.randint(40,180)/100
        
        image = sample['im']

        #scimage        
        sample['im'] = exposure.adjust_gamma(image,gamma=gamma_factor)
        return sample



#######################################################################################################
#                                   Here segmentation trasnformation   
#######################################################################################################

class RemoveBackgroundPIL(object):
    """remove bakcground from PIL image using skimage and numpy. Kernel should be used same as for the training. 
    
    Channel segmentation transformation"""
    
    def __init__(self, kernel):
        self.kernel = kernel
    
    def __call__(self, pic):
        im1 = np.asarray(pic)
        im2 = filters.gaussian(im1,self.kernel)
        im_out = im1 - im2
        im_out = Image.fromarray(im_out)
    
        return im_out

class PILPadTo16(object):
    """ pad image with 0 to match the size for Unet"""

    def __init__(self):
        pass
    
    def __call__(self,pic):
        im = np.asarray(pic)
        sz = im.shape
        pad_with = np.ceil(np.array(sz)/16)*16 - sz
        pad_with = pad_with.astype('int')
        im_pad = util.pad(im, pad_width=((0,pad_with[0]),(0,pad_with[1])),mode='constant')
        # im_pad = im_pad
        pic_out = Image.fromarray(im_pad)

        return pic_out

class PILRescale(object):
    """ rescale image by scale_factor"""

    def __init__(self,scale_factor):
        self.scale_factor = scale_factor
    
    def __call__(self,pic):
        
        w, h = pic.size
        
        w = np.round(w * self.scale_factor).astype('int')
        h = np.round(h * self.scale_factor).astype('int')

        pic_out = pic.resize((w,h))

        return pic_out
    
class SegmentationNormalize(object):
    """normalize tensor for segmentation script"""
    """ There is normalization problem if there are too few cells in the image. 
    This script forces a minimum std and prevents segmenting background"""

    def __init__(self):
        pass
    
    def __call__(self,tens):
        # normalization without any adjustments
        # tens = (tens - tens.mean()) / tens.std()
        
        # this normalizes the phase contrast data according to the training data
        # subtracts mean and estimates std by average of X dimmest pixels
        # assumes that dimmest pixels area from cell area
        
        # norm_cell_area is an estimate of normalized intensity of cell pixels in training data
        # no_pixels defines how many pixels are included into calculation of cell area brightness
        no_pixels = 1000
        norm_cell_area = -3;
        
        # subtract the mean pixel intensity
        tens = (tens - tens.mean())
        
        # convert tensor into array and calculate mean of X dimmest pixels 
        array = tens.numpy()
        vector = np.reshape(array,-1) # into vector
        vector = vector[vector != 0] # remove zeros
        vector_sorted = np.sort(vector)
        # calculate what value should be for X dimmest pixels to be at norm_cell_area
        value = np.mean(vector_sorted[1:no_pixels]) / norm_cell_area
               
        # divide by value (instead of std) to normalize the data
        tens = tens / value       
            
        return tens
