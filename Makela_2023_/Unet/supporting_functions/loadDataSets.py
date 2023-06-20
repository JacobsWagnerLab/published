import torch
from torch.utils.data import Dataset
import numpy as np
from skimage import io, transform, filters

class mmDataSetTrain(Dataset):
 
    """Dataset loader for training images.
    args: 
        phase_ims = list of phase contrast directories as full path: path/to/image/im.tiff
        mask_ims = list of masks to train segmentation as full path: path/to/mask/mask.tiff
        important! the list must be in the same order."""
    
    def __init__(self, phase_ims, mask_ims, transform=None):
        self.phase_ims = phase_ims
        self.mask_ims = mask_ims
        self.transform = transform
            
    def __len__(self):
        return len(self.phase_ims)
    
    def __getitem__(self, idx):
        
        #i dont know what it does but it was there in tutorial
        if torch.is_tensor(idx):
            idx = idx.tolist()
            
        #send to float32 type, PyTorch supposedly likes it
        phase_im = io.imread(self.phase_ims[idx],as_gray=True).astype('float32')
        
        mask_im = io.imread(self.mask_ims[idx],as_gray=True).astype('bool')
        
        sample = {'im' : phase_im, 'mask' : mask_im}
        
        if self.transform:
            sample = self.transform(sample)
            
        return sample
        
class mmDataSetSegment(Dataset):
    
    """Dataset loader for image segmentation.
    args: 
        phase_ims = list of phase contrast directories as full path: path/to/image/im.tiff
    returns normalised image and a path to save image to"""
    
    def __init__(self, phase_ims, transform=None):

        self.phase_ims = phase_ims
        self.transform = transform
            
    def __len__(self):
        return len(self.phase_ims)
    
    def __getitem__(self, idx):
        
        #i dont know what it does but it was there in tutorial
        if torch.is_tensor(idx):
            idx = idx.tolist()
            
        #send to float32 type, PyTorch supposedly likes it
        phase_im = io.imread(self.phase_ims[idx], as_gray=True).astype('float32')
        
        path_to_save = self.phase_ims[idx].replace('.tif', '_mask.tif')
        
        if self.transform:
            phase_im = self.transform(phase_im)
            
        return {'im' : phase_im, 'path_to_save' : path_to_save}