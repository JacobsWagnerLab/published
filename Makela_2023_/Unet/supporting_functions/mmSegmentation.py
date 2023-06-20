import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader
import numpy as np
from skimage import io, transform, filters, segmentation, measure, feature, morphology, util
from scipy import ndimage as ndi
import matplotlib.pyplot as plt

class mmSegmentation:
    
    """Segmentation loader script for Unet. Image can be segmented as a single
    image or splitted into square patches. Image is unfolded into smaller patches
    that are segmented and folded back into original shape. There is option for
    overlapping patches, overlap area is averaged between different patches and
    normalized to have same weight. Image batches are run in parellel and saved
    separately.
    
    Jarno Makela, 3/25/2021"""
    
    def __init__(self,net,threshold,sz,win_size,gpu):
        
        self.net = net
        self.use_cuda = torch.cuda.is_available()
        self.threshold = threshold
        self.sz = sz
        self.win_size = win_size
        self.gpu = gpu

    def segment(self, seg_loader: DataLoader):

        self.net.eval()
        
        # loops over individual images
        for sample in seg_loader:
            # image as tensor and path
            im, out_path = sample['im'], sample['path_to_save']
            
            # overlap creates problems at the moment, so disabled currently
            # -> kernel_size == stride
            kernel_size = self.win_size # patch size
            stride = kernel_size # defines how many actually are patches
            
            # if kernel size is different from image size -> make patches
            if kernel_size != 0 & kernel_size != self.sz[0]:
                
                # dimensions of tensor [batch, colors, width, height]
                B, C, W, H = im.shape
                
                # padding overlap mismatch doesn't work well at the moment
                
                # # number of pixels missing in each dimension:
                # pad_w = W % kernel_size
                # pad_h = H % kernel_size
                
                # # Padding the image with missing pixels:
                # im = F.pad(input=im, pad=(pad_w//2, pad_w-pad_w//2,
                #                                 pad_h//2, pad_h-pad_h//2), mode='constant', value=0)
                # # UPDATE the shape information to account for padding
                # B, C, H, W = im.shape
                
                # unfold tensor into patches and permute width and height positions
                # [B, C, no_patches_h, no_patches_w, kernel_size, kernel_size]
                patches = im.unfold(3, kernel_size, stride).unfold(2, kernel_size, stride).permute(0,1,2,3,5,4)
                
                # use zero tensor to save predictions because kernels overlap in patches and bleed to overlapping area
                patches_pred = torch.zeros(patches.size())
                
                # # loop over patches and segment each patch separately
                patches_dim = patches.size()
                for ii in range(patches_dim[2]):
                    for jj in range(patches_dim[3]):
                        patch = patches[:,:,ii,jj,:,:]
                        
                        # returns patch in either cpu or cuda memory
                        if self.gpu == 1:
                            patch = patch.cuda()
                        else:
                            patch = patch.cpu()  

                        # define not to store computations to free memory
                        with torch.no_grad():
                            # predictions for specific patch
                            pred_single = self.net(patch)
                        if len(pred_single)>1: #if deep supervision, use final prediction
                            pred_single = pred_single[-1]
                            
                        # save predictions to patches
                        patches_pred[:,:,ii,jj,:,:] = pred_single

                # reshape output to match input, see fold function
                patches_pred = patches_pred.contiguous().view(B, C, -1, kernel_size*kernel_size) # [B, C, no_patches_all, kernel_size*kernel_size]
                patches_pred = patches_pred.permute(0, 1, 3, 2) # [B, C, kernel_size*kernel_size, nb_patches_all]
                patches_pred = patches_pred.contiguous().view(B, C*kernel_size*kernel_size, -1) # [B, C*prod(kernel_size), L]
                pred = F.fold(patches_pred, output_size=(H, W), kernel_size=kernel_size, stride=stride)
                
                # control mask that mimics the original folding, in case overlapping regions
                control = F.fold(torch.ones_like(patches_pred), output_size=(H, W), kernel_size=kernel_size, stride=stride)
                
                # normalize the pred by the control mask
                pred = pred/control
                
                # sigmoid activation function to normalize values between 0 and 1
                pred = torch.sigmoid(pred)
                
                # empty cuda memory
                if self.gpu == 1:
                    torch.cuda.empty_cache()    
                
                # loop over image batch (tensor 1st dimension)
                pred_dim = pred.size()    
                for n in range(pred_dim[0]):
                    # single image
                    pred_single = pred[n]
                    
                    # reshape data into 2D numpy array, also ensure right size
                    pred_single = pred_single.to('cpu').detach().numpy().squeeze(0)
                    pred_single = pred_single[0:self.sz[0], 0:self.sz[1]]
                    
                    # treshold for cell areas
                    thresh = 0.5
                    pred_single = pred_single > thresh
                    
                    # remove small objects
                    pred_single = morphology.remove_small_objects(pred_single,100)
        
                    # save labels
                    pred_labels = measure.label(pred_single,connectivity=1).astype('uint16')
                    io.imsave(out_path[n],pred_labels,compress=6,check_contrast=False)
 
            else: # image segmentation without splitting
            
                # GPU or CPU
                if self.gpu == 1:
                    im = im.cuda()
                else:
                    im = im.cpu() 

                # get predictions from network (somehow multiple tensors)
                with torch.no_grad():
                    pred = self.net(im)
                if len(pred)>1: #if deep supervision, use final prediction
                    pred = pred[-1]
    
                # sigmoid activation function to normalize values between 0 and 1
                pred = torch.sigmoid(pred)
                
                # loop over image batch (tensor 1st dimension)
                pred_dim = pred.size()    
                for n in range(pred_dim[0]):
                    # single image
                    pred_single = pred[n]
                    
                    # reshape data into 2D numpy array, also ensure right size
                    pred_single = pred_single.to('cpu').detach().numpy().squeeze(0)
                    pred_single = pred_single[0:self.sz[0], 0:self.sz[1]]
                    
                    # treshold for cell areas
                    thresh = 0.5
                    pred_single = pred_single > thresh
                    
                    # remove small objects
                    pred_single = morphology.remove_small_objects(pred_single,100)
        
                    # save labels
                    pred_labels = measure.label(pred_single,connectivity=1).astype('uint16')
                    io.imsave(out_path[n],pred_labels,compress=6,check_contrast=False)
                        