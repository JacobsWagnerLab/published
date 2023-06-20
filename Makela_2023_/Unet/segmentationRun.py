from __future__ import print_function, division
import torch
from torch.utils.data import DataLoader
from torchvision import transforms
from supporting_functions.loadDataSets import mmDataSetSegment
import supporting_functions.customTransformations as ct
from supporting_functions import fileLoaders as fl
from supporting_functions.mmSegmentation import mmSegmentation
from supporting_functions.postProcess import par_watershed
import numpy as np
from skimage import io
import time
import warnings

# import unet
from networks.UnetPlusPlus2 import NestedUNet, UNet, ShallowNestedUNet

"""
A Deep Convolutional Neural network algorithm for cell segmentation. 
The U-net architecture allows pixel-based segmentation of high-level features 
from microscopy images. With a proper GPU it is also very fast, although memory 
intensive. Normal CPU based segmentation gives similar results but is slower.
Saves predictions as *_mask.tif to the same folder. Images can be split into
smaller parts if RAM/GPU memory is not enough for large images.
Additionally contains post-processing to clean up masks after segmentation.

Requires PyTorch, Numpy and Skimage

Set path to Unet folder that contains supporting_functions

Jarno Makela, 3/25/2021
"""

# folder path 
exp_dir = 'C:\\Work\\Projects\\Gene Dosage\\Data\\210322_JM169_TS\\210322_JM169_M9gluCAAT_38C_TS'

# file ending to images to be segmented
file_ending = 'c1.tif'

# unet model file
saved_model = 'Unet_11px_v1.pth'

# splitting image to square patches of win_size (e.g. 512); 0 no splitting
win_size = 0

# use GPU for segmentation if 1 (needs CUDA); 0 for CPU
use_gpu = 0

# number of images processed simultaneously; adjust according to available RAM/GPU memory
batch_size = 1

# postprocessing parameters
max_hole_fill = 50;     # smaller holes than value (in pixels) in masks are filled
max_area_remove = 100;  # smaller mask areas than value (in pixels) are removed
max_ws_length = 3;      # maximum length of a split between cells (in pixels)
h = 2;                  # minima transform value for watershed


###############################################
# Cell segmentation

warnings.filterwarnings("ignore")

# image file list in folder, define file ending to choose right channel
ims = sorted(fl.getFileList(exp_dir, first_only=False, file_ending=file_ending))

# tranformations for images
seg_transformations = transforms.Compose([transforms.ToPILImage(),
                                          ct.PILPadTo16(),
                                          transforms.ToTensor(),
                                          ct.SegmentationNormalize()])
# image set for pytorch data loader
ims_to_segment = mmDataSetSegment(phase_ims=ims,
                                  transform=seg_transformations)

# dataloader - define batch size and sampling
segmentation_loader = DataLoader(ims_to_segment, batch_size=batch_size, shuffle=False, num_workers=0)

# import unet for pytorch
net = ShallowNestedUNet(deep_supervision=True)
saved_net = torch.load(saved_model,map_location=torch.device('cpu'))
net.load_state_dict(saved_net['model_state_dict'])

# choose CPU or GPU (cuda)
if use_gpu == 1:
    net.cuda()
else:
    net.cpu()

# classifier
threshold = 0.5
sz = np.array(io.imread(ims[0]).shape).astype('int')
classifier = mmSegmentation(net, threshold, sz, win_size, use_gpu)

# segmentation
print('Starting segmentation.')
start = time.time()
classifier.segment(segmentation_loader)
end = time.time()
print(f"Segmentation took {round((end-start)/len(ims),2)} s per image: total {round(end-start,1)} s.")

###############################################
# Postprocessing to fill holes, remove too small cells, and watershed

print('Starting postprocessing.')
start = time.time()
ims_pred = sorted(fl.getFileList(exp_dir, first_only=False, file_ending='mask.tif'))
# loop over images with 'file_ending'
for im_pred in ims_pred:
    par_watershed(im_pred, h, max_ws_length, max_area_remove, max_hole_fill)

end = time.time()
print(f"Postprocessing took {round((end-start)/len(ims_pred),2)} s per image: total {round(end-start,1)} s.")
