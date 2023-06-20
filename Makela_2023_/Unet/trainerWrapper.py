from __future__ import print_function, division
import os
import torch
import re
import time
from skimage import io, transform, filters
import numpy as np
import matplotlib.pyplot as plt
from torch.utils.data import Dataset, DataLoader
from torchvision import transforms, utils
from supporting_functions.loadDataSets import mmDataSetTrain
import supporting_functions.customTransformations as ct
from supporting_functions.fileLoaders import getFileList
from supporting_functions.mmClassifier import mmClassifier

#import unets
from networks.PraNetUnet import PraNetUnet, ShallowPraNetUnet
from networks.UnetPlusPlus2 import NestedUNet, UNet, ShallowNestedUNet

# Ignore warnings
import warnings
warnings.filterwarnings("ignore")

def TrainWholeNet(save_name):

    # function wraps the training parameters for Unet training usign PyTorch. 
    # The data for the tranining has to be stored in the folders with phase / brightfield images
    # and a folder with ground truth data (masks of cells). 
    #    
    # images for training for channels     
    phase_dir = 'C:\\Work\\Software\\PyModules\\KubaUnet\\Example training data\\exp3_agarosePad_phase'
    phase_ims = sorted(getFileList(phase_dir))
    print(len(phase_ims))
    
    mask_dir = 'C:\\Work\\Software\\PyModules\\KubaUnet\\Example training data\\exp3_agarosePad_gt'
    mask_ims = sorted(getFileList(mask_dir))
    print(len(mask_ims))

    #set up custom transformations

    #crop window can be balanced for good batch size, depending on amount of memory. 
    #Important: Unet will not segment images smaller that this, only same or bigger
    crop_window = 512

    composed_train = transforms.Compose([ct.RandomCrop(crop_window),
                                         ct.RandomFlip90(),
                                         ct.RandomStretch(),
                                         ct.RandomRotate(),
                                         ct.CustomToTensor(),
                                         ct.CustomNormalize()])
    
    composed_val = transforms.Compose([ct.RandomCrop(crop_window),
                                       ct.CustomToTensor(),
                                       ct.CustomNormalize()])
    
    train_dataset = mmDataSetTrain(phase_ims = phase_ims, 
                                mask_ims  = mask_ims,
                                transform = composed_train)


    valid_dataset = mmDataSetTrain(phase_ims = phase_ims, 
                                mask_ims  = mask_ims,
                                transform = composed_val)

    #check if GPU is available device
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    #load the unet architecture
    # net = NestedUNet()
    net = ShallowNestedUNet(deep_supervision=True) #deep supervision makes it longer to train but may have better results

    # saved_model='/hdd/RecPAIR/trained_unets/ShNestedUnet_mixed_brightnessAdj_Adam_BCE_256px_channels.pth'
    # saved_net = torch.load(saved_model)
    # net.load_state_dict(saved_net['model_state_dict'])

    #send the network to the GPU
    net.cuda()

    #initialise PyToch dataloader. Change batch size dependent on the GPU memory and crop window. 
    #num workers pre-loads data but reserves memory
    train_loader = DataLoader(train_dataset, batch_size=2, shuffle=True, num_workers = 0)

    #validation is optional.
    validation_loader = DataLoader(valid_dataset, batch_size=1, shuffle=True)

    #set up training parameters
    # adaptive learning rate optimization algorithm for faster learning
    # schedular reduces learning rate every step_size by gamma
    optimizer = torch.optim.Adam(net.parameters(),  lr=0.0001)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size = 10, gamma = 0.5)

    #how long will the network be trained
    num_epochs = 30

    #trainier will save the network to a 'save_name'. May overwrite old network.
    classifier = mmClassifier(net, optimizer, scheduler, 0.9, num_epochs, save_name)
    classifier.train(train_loader, validation_loader)

start = time.time()

# train network
TrainWholeNet(save_name='D:\\Python\\Kuba Unet\\test_train\\test_512px_agarosePad.pth')

end = time.time()
print(f"Segmentation took {round(end-start,1)} s.")