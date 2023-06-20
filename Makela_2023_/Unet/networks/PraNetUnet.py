
# Project: pipeline v0.1# Author: Praneeth Karempudi
# Email: praneeth.karempudi@gmail.com

##############################################################################
# This is the basic U net architecture, almost straight from the paper
# titled "U-Net: Convolutional Networks for Biomedical Image segmentation"
# by Olaf Ronneberger, Philipp Fischer and Thomas Brox, 2015
# PyTorch impelementation

#stitched by kuba from 2 functions written by Praneeth

import torch
import torch.nn as nn
import torch.nn.functional as F

class PraNetUnet(nn.Module):
    '''
    Basic Unet architecture, full size.
    '''

    def __init__(self):
        # initialize the superclass initializer from the library modules
        super(PraNetUnet, self).__init__()
        # 1 because of the number of input channels is 1
        self.initial = double_conv(1, 64)
        self.down1 = double_conv(64, 128)
        self.down2 = double_conv(128, 256)
        self.down3 = double_conv(256, 512)
        self.down4 = double_conv(512, 512)

        self.up1 = up_conv_cat(1024, 256)
        self.up2 = up_conv_cat(512, 128)
        self.up3 = up_conv_cat(256, 64)
        self.up4 = up_conv_cat(128, 64)
        self.out = nn.Conv2d(64, 1, 1) # 2 because of no_classes

    def forward(self, x):
        # x will be the image batch tensor, that will be propagated across
        # toward the end and the 
        x1 = self.initial(x)
        x2 = self.down1(F.max_pool2d(x1, (2, 2)))
        x3 = self.down2(F.max_pool2d(x2, (2, 2)))
        x4 = self.down3(F.max_pool2d(x3, (2, 2)))
        x5 = self.down4(F.max_pool2d(x4, (2, 2)))

        # copied code
        x = self.up1(x5, x4)
        x = self.up2(x, x3)
        x = self.up3(x, x2)
        x = self.up4(x, x1)
        x = self.out(x)

        return x

class ShallowPraNetUnet(nn.Module):
    '''
    Basic Unet architecture, full size.
    '''

    def __init__(self):
        # initialize the superclass initializer from the library modules
        super(ShallowPraNetUnet, self).__init__()
        # 1 because of the number of input channels is 1
        self.initial = double_conv(1, 64)
        self.down1 = double_conv(64, 128)
        self.down2 = double_conv(128, 256)
        self.down3 = double_conv(256, 256)
        #self.down4 = double_conv(512, 512)

        #self.up1 = up_conv_cat(1024, 256)
        self.up1 = up_conv_cat(512, 128)
        self.up2 = up_conv_cat(256, 64)
        self.up3 = up_conv_cat(128, 64)
        self.out = nn.Conv2d(64, 1, 1) # 2 because of no_classes

    def forward(self, x):
        # x will be the image batch tensor, that will be propagated across
        # toward the end and the 
        x1 = self.initial(x)
        x2 = self.down1(F.max_pool2d(x1, (2, 2)))
        x3 = self.down2(F.max_pool2d(x2, (2, 2)))
        x4 = self.down3(F.max_pool2d(x3, (2, 2)))
        #x5 = self.down4(F.max_pool2d(x4, (2, 2)))

        # copied code
        x = self.up1(x4, x3)
        x = self.up2(x, x2)
        x = self.up3(x, x1)
        #x = self.up4(x, x1)
        x = self.out(x)

        return x
        
#leyers for UNet

class double_conv(nn.Module):
    '''
    Combining conv1- batch - relu - conv batch relu parts of UNet
    
    Features: No change in size of the image (maintained using padding)
    '''
    
    def __init__(self, input_channels, output_channels):
        ''' Just initializes the conv - batch - relu layers twice 
            with no reduction in image size (don't change padding , stride)
        '''
        super(double_conv, self).__init__()
        self.conv = nn.Sequential(
            nn.Conv2d(input_channels, output_channels, 3, stride = 1, padding = 1),
            nn.BatchNorm2d(output_channels),
            nn.ReLU(inplace = True),
            nn.Conv2d(output_channels, output_channels, 3, stride = 1, padding =1),
            nn.ReLU(inplace = True)
        )

    def forward(self, x):
        x = self.conv(x)
        return x



class up_conv_cat(nn.Module):
    '''
    UpConv + concatenation of features during downscaling
    '''

    def __init__(self, input_channels, output_channels):
        super(up_conv_cat, self).__init__()

        #self.up = nn.Upsample(scale_factor=2, mode='bilinear', align_corners=True)
        self.up = nn.ConvTranspose2d(input_channels//2, input_channels//2, 2, stride=2)

        self.conv = double_conv(input_channels, output_channels)

    def forward(self, from_bottom, from_side):

        x1 = self.up(from_bottom)

        x2 = torch.cat([from_side, x1], dim = 1)

        x2 = self.conv(x2)
 
        return x2
