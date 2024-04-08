
function [binary_output] = Mask_adap_thr(input_image, param)

%%% Oct18, 2018
% This is a segmentation algorithm using adaptive thresholding
% depended on local intensity. Useful to segment the fluorescence image  
% with objects having quite different intensities.

% Input: an image matrix 
% Output: a binary matrix indicating the thresholded image

% Parameters: 
% (1) local_radius: specify the radius (in pixel) of local region R used for thresholding
% (2) minimal_level:  specify the minimal thresholding level
% (3) std_level: the local thresholding criteria is (mean(R) + (std_level)*std(R)) 

input_image = double(input_image);
local_radius = param(1);
minimal_level = param(2);
std_level = param(3);

%%% 

% (1) Generate a 2D gaussian distribution on a sqaure for local sampling

L = 2*local_radius+1;
local_sampling_region = PM_2D_gaussian(L, L);
local_sampling_region = local_sampling_region / sum(sum( local_sampling_region ));

% (2) For each pixel, generating a local neighborhood around this pixel
% and determine this pixel pass thresholding or not. 

% Practiacally, this code is performed by convolution of input_image
% to local_sampling_region to increase performance. 

mean_matrix = conv2(input_image, local_sampling_region, 'same');
var_matrix = ( conv2((input_image.^2), local_sampling_region, 'same' ) ) - (mean_matrix.^2);
std_matrix = sqrt(var_matrix);

local_level = mean_matrix + (std_level * std_matrix) ;

% For every pixel, it will pass the thresolding only if following two criteria 
% are satisfied:              
% (i) it is larger than the local threshold level (local_level)
% (ii) it is larger than the global minimal level (minimal_level)

binary_output = ( (input_image > local_level ) .* (input_image > minimal_level ) );


