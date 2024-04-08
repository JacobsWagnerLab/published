
function [output] = Img_GS_conv(input, param)

% Creating 2D normal distribution (symmetric) for image convolution
%width = 3;
%sigma = 1;

width = param(1);
sigma = param(2);

f = PM_2D_gaussian(width, sigma);

output= conv2(input, f, 'same'); 

%imagesc(f2)



