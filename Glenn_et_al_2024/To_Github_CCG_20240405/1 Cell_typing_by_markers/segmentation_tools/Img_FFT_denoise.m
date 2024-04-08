
function [output] = Img_FFT_denoise(input, wave_cut)

% Feb 1, 2019

% Reomve high-frequency noise.
% input: Image 
% wave_cut: cutoff-wavelength (in terms of pixel)   

% Generate Fourier transformed matrix in (s1,s2) frequency domain
FFT_mat = fft2(input);

% Filter out the high frequency (Only take component with longer wavelength)
% by setting some of the entry of FFT_filtered into 0.

s1_max = size(FFT_mat,1);
s2_max = size(FFT_mat,2);
s1_cut = floor(s1_max/wave_cut);
s2_cut = floor(s2_max/wave_cut);

FFT_filtered = FFT_mat;
FFT_filtered( s1_cut + 1 : s1_max - s1_cut ,: ) = zeros(s1_max - (2*s1_cut) , s2_max);
FFT_filtered( :, s2_cut+1 : s2_max - s2_cut ) = zeros(s1_max, s2_max-(2*s2_cut));
output = abs(ifft2(FFT_filtered));

