
% Use matlab function to perform wavelet analysis, 

function [wavelet] = wavelet_analysis(x, frequency_Hz)

wavelet = {};

% Wavelet analysis for Rsc data

[wavelet.coef, wavelet.freqHz] = cwt(x, 'amor', frequency_Hz);
wavelet.freqHr = wavelet.freqHz*3600; 
wavelet.bootstrap = wavelet_bootstrap(x, wavelet.coef, frequency_Hz);
wavelet.PS = mean(abs(wavelet.coef), 2);
wavelet.time = (1:length(x))/(frequency_Hz*3600);

%{
figure;
subplot(211);  imagesc( abs(wavelet.coef) );
subplot(212);  imagesc( wavelet.bootstrap.ROI );
colormap(jet);
%}

end