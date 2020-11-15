function [yHat] = center_interp_signal(signal, center)
%a signal is shifted to its center and uniformly interpolated. This is a
%helper function not meant to be used elsewhere.
%
% signal is any 1D vector
% center is an index location in signal
%
% yhat is a shifted and interpolated version of signal that has been
% sampled 300 times. locations missing after the shift operation are
% represented as NaN
%
% Brad Parry, Christine Jacobs-Wagner lab; 2016 April


signal = signal(:)';
x = (1:length(signal)) / length(signal) * 100;
%the interpolant center
peakHat = center / length(signal) * 100;

xhat = (1/3):(1/3):100;
y_interp = pchip(x, signal, xhat);
[~,ix] = min(abs(peakHat-xhat));

shift = round(length(xhat)/2) - ix;
xhatIndex = (1:length(xhat)) + shift;
xhatIndex(xhatIndex < 1) = [];
xhatIndex(xhatIndex > length(xhat)) = [];

yHat = NaN(1,length(xhat));
yHat(xhatIndex) = y_interp(xhatIndex-shift);