%Written by Skye Glenn in the Jacobs-Wagner lab at Stanford University,
%published in 2023.
%
%ABOUT
%This code returns the slope of linear regressions fit to growth data 
%taken from a microplate reader (optical density, OD) over a sliding window
%to find the maximum growth rate. 
%
%INPUTS
%OD: A matrix in which each column is the OD over time of a single well on
%   a microplate reader. 
%frame_rate: The time interval between OD acquisitions on the microplate 
%   reader (min)
%wind: The window over which to smooth cell size (window)
%
%OUTPUTS
%slopes: Matrix (of the same size as OD) of slopes of a linear regression 
%   fit to log10 transformed data in OD
%max_gr: the maximum value for each column of slopes
%
%DEPENDENCIES
%none
%

%% Define the inputs
% Name matrix of ODs taken from microplate reader 'OD'
frame_rate = 10;           % time interval of OD acquisition on microplate reader (min)
temp_time = (frame_rate:frame_rate:length(OD)*frame_rate)';
wind = 25;                %size of the sliding window for smoothing (# frames)

%% Calculate slope of linear regression to semilog growth curve data
slopes = NaN(size(OD));
for c = 1:size(OD,2)
    for r = 1:size(OD,1)-wind
        temp_OD = log10(OD(:,c));
        temp_p = polyfit(temp_time(r:r+wind),temp_OD(r:r+wind),1);
        slopes(r,c) = temp_p(1);
    end
end
clear temp_OD;clear temp_p;

[max_gr,I_gr] = max(slopes);
max_gr_OD = OD(I_gr);