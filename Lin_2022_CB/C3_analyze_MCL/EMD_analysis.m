
function [output] = EMD_analysis(x, frame_interval)
% x: time series data
% frame_interval: sample interval of the data 

% =====================================================================
% Perform EMD analysis using Matlab built-in function
[IMF, res] = emd(x);

% Calculate the mean wavelength and STD for each IMF
num  = size(IMF,2);
wave_length = NaN(num,2);

for j = 1:num
    
    % Peaks
    [valP,indP] = findpeaks(IMF(:,j));
    tempP = indP(2:end) - indP(1:end-1);    
    wave_length(j,1) = frame_interval *mean(tempP);
    wave_length(j,3) = frame_interval *std(tempP);
    
    % Valleys
    [valV,indV] = findpeaks(-IMF(:,j));
    tempV = indV(2:end) - indV(1:end-1);    
    wave_length(j,2) = frame_interval *mean(tempV);
    wave_length(j,4) = frame_interval *std(tempV);
    
    % Mean wave length
    wave_length(j,5) = nanmean([wave_length(j,1:2)]);    
    wave_length(j,6) = nanmean([wave_length(j,3:4)]);    

    % Contributed STD on magnitude
    wave_length(j,7) = std(IMF(:,j));

end
% =====================================================================

output = {};
output.IMF = IMF;
output.res = res;
output.wave_length = wave_length;


