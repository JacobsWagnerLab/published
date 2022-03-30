
% Power spectrum analysis with shuffle bootstrap statistics

function [PS] = PS_analysis(x, frame_interval)

% (1) Analyze for time series x
PS = {};
PS.x = x;
[PS.freq, PS.PS] = power_spectrum(PS.x, frame_interval);

% (2) Shuffle the original data to obtain the "null hypothesis" 
shuffle_num = 100;
rec = NaN(size(PS.PS,1), shuffle_num);

parfor j = 1:shuffle_num   
    xs = shuffle_data(x);
    [freq, rec(:,j)] = power_spectrum(xs, frame_interval);    
end

PS.SH_mean = mean(rec,2);
PS.SH_ste = std(rec,0,2)/sqrt(shuffle_num);

end

%{
figure; 
subplot(121);
plot(PS.freq, PS.PS, 'o-'); 

subplot(122);
errorbar(PS.freq, PS.SH_mean, PS.SH_ste, 'o-'); 
%}
