
function [output] = get_wavelet_stat(lineage_data, MCL_data)

% Get wavelet statistics
%m = 3;

output = {};

for m = 1:length(MCL_data)
    
temp = MCL_data{m}.Rsc_data.wavelet;
freq = temp.freqHr;
ROI = temp.bootstrap.ROI;

ccT = MCL_data{m}.ccT;  % in minute;
ccT = ccT/60;  % in hour

param = {};
param.short = [0.1 0.75];
param.cc = [0.75 1.5];
param.long = [1.5 3];

% Setting up short, cell cycle, and long-term time scale
% using cell cycle time
time_scale = {};
time_scale.short = ccT* param.short;
time_scale.cc = ccT* param.cc;
time_scale.long = ccT* param.long;

% Convert time scales into frequency (1/hour)
frequency_scale = {};
frequency_scale.short = 1./time_scale.short;
frequency_scale.cc = 1./time_scale.cc;
frequency_scale.long = 1./time_scale.long;

% Collect significant statistics on wavelet power spectrum

WL_sig_region = zeros(3, size(ROI,2));

for t = 1:size(ROI,2)
    
    for fq = 1:size(ROI,1)
        
        freq_temp = freq(fq);
        
        if (freq_temp < frequency_scale.short(1)) &&  (freq_temp >= frequency_scale.short(2))            

            WL_sig_region(1,t) =  WL_sig_region(1,t) + ROI(fq,t);
            
        end
        
        if (freq_temp < frequency_scale.cc(1)) &&  (freq_temp >= frequency_scale.cc(2))            

            WL_sig_region(2,t) =  WL_sig_region(2,t) + ROI(fq,t);
            
        end
        
        if (freq_temp < frequency_scale.long(1)) &&  (freq_temp >= frequency_scale.long(2))            

            WL_sig_region(3,t) =  WL_sig_region(3,t) + ROI(fq,t);
            
        end
        
    end
    
end

% Calculate the fraction of time that each frequency scale having 
minimal_freq_num = 3;
WL_sig_regionF = (WL_sig_region > minimal_freq_num);
WL_sig_fraction = sum(WL_sig_regionF,2)/size(ROI,2);

% Write into output 
output{m}.ccT = ccT;
output{m}.param = param;
output{m}.time_scale = time_scale;
output{m}.frequency_scale = frequency_scale;
output{m}.WL_sig_region = WL_sig_region;
output{m}.WL_sig_regionF = WL_sig_regionF;
output{m}.WL_sig_fraction = WL_sig_fraction;
output{m}.minimal_freq_num = minimal_freq_num;

end

%{
figure;
subplot(211); imagesc(ROI);
subplot(212); imagesc(WL_sig_regionF);

figure; plot(WL_sig_regionF');
%}
