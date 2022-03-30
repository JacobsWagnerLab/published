
% wavelet PS summary

function [PS, PSSH] = get_waveletPS_summary(MCL_data, tag)
%tag = 'Rsc_data';

num = length(MCL_data);

p = 60;
freq_mat = NaN(p, num);
PS_mat = NaN(p, num);
PS_mat_SH = NaN(p, num);

% (1) Examine the frequency scale (should be all identical)
for m = 1:num
    
    Rsc_WLfreq = MCL_data{m}.(tag).wavelet.freqHr;    
    M = min(p, length(Rsc_WLfreq));
    freq_mat(1:M,m) = Rsc_WLfreq(1:M);    
    
end

% (2) Collect the Power spectrum from wavelet 2D map
for m = 1:num
    
    Rsc_WLPS = MCL_data{m}.(tag).wavelet.PS;    
    M = min(p, length(Rsc_WLPS));
    PS_mat(1:M, m) = Rsc_WLPS(1:M);    
    
    Rsc_WLPSSH = MCL_data{m}.(tag).wavelet.bootstrap.PSmean;  
    PS_mat_SH(1:M, m) = Rsc_WLPSSH(1:M);
    
end

PS = {};
PS.mat = PS_mat;
PS.mean = nanmean(PS.mat, 2);
PS.se = nanstd(PS.mat,0,2)/sqrt(num);
PS.freq = freq_mat(:,1);  % Note that all wavelet matrix has same frequency scale

PSSH = {};
PSSH.mat = PS_mat_SH;
PSSH.mean = nanmean(PSSH.mat, 2);
PSSH.se = nanstd(PSSH.mat,0,2)/sqrt(num);
PSSH.freq = PS.freq;

%{
figure; 
subplot(211); imagesc(PS.mat);
subplot(212); imagesc(PSSH.mat); 
colormap(jet);

figure; 
errorbar(PS.freq, PS.mean, PS.se); hold on;
errorbar(PSSH.freq, PSSH.mean, PSSH.se); hold on;
%}

end

