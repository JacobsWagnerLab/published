
% Using bootstrap method to test statistical significant of wavelet analysis

function [null_dist] = wavelet_bootstrap(x0, c0, frequency_Hz)

%x0 = MCL_data{1}.PS_data.Rsc.mc;
%c0 = abs(MCL_data{1}.wavelet_data.Rsc{2}.coef);

shuffle_num = 100;

% Shuffle the time series to get random wavelet coefs.

rec = NaN(size(c0,1), size(c0,2), shuffle_num);
recPS = NaN( size(c0,1), shuffle_num );

parfor j = 1:shuffle_num
    
    temp = shuffle_data(x0);
    [coef, freq] = cwt(temp, 'amor', frequency_Hz);    
    
    rec(:,:,j) = abs(coef);
    recPS(:,j) = mean(abs(coef), 2);
    
end

% Calculate the null distribution

null_dist = {};
null_dist.c0 = abs(c0);
null_dist.mean = nanmean(rec, 3);
null_dist.std = nanstd(rec, 0, 3);
null_dist.PSmean = nanmean(recPS, 2);
null_dist.PSstd = nanstd(recPS, 0, 2);

sigma = 2;
null_dist.upperbound = null_dist.mean + sigma* null_dist.std;
null_dist.ROI = (null_dist.c0 > null_dist.upperbound);

end

%{
figure;
subplot(311); imagesc(c0);
subplot(312); imagesc(null_dist.upperbound);
subplot(313); imagesc(null_dist.ROI);
colormap(jet);
%}