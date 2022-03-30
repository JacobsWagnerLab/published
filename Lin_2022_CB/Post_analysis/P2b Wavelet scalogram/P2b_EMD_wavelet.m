
% Plot scalogram for ATP dynamics in mother cell lineage

% (1) Load MCL-PS and MCL-summary data
path = {};
path.load = '/Volumes/Data_04/Wei-Hsiang Lin/uF_dataset/exponential/MCL/';
dataset_nametag = 'M9GlcCA_lac';

cd(path.load);   
dataset_name0 = strcat('MCL_combined_', dataset_nametag, '.mat');
dataset_name1 = strcat('MCL_PS_', dataset_nametag, '.mat');
data0 = load(dataset_name0).lineage_data;
data1 = load(dataset_name1).MCL_data;

% ==============================================================

% (2) Plot scalogram

m = 4;  % Choose the MCL lineage (m is index number) for demonstration

tag1 = 'Rsc_data';
demo0 = data1{m}.MC;
demo1 = data1{m}.(tag1);
x = demo1.time;
y = demo1.x0;
ybase = demo1.baseline;
yDT = demo1.xDT;

% Get the frequency axis from scalogram
wavelet = abs(demo1.wavelet.coef);
wavelet_freq = demo1.wavelet.freqHr;

% Calculate the average cell cycle length, and the 75% to 150% range
ccT = (data1{m}.ccT/60);
ccT_freq = 1/ccT;
ccT_freq_range = [1/(1.5*ccT) 1/(0.75*ccT)]; 

% Generate plot

figure('position', [1 1 900 350]);

ccT_freq_ind = interp1(wavelet_freq, 1:size(wavelet), ccT_freq);
ccT_freq_range_ind = interp1(wavelet_freq, 1:size(wavelet), ccT_freq_range);

qx = x';
qc = ones(length(qx),1);
qymean = qc* ccT_freq_ind;
qyR1 = qc* ccT_freq_range_ind(1);
qyR2 = qc* ccT_freq_range_ind(2);
vx = [ 1:length(demo0.flag) ]';
vy = 20*(demo0.flag - 0.5);

% (i) Plot the detrended ATP trajectory

subplot(3,1,1);
plot(x, yDT, 'k-');     hold on;
plot(vx, vy, 'r-');  hold off;

ylim([-1.5 1.5]);
 
% (ii) Scalogram

subplot(3,1,2:3);

imagesc(wavelet, [0,1]); hold on;
yticks = 1:5:size(wavelet,1);
yticklabels(wavelet_freq(yticks));
set(gca,'xtick',[]);

RYBmap = customcolormap(linspace(0,1,11), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd','#def4f9','#abd9e9','#73add2','#4873b5','#313691'}, 128);
colormap(RYBmap);

% (iii) Overlay the cell division frequecy 150% and 75%

plot(qx, qyR1, 'w-', 'linewidth', 2);  hold on;
plot(qx, qyR2, 'w-', 'linewidth', 2);  hold off;

