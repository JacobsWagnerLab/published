%% Order of image analysis for smFISH
%% 1. Obtain cell meshes by running MicrobeTracker

%% 2. Identify spots and obtain their info (e.g., magnitude, from gaussian fitting) by using spotFinder in MicrobeTracker
% spots in channel 1 (Cy5, Z5) and spots in channel 2 (Cy3, Z3) are saved
% in different cellList and hence mat files.

%% 3. Obtain threshold for spot curation
SegHist(cellList_timezero); % At time zero, segment intensities indicate background level. Get those and use as a threshold.
% From the histogram plot, I take the end of the Gaussian curve as the threshold.

%% 4. Find SM normalization factor
ListSpotmag(folderPath,channel,timePoints,thresh1,thresh2)
% Save the result (Xans, nx1 array) in Excel. Have results from different
% days as independent columns. Top of the column goes the description of the column e.g.,
% "day1z5spot" "day1z3spot". 
% Run SJ_Gaussian_Fitting_Clean in python to run gmm (gaussian mixture model)
% system('python SJ_Gaussian_Fitting_Clean.py');

%% 5. smFISH based on spot intensity only
masterfolder = 'J:\Shares\Data_02\Sangjin Kim\Microscope Data\FISH trials\12early_lacZ\';
dayList(1).folderpath = strcat(masterfolder,'140505-lacZ 2CFISH IPTG glucose pulse burstcheck-NSTORM\B.MGhighIPTGearly\');
dayList(2).folderpath = strcat(masterfolder,'140611-lacZ 2CFISH IPTG glucose pulse-NSTORM\B.MGhighearly\');
dayList(3).folderpath = strcat(masterfolder,'140613-lacZ 2CFISH IPTG glucose pulse high-NSTORM\B.MGhighearly\');
dayList(4).folderpath = strcat(masterfolder,'140716-lacZ 2CFISH IPTG glucose pulse high-NSTORM\A.MGhighearly\');

% SM normalization factor from above process (step #4)
daily_scaleZ5(1) = 0.36494524; daily_scaleZ3(1) = 0.49707307;
daily_scaleZ5(2) = 0.28242603; daily_scaleZ3(2) = 0.55073717;
daily_scaleZ5(3) = 0.29541658; daily_scaleZ3(3) = 0.49020809;
daily_scaleZ5(4) = 0.34345343; daily_scaleZ3(4) = 0.67065377;
daily_scale = [daily_scaleZ5', daily_scaleZ3'];

% threshold for spot curation (step #3)
daily_thresh1(1) = 0.003; daily_thresh2(1) = 0.005; 
daily_thresh1(2) = 0.0025; daily_thresh2(2) = 0.007;
daily_thresh1(3) = 0.003; daily_thresh2(3) = 0.006;
daily_thresh1(4) = 0.003; daily_thresh2(4) = 0.006;
daily_thresh = [daily_thresh1', daily_thresh2'];

% time course information
% time (minutes) of FISH 
timecourseInfo.time = [0,1,2,3,4,5,6,7,8,9,10,12]; 
% days for each time point
timecourseInfo.day = {[1,2,3,4],[1,2,3,4],[1,2,3,4],[2,3,4],[2,3,4],[2,3,4],[2,3,4],[2,3,4],[2,3,4],[3,4],[2,3],[2,3,4]};
% file name corresponding to the timepoint in each day
timecourseInfo.file = {[1,1,1,1],[2,2,2,2],[3,3,3,3],[4,4,4],[5,5,5],[6,6,6],[7,7,7],[8,8,8],[9,9,9],[10,10],[11,11],[12,12,12]}; 

Xans = RNApercell_by_spotmag_bootstrap(dayList,daily_scale,daily_thresh,timecourseInfo);

%% 5 (opt). smFISH based on total fluorescence signal
% Get avgbasePix from gmm fitting of 
load(cellListfile);
x(:,1) = meaninthist(cellList,'signal1');
x(:,2) = meaninthist(cellList,'signal2');
masterfolder = 'J:\Shares\Data_02\Sangjin Kim\Microscope Data\FISH trials\12early_lacZ\';
dayList(1).folderpath = strcat(masterfolder,'140505-lacZ 2CFISH IPTG glucose pulse burstcheck-NSTORM\B.MGhighIPTGearly\');
dayList(2).folderpath = strcat(masterfolder,'140611-lacZ 2CFISH IPTG glucose pulse-NSTORM\B.MGhighearly\');
dayList(3).folderpath = strcat(masterfolder,'140613-lacZ 2CFISH IPTG glucose pulse high-NSTORM\B.MGhighearly\');
dayList(4).folderpath = strcat(masterfolder,'140716-lacZ 2CFISH IPTG glucose pulse high-NSTORM\A.MGhighearly\');

% SM normalization factor from above process (step #4)
daily_scaleZ5(1) = 0.36494524; daily_scaleZ3(1) = 0.49707307;
daily_scaleZ5(2) = 0.28242603; daily_scaleZ3(2) = 0.55073717;
daily_scaleZ5(3) = 0.29541658; daily_scaleZ3(3) = 0.49020809;
daily_scaleZ5(4) = 0.34345343; daily_scaleZ3(4) = 0.67065377;
daily_scale = [daily_scaleZ5', daily_scaleZ3'];

% average pixel value at time zero.
% get this by getting mean pixcel values zt time zero and use gmm code (in
% step #4) with k = 1, one gaussian fit
% x(:,1) = meaninthist(cellList,'signal1'); x(:,2) = meaninthist(cellList,'signal2');
daily_avgbasePix1(1) = 0.0005091; daily_avgbasePix2(1) = 0.00071035;
daily_avgbasePix1(2) = 0.0005255; daily_avgbasePix2(2) = 0.00079067;
daily_avgbasePix1(3) = 0.00061279; daily_avgbasePix2(3) = 0.00116339;
daily_avgbasePix1(4) = 0.00049022; daily_avgbasePix2(4) = 0.00062315;
daily_avgbasePix = [daily_avgbasePix1', daily_avgbasePix2'];

% time course information
% time (minutes) of FISH 
timecourseInfo.time = [0,1,2,3,4,5,6,7,8,9,10,12];
% days for each time point
timecourseInfo.day = {[1,2,3,4],[1,2,3,4],[1,2,3,4],[1,2,3,4],[2,3,4],[2,3,4],[2,3,4],[2,3,4],[2,3,4],[3,4],[2,3],[2,3,4]};
% file name corresponding to the timepoint in each day
timecourseInfo.file = {[1,1,1,1],[2,2,2,2],[3,3,3,3],[4,4,4,4],[5,5,5],[6,6,6],[7,7,7],[8,8,8],[9,9,9],[10,10],[11,11],[12,12,12]};

Xans = RNApercell_by_totalfluorescence_bootstrap(dayList,daily_scale,daily_avgbasePix,timecourseInfo);

%% 6. To estimate premature termination
% example data: Promoter turn-off at t = 90 s after 0.2 mM IPTG induction
% (one day example)
% mRNA lifetime (mean) was 1.52 min (Z5) and 1.66 min (Z3)
X = [0,0.0134487000000000,0.0181838000000000;1,1.79951000000000,0.0261977000000000;2,4.73399000000000,0.122877000000000;3,3.79987000000000,0.287993000000000;4,3.08398000000000,1.61866000000000;5,1.32696000000000,1.80267000000000;6,0.976509000000000,1.54928000000000;7,0.587405000000000,0.963588000000000;8,0.288642000000000,0.430795000000000;9,0.206146000000000,0.236055000000000;10,0.172012000000000,0.194890000000000;12,0.104088000000000,0.0946338000000000];
intP = IntegrationFISH(X,[1.52,1.66]);