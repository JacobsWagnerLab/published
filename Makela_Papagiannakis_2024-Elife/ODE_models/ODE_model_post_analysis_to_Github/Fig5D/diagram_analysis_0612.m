
% Plot the 2D diagram of RNAP and ribosome active fraction
% Used for Figure 5E, and Figure 5E-supplementary figure 1

clear all;

common_dir = '/Users/wei-hsiang/Desktop/DNA content model/';
data_directory = strcat(common_dir, 'A1_exp_data/raw_data/');
simu_directory_ModelA = strcat(common_dir, 'A2_simu_data/0329/');
simu_directory_ModelB = strcat(common_dir, 'A2_simu_data/0523A/');

RNAP_name = {};
ribo_name = {};
simu_name = {};

RNAP_name.W1 = 'RNAP_WT_fracs_raw';
ribo_name.W1 = 'Ribosome_WT_fracs_raw';

RNAP_name.Z1 = 'RNAP_ftsZ_fracs_raw';
ribo_name.Z1 = 'Ribosome_ftsZ_fracs_raw';

RNAP_name.C1 = 'RNAP_oriC_fracs_raw';
ribo_name.C1 = 'Ribosome_oriC_fracs_raw';

RNAP_name.W2 = 'M9gly_RNAP_WT_fracs_raw';
ribo_name.W2 = 'M9gly_Ribosome_WT_fracs_raw';

RNAP_name.C2 = 'M9gly_RNAP_oriC_fracs_raw';
ribo_name.C2 = 'M9gly_Ribosome_oriC_fracs_raw';

RNAP_name.W3 = 'M9ala_RNAP_WT_fracs_raw';
ribo_name.W3 = 'M9ala_Ribosome_WT_fracs_raw';

RNAP_name.C3 = 'M9ala_RNAP_oriC_fracs_raw';
ribo_name.C3 = 'M9ala_Ribosome_oriC_fracs_raw';

simu_name.S1 = 'data_summaryC1';   % for ModelA, dataset1
simu_name.S2 = 'data_summaryC2';   % for ModelA, dataset2
simu_name.S3 = 'data_summaryC3';   % for ModelA, dataset3
simu_name.SPQ = 'data_summaryC';   % for ModelB, dataset1


% (1) Read experimental data and perform binning

bin_data_rec = {};
exp_repeat_num = 3;  %    
    
[WT_data1] = bin_data(RNAP_name.W1, ribo_name.W1, exp_repeat_num);
[ftsZ_data1] = bin_data(RNAP_name.Z1, ribo_name.Z1, exp_repeat_num);    
[oriC_data1] = bin_data(RNAP_name.C1, ribo_name.C1, exp_repeat_num);
    
bin_data_rec{1}.wt = WT_data1;    
bin_data_rec{1}.ftsZ = ftsZ_data1;    
bin_data_rec{1}.oriC = oriC_data1;

[WT_data2] = bin_data(RNAP_name.W2, ribo_name.W2, exp_repeat_num);  
[oriC_data2] = bin_data(RNAP_name.C2, ribo_name.C2, exp_repeat_num);
    
bin_data_rec{2}.wt = WT_data2;        
bin_data_rec{2}.oriC = oriC_data2;

[WT_data3] = bin_data(RNAP_name.W3, ribo_name.W3, exp_repeat_num);  
[oriC_data3] = bin_data(RNAP_name.C3, ribo_name.C3, exp_repeat_num);
    
bin_data_rec{3}.wt = WT_data3;        
bin_data_rec{3}.oriC = oriC_data3;

% (2) Read simulation data

cd(simu_directory_ModelA);

temp = {};
diagram = {};

temp{1} = load(simu_name.S1).data_summaryC.simudata{1};
temp{2} = load(simu_name.S2).data_summaryC.simudata{1};
temp{3} = load(simu_name.S3).data_summaryC.simudata{1};

for DB = 1:3
    
    diagram{DB}.F.AF_RNAP = temp{DB}.AF.F.bin_AF_RNAP(:,2);
    diagram{DB}.F.AF_ribo = temp{DB}.AF.F.bin_AF_ribo(:,2);

    diagram{DB}.G.AF_RNAP = temp{DB}.AF.G.bin_AF_RNAP(:,2);
    diagram{DB}.G.AF_ribo = temp{DB}.AF.G.bin_AF_ribo(:,2);

end


tempPQ = {};
diagramPQ = {};

cd(simu_directory_ModelB);
tempPQ{1} = load(simu_name.SPQ).data_summaryC.simudata{1};


for DB = 1:1
    
    diagramPQ{DB}.F.AF_RNAP = tempPQ{DB}.AF.F.bin_AF_RNAP(:,2);
    diagramPQ{DB}.F.AF_ribo = tempPQ{DB}.AF.F.bin_AF_ribo(:,2);
    
    diagramPQ{DB}.G.AF_RNAP = tempPQ{DB}.AF.G.bin_AF_RNAP(:,2);
    diagramPQ{DB}.G.AF_ribo = tempPQ{DB}.AF.G.bin_AF_ribo(:,2);

    
end


% =====================================================================
% Scatter plot on diagram
% =====================================================================

plot_flag = 1;

color = {};
color.G0 = [50 50 50]/255;
color.B1 = [71 88 248]/255;
color.B2 = [68 52 206]/255;
color.Y1 = [251 210 47]/255;
color.Y2 = [245 172 80]/255;
color.G1 = [166 168 171]/255;
color.G2 = [48 48 48]/255;

% Calculate DNA concentration for multi-N and 1N cells

DNA_conc_multiN = 1.4;       % copy/um^3
area_to_volume = 0.6;

area_list = temp{1}.AF.area_axis;  % um^3 
volume_list = area_list* area_to_volume;
DNA_conc_1N = 1./volume_list;
ind_gap = 100;

ctemp = parula;
parulaB = ctemp(1:end-30,:);
color_mat = flip(parulaB);

close all;


figure;
imagesc(1:256-30);
colormap(color_mat);


figure('position', [1 1 300 350]);

% (1)  oriC data

RNAP = oriC_data1.RNAP;
ribo = oriC_data1.ribo;
bin_num = 12;
sigma = 2.5;  % 95% condidence interval

for b = 1:bin_num
    
    DNA_conc = DNA_conc_1N(2*b);  
    % Factor of 2 is for adjusting the scale. 
    % Each bin in simulation data is 0.5 um^3,
    % each bin in experiment data is 1.0 um^3
    
    DNA_conc_ind = floor(DNA_conc *ind_gap);
    
    errorbar(RNAP.stat(b,2), ribo.stat(b,2), ...
             sigma *ribo.stat(b,4), sigma *ribo.stat(b,4), ...
             sigma *RNAP.stat(b,4), sigma *RNAP.stat(b,4), ...
             '', 'color', color_mat(DNA_conc_ind,:));
    
    hold on;
     
end
%color(color_mat);

% (2) ftsZ data

RNAP = ftsZ_data1.RNAP;
ribo = ftsZ_data1.ribo;
bin_num = 12;
sigma = 2.5;  % 95% condidence interval

for b = 1:bin_num
    
    DNA_conc_ind0 = DNA_conc_multiN * ind_gap;
    
    errorbar(RNAP.stat(b,2), ribo.stat(b,2), ...
             sigma *ribo.stat(b,4), sigma *ribo.stat(b,4), ...
             sigma *RNAP.stat(b,4), sigma *RNAP.stat(b,4), ...
             '', 'color', color_mat(DNA_conc_ind0,:));    
    hold on;
     
end


% (3) Plot simulation data of 1N cells

for b = 1:bin_num
    
    ind = 2*b;
    DNA_conc = DNA_conc_1N(ind);
    
    DNA_conc_ind = floor(DNA_conc *ind_gap);   % convert into the same scale as the experimental data
    
    % Model A
    plot(diagram{1}.G.AF_RNAP(ind), diagram{1}.G.AF_ribo(ind), '.', 'color', color_mat(DNA_conc_ind,:), 'MarkerSize', 20 );  hold on;    
    
    % Model B
    %plot(diagramPQ{1}.G.AF_RNAP(ind), diagramPQ{1}.G.AF_ribo(ind), '.', 'color', color_mat(DNA_conc_ind,:), 'MarkerSize', 20 ); hold on;
    
end

% (4) Plot simulation data of multi-N cells

DNA_conc_ind0 = DNA_conc_multiN * ind_gap;


% Model A
plot(diagram{1}.F.AF_RNAP(end), diagram{1}.F.AF_ribo(end), '.','color',  color_mat(DNA_conc_ind0,:), 'MarkerSize', 20); hold on;

% Model B 
%plot(diagramPQ{1}.F.AF_RNAP(end), diagramPQ{1}.F.AF_ribo(end), '.','color',  color_mat(DNA_conc_ind0,:), 'MarkerSize', 20); hold on;

xlim([0 0.6]);
ylim([0 1.0]);

%
% ======================================================================

function [output] = bin_data(RNAP_name, ribo_name, n)

%(1) Load RNAP data

RNAP = {};
RNAP.data = load(RNAP_name);
RNAP.table = {};

for j = 1:n
    
    temp1 = RNAP.data.areas{j};
    temp2 = RNAP.data.fracs{j};
    RNAP.table{j} = [temp1' temp2'];
    
end

RNAP.tableC = [];

for j = 1:n

    RNAP.tableC = [RNAP.tableC; RNAP.table{j}];
    
end

% (2) Load ribo data

ribo = {};
ribo.data = load(ribo_name);
ribo.table = {};

for j = 1:n
    
    temp1 = ribo.data.areas{j};
    temp2 = ribo.data.fracs{j};
    ribo.table{j} = [temp1' temp2'];
    
end

ribo.tableC = [];

for j = 1:n

    ribo.tableC = [ribo.tableC; ribo.table{j}];
    
end


% ----------------------------------------------------------------------
% (0) Determine cell size bin
% ----------------------------------------------------------------------

bin_size = 1;  % cell area
bin_num = 12;  % number of bins

data_range = NaN(bin_num, 2); % min and max for each bin

for b = 1:bin_num
    
    data_range(b,1) = bin_size *(b-1);
    data_range(b,2) = bin_size * b;
    
end

% ----------------------------------------------------------------------
% (i) RNAP: combine data from all table and partition by cell area bin
% ----------------------------------------------------------------------

for b = 1:bin_num
    
    bmin = data_range(b,1);
    bmax = data_range(b,2);

    temp = RNAP.tableC;
    temp( (temp(:,1) < bmin), : ) = [] ;
    temp( (temp(:,1) >= bmax), : ) = [] ;
        
    RNAP.bindata{b} = temp;
    
end

RNAP.stat = NaN(bin_num,3);  % num, mean, SD, SE

for b = 1:bin_num
    
    RNAP.stat(b,1) = size(RNAP.bindata{b}, 1);
    RNAP.stat(b,2) = nanmean(RNAP.bindata{b}(:,2));
    RNAP.stat(b,3) = nanstd(RNAP.bindata{b}(:,2));    
    RNAP.stat(b,4) = RNAP.stat(b,3)./sqrt(RNAP.stat(b,1));
    
end


% ----------------------------------------------------------------------
% (ii) Ribosome: combine data from all table and partition by cell area bin
% ----------------------------------------------------------------------

for b = 1:bin_num
    
    bmin = data_range(b,1);
    bmax = data_range(b,2);

    temp = ribo.tableC;
    temp( (temp(:,1) < bmin), : ) = [] ;
    temp( (temp(:,1) >= bmax), : ) = [] ;
        
    ribo.bindata{b} = temp;
    
end

ribo.stat = NaN(bin_num,3);  % num, mean, SD, SE

for b = 1:bin_num
    
    ribo.stat(b,1) = size(ribo.bindata{b}, 1);
    ribo.stat(b,2) = nanmean(ribo.bindata{b}(:,2));
    ribo.stat(b,3) = nanstd(ribo.bindata{b}(:,2));    
    ribo.stat(b,4) = ribo.stat(b,3)./sqrt(ribo.stat(b,1));
    
end

output = {};
output.RNAP = RNAP;
output.ribo = ribo;

end

% ======================================================================


