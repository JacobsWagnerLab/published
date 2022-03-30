
% Perform sibling statistics (correlation plot) and compare with shuffled
% statistics. Manaully load 'ATP_stat_type.mat'

path = {};
path.save = '/Users/wei-hsiang/Desktop/uF_local/Tree_data/sibling_stat_corr/';

dataset_name = {'M9GlcCA_lac', 'M9GlcCA_ldhA', 'M9GlcCA_adhE', 'M9GlcCA_pta', 'M9GlcCA_ackA', ...
                'M9Glc_lac', 'M9Xyl_lac', 'M9Gly_lac'};

% -----------------------------------------------------------------------            
DB_num = 8;
mtype = 1;
datatag1 = 'Rsc_stat';
datatag2 = 'mean';

for DB = 1:DB_num
    
dataA = sibling_stat{DB}.data.(datatag1).(datatag2){mtype};
dataASH = shuffle_data( dataA(:,1:2), 2 );

% Calculate correlation coefficient
corr_stat = {};
corr_stat.P = corr(dataA(:,1), dataA(:,2), 'type','Pearson');
corr_stat.K = corr(dataA(:,1), dataA(:,2), 'type','Kendall');
corr_stat.S = corr(dataA(:,1), dataA(:,2), 'type','Spearman');
corr_stat.P2 = corr(dataASH(:,1), dataASH(:,2), 'type','Pearson');
corr_stat.K2 = corr(dataASH(:,1), dataASH(:,2), 'type','Kendall');
corr_stat.S2 = corr(dataASH(:,1), dataASH(:,2), 'type','Spearman');

% -----------------------------------------------------------------------

rg = [2 8];

h1 = figure('position', [1 1 600 250]); 

subplot(121);
plot (dataA(:,1), dataA(:,2), '.'); 
title_string1 = strcat('corrP=', num2str(corr_stat.P), ' corrK=', num2str(corr_stat.K), ' corrS=', num2str(corr_stat.S) );
title(title_string1);

xlim(rg);
ylim(rg);

subplot(122);
plot (dataASH(:,1), dataASH(:,2), '.'); 
title_string2 = strcat('corrP=', num2str(corr_stat.P2), ' corrK=', num2str(corr_stat.K2), ' corrS=', num2str(corr_stat.S2) );
title(title_string2);

xlim(rg);
ylim(rg);

cd(path.save);
name_tag = strcat(datatag1, '_', datatag2);
savename = strcat(name_tag, '_', dataset_name{DB}, '_', 'type_', num2str(mtype) );
saveas(h1, savename, 'epsc');

end

close all;


%=========================================================================

function [output] = shuffle_data(input, col)

shuffle_ind = randperm(size(input,1))';

output = input;

for r = 1:size(output,1)
    output(r,col) = input(shuffle_ind(r), col);
end

end

