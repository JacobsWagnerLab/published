
% Perform sibling statistics (correlation plot) and compare with shuffled
% statistics. Manaully load 'ATP_stat_type.mat'

path = {};
path.save = '/Users/wei-hsiang/Desktop/uF_local/Tree_data/parent_daughter_stat_corr/';

dataset_name = {'M9GlcCA_lac', 'M9GlcCA_ldhA', 'M9GlcCA_adhE', 'M9GlcCA_pta', 'M9GlcCA_ackA', ...
                'M9Glc_lac', 'M9Xyl_lac', 'M9Gly_lac'};

% -----------------------------------------------------------------------            
DB_num = 8;
mtype = 1;

datatag1 = 'Rsc_stat';
datatag2 = 'std';

for DB = 1:DB_num
    
data0 = triad_stat{DB}.data.(datatag1).(datatag2){mtype};
dataA = [data0(:,1) data0(:,2)]; 
dataASH = shuffle_data( dataA(:,1:2), 2 );

dataB = [data0(:,1) data0(:,3)]; 
dataBSH = shuffle_data( dataB(:,1:2), 2 );

% Calculate correlation coefficient
corr_statA = {};
corr_statA.P = corr(dataA(:,1), dataA(:,2), 'type','Pearson');
corr_statA.K = corr(dataA(:,1), dataA(:,2), 'type','Kendall');
corr_statA.S = corr(dataA(:,1), dataA(:,2), 'type','Spearman');
corr_statA.P2 = corr(dataASH(:,1), dataASH(:,2), 'type','Pearson');
corr_statA.K2 = corr(dataASH(:,1), dataASH(:,2), 'type','Kendall');
corr_statA.S2 = corr(dataASH(:,1), dataASH(:,2), 'type','Spearman');

corr_statB = {};
corr_statB.P = corr(dataB(:,1), dataB(:,2), 'type','Pearson');
corr_statB.K = corr(dataB(:,1), dataB(:,2), 'type','Kendall');
corr_statB.S = corr(dataB(:,1), dataB(:,2), 'type','Spearman');
corr_statB.P2 = corr(dataBSH(:,1), dataBSH(:,2), 'type','Pearson');
corr_statB.K2 = corr(dataBSH(:,1), dataBSH(:,2), 'type','Kendall');
corr_statB.S2 = corr(dataBSH(:,1), dataBSH(:,2), 'type','Spearman');

% -----------------------------------------------------------------------

rg = [0 1.2];

h1 = figure('position', [1 1 600 250]); 

subplot(121);
plot (dataA(:,1), dataA(:,2), 'b.'); hold on;
plot (dataB(:,1), dataB(:,2), 'g.'); hold off;

title_stringA1 = strcat('cP=', num2str(corr_statA.P), ' cK=', num2str(corr_statA.K), ' cS=', num2str(corr_statA.S) );
title_stringB1 = strcat('cP=', num2str(corr_statB.P), ' cK=', num2str(corr_statB.K), ' cS=', num2str(corr_statB.S) );

title( {title_stringA1, title_stringB1});

xlim(rg);
ylim(rg);

subplot(122);
plot (dataASH(:,1), dataASH(:,2), 'b.'); hold on;
plot (dataBSH(:,1), dataBSH(:,2), 'g.'); hold off;

title_stringA2 = strcat('cP=', num2str(corr_statA.P2), ' cK=', num2str(corr_statA.K2), ' cS=', num2str(corr_statA.S2) );
title_stringB2 = strcat('cP=', num2str(corr_statB.P2), ' cK=', num2str(corr_statB.K2), ' cS=', num2str(corr_statB.S2) );

title( {title_stringA2, title_stringB2});

xlim(rg);
ylim(rg);

cd(path.save);
name_tag = strcat(datatag1, '_', datatag2);
savename = strcat(name_tag, '_', dataset_name{DB}, '_', 'mtype_', num2str(mtype) );
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

