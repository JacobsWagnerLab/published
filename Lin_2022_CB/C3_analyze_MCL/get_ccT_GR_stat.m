
function [ccT_stat, GR_stat] = get_ccT_GR_stat(lineage_data)

% Statistics on ccT and GR on MCL

num = length(lineage_data);

ccT_total = [];
GR_total = [];

ccT_per_MCL = NaN(num,2);  % mean, std
GR_per_MCL = NaN(num,2);   % mean, std


for m = 1:num
    
    % Per cell cycle averaging
    temp_ccT = lineage_data{m}.MCstat.cc_data(2:end-1,3);  % ccT; cell cycle time
    temp_GR = lineage_data{m}.MCstat.cc_data(2:end-1,4);   % GR; cell cycle growth rate    
    
    ccT_total = [ccT_total; temp_ccT];
    GR_total = [GR_total; temp_GR];
    
    % Per MCL averaging
    
    ccT_per_MCL(m,:) = [lineage_data{m}.MCstat.mean_ccT lineage_data{m}.MCstat.std_ccT];
    GR_per_MCL(m,:) = [lineage_data{m}.MCstat.mean_GR lineage_data{m}.MCstat.std_GR];
    
end

ccT_stat.num = length(ccT_total);
ccT_stat.mean = mean(ccT_total);
ccT_stat.std = std(ccT_total);

GR_stat.num = length(GR_total);
GR_stat.mean = mean(GR_total);
GR_stat.std = std(GR_total);

ccT_stat.per_MCL_mean = mean(ccT_per_MCL(:,1));
ccT_stat.per_MCL_std = mean(ccT_per_MCL(:,2));

GR_stat.per_MCL_mean = mean(GR_per_MCL(:,1));
GR_stat.per_MCL_std = mean(GR_per_MCL(:,2));


