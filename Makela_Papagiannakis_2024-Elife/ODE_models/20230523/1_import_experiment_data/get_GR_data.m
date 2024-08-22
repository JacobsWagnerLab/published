
function [GR_data] = get_GR_data()

GR_data = {};
addpath("./GR_data/");

file_name = {};
file_name{1} = 'CRISPRi_FtsZ_data.mat';
file_name{2} = 'CRISPRi_oriC_data.mat';
file_name{3} = 'M9gly_CRISPRi_WT_data.mat';
file_name{4} = 'M9gly_CRISPRi_oriC_data.mat';
file_name{5} = 'M9ala_CRISPRi_WT_data.mat';
file_name{6} = 'M9ala_CRISPRi_oriC_data.mat';

% Merging statistics.

GR_data = {};

for m = 1:6
    GR_data{m}.name = file_name{m};
end

for m = 1:6
    GR_data{m}.rawdata = combine_data(file_name{m});
end

% Binning statistics.

bin_size = 0.5;       % um^2
bin_range = [1 18];   % um^2

for m = 1:6
    GR_data{m}.bindata = binning(GR_data{m}.rawdata, bin_size, bin_range);
end

end


%{
%-------------------------------------------------------------------------

figure;

subplot(131);

for m = 1:2
    plot(GR_data{m}.rawdata(:,1), GR_data{m}.rawdata(:,2), '.'); hold on; 
end

subplot(132);

for m = 3:4
    plot(GR_data{m}.rawdata(:,1), GR_data{m}.rawdata(:,2), '.'); hold on; 
end

subplot(133);

for m = 5:6
    plot(GR_data{m}.rawdata(:,1), GR_data{m}.rawdata(:,2), '.'); hold on; 
end

%-------------------------------------------------------------------------

YL = [0 0.1];

figure;

subplot(131);

for m = 1:2
    errorbar(GR_data{m}.bindata(:,2), GR_data{m}.bindata(:,3), GR_data{m}.bindata(:,4), '.'); hold on; 
end
ylim(YL);

subplot(132);

for m = 3:4
    errorbar(GR_data{m}.bindata(:,2), GR_data{m}.bindata(:,3), GR_data{m}.bindata(:,4), '.'); hold on; 
end
ylim(YL);

subplot(133);

for m = 5:6
    errorbar(GR_data{m}.bindata(:,2), GR_data{m}.bindata(:,3), GR_data{m}.bindata(:,4), '.'); hold on; 
end
ylim(YL);
%}


% ======================================================================

function [stat] = binning(data, bin_size, bin_range)

% Create bin edges

bin_edge = bin_range(1):bin_size:bin_range(2);
bin_num = length(bin_edge)-1;

% Bin statistics

stat = NaN(bin_num, 4);

for j = 1:bin_num
    
    % filter the data by area
    
    area_min = bin_edge(j);
    area_max = bin_edge(j+1);
    
    temp = data;
    temp( (temp(:,1) < area_min) , : ) = [];
    temp( (temp(:,1) > area_max) , : ) = [];
    
    stat(j,1) = size(temp,1);
    stat(j,2) = (area_min + area_max)/2;
    stat(j,3) = mean(temp(:,2));
    stat(j,4) = std(temp(:,2));
    %stat(j,5) = mean(temp(:,3));
    %stat(j,6) = std(temp(:,3));

end


end

% ======================================================================

function [temp_table] = combine_data(file_name)

temp = load(file_name);

temp_table = [];  
% col1: cell size
% col2: step size
% col3: rel step size

max_num = 10000;

for c = 1:max_num
    
    temp1 = temp.cell_size{c};
    temp2 = temp.step_size{c};
    temp3 = temp.rel_step_size{c};

    temp4 = c* ones(1,length(temp1));
    
    write = [temp1' temp2' temp3' temp4'];
    
    temp_table = [temp_table; write];
    
end


end

% ======================================================================

