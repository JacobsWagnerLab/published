
% [Input]: 
% x (n-by-1 array), value will be binned
% y (n-by-1 array)
% bin_size: size of x_data
% bin_endpoints: range for binning
% min_data_point: minimal data point used in the bin, otherwise return NaN

% [Output]: 
% data: a structure with following fields
% [bin_index] : index for this bin
% [x]         : binned x axis 
% [y.num]     : number of data point in this bin
% [y.mean]    : data mean in this bin
% [y.std]     : data std in this bin
% [y.se]      : data se in this bin

function [data] = x_binning_V5(x, y, bin_size, bin_endpoints, min_data_point)

data = {};

% (1) Binning data for x-axis
data.bin_index = floor(x/bin_size)+1;  % binning index

% (2) Compute mean and std of each x-bin

bin_min = bin_endpoints(1)/bin_size;
bin_max = bin_endpoints(2)/bin_size;
bin_num = bin_max - bin_min + 1;

data.y.num = zeros(bin_num, 1); 
data.y.mean = zeros(bin_num, 1); 
data.y.std = zeros(bin_num, 1); 
data.y.se = zeros(bin_num, 1); 

for j = bin_min:bin_max
    
    temp = y;
    temp( (data.bin_index ~= j) ,:) = [];  % Remove temporal data points that not in this bin
    
    ind = j - bin_min + 1;
    
    if (size(temp,1) > min_data_point)
    
        data.y.num(ind) = size(temp,1);
        data.y.mean(ind) = nanmean(temp);
        data.y.std(ind) = nanstd(temp);
        
    else
               
        data.y.num(ind) = NaN;
        data.y.mean(ind) = NaN;
        data.y.std(ind) = NaN;
        
    end
    
end

data.x = bin_size *(bin_min:bin_max)' - (bin_size/2);
data.y.se = data.y.std ./ sqrt(data.y.num);


