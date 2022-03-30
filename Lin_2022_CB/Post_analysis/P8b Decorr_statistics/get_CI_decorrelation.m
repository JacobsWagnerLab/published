
% Construct condifence interval by bootstrap method. 
% Manaully load "corr_data_ensemble"

function [CI95_data] = get_CI_decorrelation(data, Tmax)

%data = corr_data_ensemble{DB}{mtype}.corr_data; 
%Tmax = 60;

bootstrap_repeat = 1000;
sample_size = 1000;
min_data_num = 10;

CI95_data = NaN(Tmax, 3);

for t = 1:Tmax

    data_t = data{t};    
    data_num = size(data_t,1); 
    
    if (data_num > min_data_num)
        
        [CI95] = get_CI_single_timepoint(data_t, bootstrap_repeat, sample_size);
        CI95_data(t,:) = CI95;

    end
    
end

%figure; plot(CI95_data);

end


% =======================================================================

function [CI95] = get_CI_single_timepoint(data_t, bootstrap_repeat, sample_size)

bootstrap_result = NaN(bootstrap_repeat,2);

for k = 1:bootstrap_repeat
    [output] = bootstrap_one_time(data_t, sample_size);
    bootstrap_result(k,1) = output.corr_original;
    bootstrap_result(k,2) = output.corr_bootstrap;
end

SD = std(bootstrap_result(:,2));
CI95_level = 2;   % 2*SD for 95% confidence interval
CI95 = [bootstrap_result(1,1)-CI95_level*SD bootstrap_result(1,1) bootstrap_result(1,1)+CI95_level*SD];


end

% =======================================================================
% Perform bootstrap by resampling

function [output] = bootstrap_one_time(data_t, sample_size)

num = length(data_t);
ind = randi(num, sample_size, 1);
temp_t = NaN(sample_size, 2);

for j = 1:sample_size
    
    temp_t(j,:) = data_t(ind(j),:);
    
end

output = {};
output.corr_original = corr(data_t(:,1), data_t(:,2), 'type', 'Pearson');
output.corr_bootstrap = corr(temp_t(:,1), temp_t(:,2), 'type', 'Pearson');

end

