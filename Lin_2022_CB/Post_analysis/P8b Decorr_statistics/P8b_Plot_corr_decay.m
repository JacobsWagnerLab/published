
% Manually load "corr_data_ensemble"
% This scripts plot the correlation decay curve between sibling cells.

frame_interval = 6;
mtype = 1;
corr_type = 1;  % Pearson
Tmax = 60;

% =======================================================================

index_list = [1];

figure;

subplot(121);

for p = 1:length(index_list)
    
    DB = index_list(p);
    data0 = corr_data_ensemble{DB}{mtype}.corr_data;
    
    data = get_CI_decorrelation(data0, Tmax);   % CI95% data    
    time = frame_interval * (1:length(data0))';        
    
    [output] = exp_decay_fit(data(:,2), time);        
    datafit = output.yfit;
    
    plot(time, data(:,2), 'o', 'color', 0*[1 1 1]); hold on;    
    plot(time, datafit, '-', 'color', 0*[1 1 1]); hold on;     
    
    plot(time, data(:,1), '-.', 'color', 0.5*[1 1 1]); hold on;
    plot(time, data(:,3), '-.', 'color', 0.5*[1 1 1]); hold on;
    
    xlim([0 150]);
    ylim([0.4 1]);
    
end
hold off;

% Plot the shuffled data

subplot(122);


for p = 1:length(index_list)
    
    DB = index_list(p);
    data0 = corr_data_ensemble{DB}{mtype}.corr_data;
    data_SH = shuffle_data(data0);
    
    data = get_CI_decorrelation(data_SH, Tmax); % CI95% data    
    time = frame_interval * (1:length(data0))';        
    
    [output] = exp_decay_fit(data(:,2), time);        
    datafit = output.yfit;
    
    plot(time, data(:,2), 'o', 'color', 0*[1 1 1]); hold on;    
    plot(time, datafit, '-', 'color', 0*[1 1 1]); hold on;     
    
    plot(time, data(:,1), '-.', 'color', 0.5*[1 1 1]); hold on;
    plot(time, data(:,3), '-.', 'color', 0.5*[1 1 1]); hold on;
    
    xlim([0 150]);
    %ylim([0.4 1]);
    
end
hold off;

% =========================================================================

function [output] = exp_decay_fit(data, time)

modelfun = @(b,x) b(1)* ( exp(-b(2)*x) ) + (1-b(1));

beta0 = [0.5 0.02];

x = time;
y = data;

opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
beta = nlinfit(x,y, modelfun, beta0, opts);

yfit = beta(1) * ( exp(-beta(2)*x) ) + (1-beta(1));

output = {};
output.beta = beta;
output.x = x;
output.y = y;
output.yfit = yfit;

end


% ======================================================================

function [corr_data_SH] = shuffle_data(corr_data0)

corr_data_SH = {};

L = length(corr_data0);

for k = 1:L
    
    temp = corr_data0{k};
    
    if (size(temp,1) > 1)
        corr_data_SH{k} = shuffle_single_time(temp);
    else
        corr_data_SH{k} = temp;
    end
    
end

end


% ======================================================================

function [temp_permuted] = shuffle_single_time(temp)

temp_permuted = NaN*temp;

    n = size(temp,1);
    indP = randperm(n)';

    for j = 1:n
        temp_permuted(j,:) = [temp(j,1) temp(indP(j),2)];
    end

end

% ======================================================================


