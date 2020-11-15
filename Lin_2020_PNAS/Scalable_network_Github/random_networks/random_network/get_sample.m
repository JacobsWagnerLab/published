
% This script generate a sampling distribution (sample number = N)
% from a cumulative distribution funtion (cdf)
% Use histogram plot for the output to visualize the sample distribution

function [sample] = get_sample(cdf, N)

sample = NaN(N,1);

for t = 1:N

    rv = rand(1,1);  % generate an uniformly-distributed random variable on [0,1]

    for j = 1:size(cdf,1)
        
        if ( rv <= cdf(j) )
            sample(t) = j;
            break;
        else
            sample(t) = size(cdf,1);
        end    
        
    end
end

% End of the script

