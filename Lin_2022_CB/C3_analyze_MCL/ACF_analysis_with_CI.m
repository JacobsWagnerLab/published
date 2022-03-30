
function [data] = ACF_analysis_with_CI(input)

% Autocorrelation with confidence interval

data = {};
data.x = input;   % This time series need to stationary

data.fft = fft(data.x);
data.power_spec = (data.fft).*conj(data.fft);

ACF = ifft(data.power_spec);
ACF = ACF/ACF(1);

L = length(data.x);
ACF_variance = NaN(L,1);

for j = 1:L
    
    temp = 0;
    
    for k = 1:j       
        temp = temp + ACF(k)^2;        
    end
    
    ACF_variance(j) = sqrt( (1 + 2*temp)/L );

end

% Define confidence interval (95%)
alpha = 0.05;
z_alpha = erfinv(1-(alpha/2)); 
CI = [ACF-z_alpha*ACF_variance ACF+z_alpha*ACF_variance];

data.ACF = ACF;
data.CI = CI;

end

