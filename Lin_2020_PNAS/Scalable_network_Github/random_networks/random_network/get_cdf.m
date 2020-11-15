
% This script computes cumulative distribution funnction (cdf)
% based on probability distribution funcion (pdf).

% Input: an (n-by-1) vector. Should summing up to 1. 
% Ouput: cumulative vector, with cdf(n) = 1.

% Example: 
% - Input: pdf = [0.2 0.4 0.4]' 
% - Ouput: cdf = [0.2 0.6 1.0]';

function [cdf] = get_cdf(pdf)
    
    cdf = zeros(size(pdf,1),1);
    cdf(1,1) = pdf(1,1);

    for j = 2:size(pdf,1)
    
        cdf(j,1) = cdf(j-1,1) + pdf(j,1);
    
    end
    
end

% End for the script

