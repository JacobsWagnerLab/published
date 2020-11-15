
% This script provide new indicies for the new variables. 

% For example,: when a new flux function, which having n=3 parameters, 
% is defined, we need to assign 3 new indicies for these parameters. 
% Suppose we already have parmaeter p(1) to p(pcount), with pcount = 20,
% then we would like to assign the new indicies as p(21), p(22), p(23) 
% to the new parameters. Hence, the output 'pind' is the array [21 22 23];

% Input:
% - n: numbers of new indicies required.

% Ouput: 
% - pind: an array of length n, with n indicies inside. 

% ==================================================================== 

function [pind] = GetNewInd(n)

global pcount;  % pcount: current index for parameters 

% Assign n new indecies to the n new parameters 
pind = pcount + [1:n];

% Increase the current index by n
pcount = pcount + n;

