function [wx, wy] = gaussianWeightX(x, y, sigma)
%{
smooth data x,y by weighting each point in x from a normal distribution of
shape deterimined by sigma.

for all x in X calculate
W = G(x-X,sigma * std(X))
where G(a,b) is the Gaussian function at points a with a standard deviation
of sigma and W is probability, or in the case of this function, weights.
Create a


for an example, run the script with no inputs.

Brad Parry, Christine Jacobs-Wagner lab; April 2016
%}

if nargin == 0
    X = 0:0.005:7;
    Y = sin(X) + exp(.3*X) + 1.5*randn(1,length(X));
    [wx, wy] = gaussianWeightX(X,Y,.2);
    figure
    plot(X,Y,'.','linewidth',2)
    hold on
    plot(wx,wy,'r','linewidth',2)
    return
end

[~,ix] = sort(x);
x = x(ix);
y = y(ix);

tile = @(x) repmat(x(:), [1, length(x)]);
%make sure the data is the correct shape
x = x(:);
y = y(:);

%normalize data for use with distance calculations
Xn = x / std(x);

dx = (repmat(Xn(:), [1, length(x)]) - repmat(Xn(:)', [length(x), 1])).^2;

%calculate weights
weights = exp( -dx ./ (2*sigma^2) );

%calculate weighted averages
wx = sum(weights.*tile(x)) ./ sum(weights);
wy = sum(weights.*tile(y)) ./ sum(weights);
end