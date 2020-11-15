%% Fit input matrix with a bivariate Gaussian distribution (dependency).
%
% -About-
%   A fitting function that fits the input matrix with a bivariate Gaussian
%   distribution.
%
% -Input-
%   - img: numeric matrix that represents the image intensity at each pixel
%
% -Output-
%   - sfit: statistics of the fit results
%   - gof:  goodness of fit
%
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University


function [sfit, gof] = gaussfit2d(img)
[r,c] = size(img);
[X,Y] = meshgrid(1:r,1:c);
X = X(:);
Y = Y(:);
Z = img(:);

% Set up for initial guesses
% The intensity peak exists approximately at the center of the image
x0 = c/2;
y0 = r/2;
% Just some guess for the std for the intensity distribution
sigmaX = c/4;
sigmaY = r/4;
% Intensity range
A = max(img(:))- min(img(:));
% Image background
B = min(img(:));
% Fit function form
gauss2d = fittype('B+A*exp(-(x-x0).^2/(2*sigmaX^2)-(y-y0).^2/(2*sigmaY^2))','independent', {'x', 'y'},'dependent', 'z', 'coefficients',{'B','A','sigmaX','sigmaY','x0','y0'});
% Perform actual fitting and generate results
[sfit, gof] = fit([X,Y],double(Z),gauss2d,...
              'StartPoint',[B, A, sigmaX, sigmaY, x0, y0],...
              'Lower',[min(img(:)),0,sigmaX,sigmaY,0,0],...
              'Upper',[mean(img(:)), max(img(:))- min(img(:)), c/2, r/2, c, r]);
   