function [sfit, gof] = gaussfit2d(img)
%Fit a 2D Gaussian surface to the matrix img
%sfit and goodness of fit metrics, as constructed by matlab, as returned
%
% Brad Parry, Christine Jacobs-Wagner lab, Yale University
[r,c] = size(img);
[X,Y] = meshgrid(1:r,1:c);
X = X(:);
Y=Y(:);
Z = img(:);
[y0, x0] = find(img == max(img(:)));
sigma = 1; 
k = max(img(:));
A = mean(img(:));

gauss2d = fittype('A+k*exp(-(x-x0).^2/(2*sigma^2)-(y-y0).^2/(2*sigma^2))','independent', {'x', 'y'},'dependent', 'z' );
[sfit, gof] = fit([X,Y],double(Z),gauss2d,'StartPoint',[A,k, sigma, mean(x0), mean(y0)],'Lower',[0,0,0,-inf,-inf],'Upper',[inf, inf, inf, inf]);
 