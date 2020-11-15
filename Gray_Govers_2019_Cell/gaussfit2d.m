function [sfit, gof] = gaussfit2d(img)
%{
-About-
fits a 2D Gaussian function with a single sigma to an input 2D matrix. 
A+k*exp(-(x-x0).^2/(2*sigma^2)-(y-y0).^2/(2*sigma^2))

-Inputs-
img: a 2D matrix

-varargin-
na

-Outputs-
sfit and gof: both of the outputs from matlab's builtin fit function

-Example-
   
-Supplementary-

-Keywords-

-Dependencies-
the matlab builtin fit function

-References-

-Author-
Brad Parry, 2014 July 28
%}

[r,c] = size(img);
[X,Y] = meshgrid(1:r,1:c);
X = X(:);
Y=Y(:);
Z = img(:);
[y0, x0] = find(img == max(img(:)));
%use sigma=1 as an initial estimate, this is ~near the diffraction limit
%for most pixel sizes on popular microscope cameras
sigma = 1; 
%estimate amplitude by the max
k = max(img(:));
%estimate im offset as the average of the input image
A = mean(img(:));

gauss2d = fittype('A+k*exp(-(x-x0).^2/(2*sigma^2)-(y-y0).^2/(2*sigma^2))','independent', {'x', 'y'},'dependent', 'z' );
[sfit, gof] = fit([X,Y],double(Z),gauss2d,'StartPoint',[A,k, sigma, mean(x0), mean(y0)],'Lower',[0,0,0,-inf,-inf],'Upper',[inf, inf, inf, inf]);

end

