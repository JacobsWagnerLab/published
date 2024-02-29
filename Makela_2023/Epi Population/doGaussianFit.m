function [gaussianApproximation, parameters] = doGaussianFit(imageToFit,meshX,meshY,parameters)
% taken from supersegger function intFindFociCurve.m for testing spot
% fitting

%parameters to be fitted are : 
%parameter(1) - Sub-pixel resolution of foci position X
%parameter(2) - Sub-pixel resolution of foci position Y
%parameter(3) - Intensity of the gaussian
%parameter(4) - sigma of gaussian
%parameter(5) - background intensity


options =  optimset('MaxIter', 1000, 'Display', 'off', 'TolX', 1/10);
% optimize parameters
[parameters] = fminsearch( @doFit, parameters, options);

% calculate final gaussian function for later use
gaussianApproximation = makeGaussianTestImage(meshX, meshY, parameters(1), parameters(2), parameters(3), parameters(5), parameters(4));


function error = doFit(parameters)
    % doFit : does the gaussian fit to the foci and calculates the error
    gaussian = makeGaussianTestImage(meshX, meshY, parameters(1), parameters(2), parameters(3), parameters(5), parameters(4));
    % difference image between gaussian and image
    tempImage = (double(imageToFit) - gaussian);
    % summed squared error between gaussian fit and image
    error = sum(sum(tempImage.^2));
end
% makes a Gaussian function from parameters
function testImage = makeGaussianTestImage(meshX, meshY, fociX, fociY, gaussianIntensity, backgroundIntensity, sigmaValue)
    testImage = backgroundIntensity + gaussianIntensity * exp( -((meshX - fociX).^2 + (meshY - fociY).^2)/(2 * sigmaValue^2) );
end

end
