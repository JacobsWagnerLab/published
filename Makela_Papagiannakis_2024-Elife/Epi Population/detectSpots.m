function [locX, locY, fitSigma, bgsubtrInt, fitScore, score] = ...
    detectSpots(im, windowSize, cellMask)
% Detects spots from an image by fitting a Gaussian.
% Image is filtered by bandpass filter and local maxima are found. 
% X-score threshold is used to remove maxima that are not above the
% normalized threshold in raw pixel values. 
% After this remaining maxima are fitted by a Gaussian function and 
% parameters for each focus are saved.
% NOTE! same median cellular background intensity is subtracted from all 
% spots to avoid problems. 

% input:
%   im - image matrix
%	windowSize - window size for filtering and fitting ROI
%   cellMask - a binary mask for cell area used for background estimation

% output:
%   locX, locY - sub-pixel coordinates of a spot
%   fitSigma - std of the gaussian fit
%   bgsubtrInt - background subtracted intensity of the filtered focus
%   fitScore - score of gaussian fit to the focus. A measure of the overlap 
%               of the two distributions, the focus and the gaussian fit.
%   score - how good of a focus it is. It combines intensity and foci 
%               fitting (Intensity * fitScore / * fitSigma)

% 5/28/21 Jarno Makela

% bandpass window size
bpasskerneldiametre = max(round(windowSize),1);

% filtered image by bandpass
imFilt = bpass(im,1,bpasskerneldiametre);

% find maxima from filtered image; im is used for threshold
pos = findLocalMax(imFilt,cellMask);

% only if local maxima found
if ~isempty(pos)
    % initialize vectors
    locX = nan(size(pos,1),1);
    locY = nan(size(pos,1),1);
    fitSigma = nan(size(pos,1),1);
    bgsubtrInt = nan(size(pos,1),1);
    fitScore = nan(size(pos,1),1);
    score = nan(size(pos,1),1);

    for jj = 1:size(pos,1) %for each spot  
        % get ROI for spots and create mesh grids
        [ROI, xStart, yStart] = getLocalizationROI(pos(jj,:),imFilt,windowSize);
        [meshX,meshY] = meshgrid(xStart:(xStart+size(ROI,1)-1), ...
            yStart:(yStart+size(ROI,1)-1));

        % for each ROI fit gaussian
        if ~isempty(ROI)
            % initialize parameters
            backgroundIntensity = 0;
            gaussianIntensity = imFilt(pos(jj,2), pos(jj,1)) - backgroundIntensity;
            sigmaValue = 1;

            % parameters to fit
            parameters(1) = pos(jj,1); % sub-pixel X position of foci 
            parameters(2) = pos(jj,2); % sub-pixel Y position of foci
            parameters(3) = gaussianIntensity; % max intensity of the gaussian
            parameters(4) = sigmaValue; % sigma of gaussian
            parameters(5) = backgroundIntensity; % background intensity

            % settings for optimizer
            options =  optimset('MaxIter', 1000, 'Display', 'off', 'TolX', 1/10);

            % optimize the fit parameters
            [parameters] = fminsearch( @doFit, parameters, options);

            % calculate the final gaussian function
            % meshX, meshY, fociX, fociY, Int, bgInt, sigma
            ROI_fitted = makeGaussianTestImage(meshX, meshY, parameters(1), ...
                parameters(2), parameters(3), parameters(5), parameters(4));

            % crop out fit gaussian from original image with threshold of
            % 0.1 of max intensity
            ROI_crop = ROI;
            ROI_crop(ROI_fitted < 0.1 * max(max(ROI_fitted))) = 0;

            % summed intensities for filtered and gaussian fit
            imageTotal = sqrt(sum(sum(ROI_crop)));
            guassianTotal = sqrt(sum(sum(ROI_fitted)));

            % measure of overlap between the two distributions
            fitQuality = sum(sum(sqrt(ROI_crop) .* sqrt(ROI_fitted))) / ...
                (imageTotal * guassianTotal);

            % save fit parameters
            locX(jj) = parameters(1);
            locY(jj) = parameters(2);
            fitSigma(jj) = parameters(4);
            bgsubtrInt(jj) = parameters(3);
            fitScore(jj) = fitQuality;
            score(jj) = parameters(3) / parameters(4) * fitQuality;
        end
    end
    % remove values outside image area from low quality fitting
    for ii = 1:length(locX)
       if locX(ii) < 0 || locX(ii) > size(im,2) || locY(ii) < 0 || locY(ii) > size(im,1)
            locX(ii) = NaN;
            locY(ii) = NaN;
            fitSigma(ii) = NaN;
            bgsubtrInt(ii) = NaN;
            fitScore(ii) = NaN;
            score(ii) = NaN; 
       end
    end

    % remove NaN
    locX = locX(~isnan(locX));
    locY = locY(~isnan(locY));
    fitSigma = fitSigma(~isnan(fitSigma));
    bgsubtrInt = bgsubtrInt(~isnan(bgsubtrInt));
    fitScore = fitScore(~isnan(fitScore));
    score = score(~isnan(score));
    
else
    locX = NaN;
    locY = NaN;
    fitSigma = NaN;
    bgsubtrInt = NaN;
    fitScore = NaN;
    score = NaN;
end



function res = bpass(image_array,lnoise,lobject)         
    % PURPOSE:      Implements a real-space bandpass filter that suppresses
    %               pixel noise and long-wavelength image variations while
    %               retaining information of a characteristic size.
    % INPUTS:
    %               image:  The two-dimensional array to be filtered.
    %               lnoise: Characteristic lengthscale of noise in pixels.
    %                       Additive noise averaged over this length should
    %                       vanish. May assume any positive floating value.
    %                       May be set to 0 or false, in which case only the
    %                       highpass "background subtraction" operation is
    %                       performed.
    %               lobject: (optional) Integer length in pixels somewhat
    %                       larger than a typical object. Can also be set to
    %                       0 or false, in which case only the lowpass
    %                       "blurring" operation defined by lnoise is done,
    %                       without the background subtraction defined by
    %                       lobject.  Defaults to false.
    %               threshold: (optional) By default, after the convolution,
    %                       any negative pixels are reset to 0.  Threshold
    %                       changes the threshhold for setting pixels to
    %                       0.  Positive values may be useful for removing
    %                       stray noise or small particles.  Alternatively, can
    %                       be set to -Inf so that no threshholding is
    %                       performed at all.
    %		      SH101104 -set this so is always zero
    %
    % OUTPUTS:
    %               res:    filtered image.
    % PROCEDURE:
    %               simple convolution yields spatial bandpass filtering.
    % NOTES:
    % Performs a bandpass by convolving with an appropriate kernel.  You can
    % think of this as a two part process.  First, a lowpassed image is
    % produced by convolving the original with a gaussian.  Next, a second
    % lowpassed image is produced by convolving the original with a boxcar
    % function. By subtracting the boxcar version from the gaussian version, we
    % are using the boxcar version to perform a highpass.
    %
    % original - lowpassed version of original => highpassed version of the
    % original
    %
    % Performing a lowpass and a highpass results in a bandpassed image.
    %
    % Converts input to double.  Be advised that commands like 'image' display
    % double precision arrays differently from UINT8 arrays.

    % MODIFICATION HISTORY:
    %               Written by David G. Grier, The University of Chicago, 2/93.
    %
    %               Greatly revised version DGG 5/95.
    %
    %               Added /field keyword JCC 12/95.
    %
    %               Memory optimizations and fixed normalization, DGG 8/99.
    %               Converted to Matlab by D.Blair 4/2004-ish
    %
    %               Fixed some bugs with conv2 to make sure the edges are
    %               removed D.B. 6/05
    %
    %               Removed inadvertent image shift ERD 6/05
    %
    %               Added threshold to output.  Now sets all pixels with
    %               negative values equal to zero.  Gets rid of ringing which
    %               was destroying sub-pixel accuracy, unless window size in
    %               cntrd was picked perfectly.  Now centrd gets sub-pixel
    %               accuracy much more robustly ERD 8/24/05
    %
    %               Refactored for clarity and converted all convolutions to
    %               use column vector kernels for speed.  Running on my
    %               macbook, the old version took ~1.3 seconds to do
    %               bpass(image_array,1,19) on a 1024 x 1024 image; this
    %               version takes roughly half that. JWM 6/07
    %
    %       This code 'bpass.pro' is copyright 1997, John C. Crocker and
    %       David G. Grier.  It should be considered 'freeware'- and may be
    %       distributed freely in its original form when properly attributed.

    if nargin < 3, lobject = false; end
    normalize = @(x) x/sum(x);
    image_array = double(image_array);

    if lnoise == 0
        gaussian_kernel = 1;
    else
        gaussian_kernel = normalize(...
            exp(-((-ceil(5*lnoise):ceil(5*lnoise))/(2*lnoise)).^2));
    end

    if lobject
        boxcar_kernel = normalize(...
            ones(1,length(-round(lobject):round(lobject))));
    end

    % JWM: Do a 2D convolution with the kernels in two steps each.  It is
    % possible to do the convolution in only one step per kernel with
    %
    % gconv = conv2(gaussian_kernel',gaussian_kernel,image_array,'same');
    % bconv = conv2(boxcar_kernel', boxcar_kernel,image_array,'same');
    %
    % but for some reason, this is slow.  The whole operation could be reduced
    % to a single step using the associative and distributive properties of
    % convolution:
    %
    % filtered = conv2(image_array,...
    %   gaussian_kernel'*gaussian_kernel - boxcar_kernel'*boxcar_kernel,...
    %   'same');
    %
    % But this is also comparatively slow (though inexplicably faster than the
    % above).  It turns out that convolving with a column vector is faster than
    % convolving with a row vector, so instead of transposing the kernel, the
    % image is transposed twice.

    gconv = conv2(image_array',gaussian_kernel','same');
    gconv = conv2(gconv',gaussian_kernel','same');

    if lobject
        bconv = conv2(image_array',boxcar_kernel','same');
        bconv = conv2(bconv',boxcar_kernel','same');

        filtered = gconv - bconv;
    else
        filtered = gconv;
    end

    % Zero out the values on the edges to signal that they're not useful.
    lzero = max(lobject,ceil(5*lnoise));
    filtered(1:(round(lzero)-1),:) = 0;
    filtered((end - lzero + 2):end,:) = 0;
    filtered(:,1:(round(lzero)-1)) = 0;
    filtered(:,(end - lzero + 2):end) = 0;

    filtered(filtered < 0) = 0;
    res = filtered;
end

function coord = findLocalMax(imFilt,cellMask)
    % finds local maxima in an image to pixel level accuracy.
    % local maxima condition is >= rather than >
    % inspired by Crocker & Griers algorithm, and also Daniel Blair & Eric 
    % Dufresne's implementation
    % input:
    %       im - input image for detection
    %       zscore - detection threshold - pixels must be strictly greater above
    %           the z-score value. Mean and Std calculated with 0.95 quantile
    % out : 
    %       x,y coordinates of local maxima
    
    %identify above threshold pixels
    [y,x] = find(imFilt>0 & cellMask);
    
    % delete pixels identified at the boundary of the image
    [imsizey, imsizex]= size(imFilt);
    edgepixel_idx = find( y==1 | y==imsizey | x==1 | x==imsizex);
    y(edgepixel_idx) = [];
    x(edgepixel_idx) = [];

    % check if pixel is a local maxima
    islocalmaxima = zeros(numel(x),1);
    for i = 1:numel(x)
        subim = imFilt([y(i)-1:y(i)+1],[x(i)-1:x(i)+1]);
        islocalmaxima(i) = all(subim(2,2)>=subim(:));
    end
    
    % assign all above threshold pixels to out initially
    coord = [x,y];
    
    % delete pixels which are not local maxima
    coord(find(~islocalmaxima),:) = [];
end

function [ROI, xstart, ystart] = getLocalizationROI(point_pos,im, windowSize)
    [sizey, sizex] = size(im);
    X0 = point_pos(1);
    Y0 = point_pos(2);
    
    % round X0, Y0 to use as matrix locations
    X0_int = round(X0);
    Y0_int = round(Y0);
    windowSize = round(windowSize); %radius should already be an integer anyway
    
    % setup the limits of the cropped image
    xstart = X0_int-windowSize;
    xfinish = X0_int+windowSize;
    ystart = Y0_int-windowSize;
    yfinish = Y0_int+windowSize;
    
    % check if any of the limits are out of bounds - if so, skip that point
    if (xstart<1) || (xstart>sizex) ||  (xfinish<1) || (xfinish>sizex) ...
            || (ystart<1) || (ystart>sizey) ||  (yfinish<1) || (yfinish>sizey)
        ROI = [];
    else
        %crop to a small area around the point
        ROI = im( ystart:yfinish, xstart:xfinish);
    end
end

% fits gaussian to the foci and calculates the error
function error = doFit(parameters)
    % gaussian fit from parameters
    gaussian = makeGaussianTestImage(meshX, meshY, parameters(1), ...
        parameters(2), parameters(3), parameters(5), parameters(4));
    
    % difference image between gaussian fit and image
    tempImage = (double(ROI) - gaussian);
    
    % summed squared error between gaussian fit and image
    error = sum(sum(tempImage.^2));
end

% generates a 2D Gaussian function from parameters
function testImage = makeGaussianTestImage(meshX, meshY, fociX, fociY, ...
        gaussianIntensity, backgroundIntensity, sigmaValue)

    testImage = backgroundIntensity + gaussianIntensity * ...
        exp(-((meshX - fociX).^2 + (meshY - fociY).^2)/(2 * sigmaValue^2));
end


end

