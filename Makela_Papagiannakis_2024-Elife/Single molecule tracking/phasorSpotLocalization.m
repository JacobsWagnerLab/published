function [locX, locY, magX, magY, scores] = phasorSpotLocalization(im, zscore_thres, windowSize)
% Phasor based localization of fluorescent spots.
% The algorithm converts the region of interest around a point spread function 
% to two phase vectors (phasors) by calculating the first Fourier coefficients 
% in both the x- and y-direction. The angles of these phasors are used to localize 
% the center of the fluorescent emitter, and the ratio of the magnitudes of the
% two phasors can be used for detection of assymetry in spots (e.g. astigmatism).

% input:
%   im - image matrix
%	zscore_thres - minimum z-score of spots above background from filtered image
%	windowSize - window size for filtering and phasor detection

% output:
%   locX, locY - coordinates of spots
%   magX, magY - magnitude of detected spots from fourier spectrum

% Cite: J. Chem. Phys. 148, 123311 (2018); https://doi.org/10.1063/1.5005899

% 3/26/21 Jarno Makela

% bandpass window size
bpasskerneldiametre = max(round(windowSize),1);

% filtered image
imFilt = bpass(im,1,bpasskerneldiametre);

% find maxima from filtered image
[pos, scores] = findLocalMax(imFilt,zscore_thres);

% initialize vectors
locX = nan(size(pos,1),1);
locY = locX;
magX = locX;
magY = locX;
for jj = 1:size(pos,1) %for each spot
    % get ROI for spots
    [ROI, xStart, yStart] = getLocalizationROI(pos(jj,:),im,windowSize);
    
    % do the phasor localization
    if ~isempty(ROI)
        [posX, posY, intX, intY] = phasorLocalization(ROI);
        locX(jj) = posX + xStart - 1;
        locY(jj) = posY + yStart - 1;
        magX(jj) = intX;
        magY(jj) = intY;
    end
end
locX = locX(~isnan(locX));
locY = locY(~isnan(locY));
magX = magX(~isnan(magX));
magY = magY(~isnan(magY));
scores = scores(~isnan(locY));

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

function [coord, scores] = findLocalMax(imFilt,zscore_thres)
    % finds local maxima in an image to pixel level accuracy.
    % local maxima condition is >= rather than >
    % inspired by Crocker & Griers algorithm, and also Daniel Blair & Eric Dufresne's implementation
    % input:
    %       im - input image for detection
    %       zscore - detection threshold - pixels must be strictly greater above
    %           the z-score value. Mean and Std calculated with 0.95 quantile
    % out : 
    %       x,y coordinates of local maxima
    
    % estimate mean and std of filtered image values below 0.95 quantile
    im_vector = double(imFilt(:));
    mean_im = mean(im_vector(im_vector<quantile(im_vector,0.95)));
    std_im = std(im_vector(im_vector<quantile(im_vector,0.95)));
    
    % normalize pixel values into z-score
    norm_im = (double(imFilt)-mean_im)./std_im;
    
    % find pixels above z-score
    [y,x] = find(norm_im > zscore_thres);
    scores = norm_im(norm_im > zscore_thres);

    % delete pixels identified at the boundary of the image
    [imsizey, imsizex]= size(imFilt);
    edgepixel_idx = find( y==1 | y==imsizey | x==1 | x==imsizex);
    y(edgepixel_idx) = [];
    x(edgepixel_idx) = [];
    scores(edgepixel_idx) = [];

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
    scores(find(~islocalmaxima),:) = [];
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
    if (xstart<1) || (xstart > sizex) ||  (xfinish<1) || (xfinish > sizex) ...
            || (ystart<1) || (ystart > sizey) ||  (yfinish<1) || (yfinish > sizey)
        ROI = [];
    else
        %crop to a small area around the point
        ROI = im( ystart:yfinish, xstart:xfinish);
    end
end

function [PositionX, PositionY, MagnitudeX, MagnitudeY] = phasorLocalization(ROI)
    % Perform a 2D Fourier transformation on the complete ROI.
    fft_values = fft2(ROI); 
    
    % Get the size of the matrix 
    WindowPixelSize = size(ROI,1); 
    
    % Calculate the angle of the X-phasor from the first Fourier coefficient in X
    angX = angle(fft_values(1,2));
    
    % Correct the angle
    if angX > 0 
        angX = angX-2*pi; 
    end

    % Normalize the angle by 2pi and the amount of pixels of the ROI
    PositionX = (abs(angX)/(2*pi/WindowPixelSize) + 1);
    
    % Calculate the angle of the Y-phasor from the first  Fourier coefficient in Y
    angY = angle(fft_values(2,1));
    
    % Correct the angle
    if angY>0
        angY=angY-2*pi; 
    end

    % Normalize the angle by 2pi and the amount of pixels of the ROI
    PositionY = (abs(angY)/(2*pi/WindowPixelSize) + 1);
    
    % Calculate the magnitude of the X and Y phasors by taking the absolute value of the first Fourier coefficient in X and Y
    MagnitudeX = abs(fft_values(1,2))/WindowPixelSize^2;
    MagnitudeY = abs(fft_values(2,1))/WindowPixelSize^2; 
end

end
