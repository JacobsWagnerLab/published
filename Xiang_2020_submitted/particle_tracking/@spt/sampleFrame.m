%% Sample the 1st frame to help determine the object intensity and size.
%
% -About-
%   The method sampleFrame() samples the first frame of the image stack. It
%   helps the user to determine the appropriate intensity and size range
%   for thresholding in later steps, so that only the particles are
%   localized. It generates a 1x3 subplots, which shows:
%       1. the distribution of spot sizes detected after filtering and
%          thresholding
%       2. the band-pass filtered image
%       3. a binary copy of the filtered image
% 
% -Input-
%   - obj: spt object
%   - intensityThreshold: the intensity threshold below which any pixels
%                         are zeroed after the filtering
%
% -Example-
%   % Sample the 1st frame with an intensity threshold of 10
%   myParticle.sampleFrame(10);
% 
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

function sampleFrame(obj,intensityThreshold)

% If no input for the intensity threshold, raise the error for the user.
if nargin < 2
    error('No intensity threshold found.');
end

% Read in the first frame of the image stack
im = double(imread(obj.stackPath, 'Index', 1));
% Filter with a band-pass filter
bandPassedIm = bpass(im,1,obj.estimatedSize,intensityThreshold);
% Create a binary copy of the filtered image
bwIm = bandPassedIm;
bwIm(bwIm>0) = 1;

% Print the number of spots detected and their sizes
bwStruct = bwconncomp(bandPassedIm);
spotNum = bwStruct.NumObjects;
spotSize = zeros(spotNum,1);
for ii = 1:spotNum
    spotSize(ii) = length(bwStruct.PixelIdxList{ii}(:));
end
fprintf('%d spots detected. Size range: %d ~ %d pixels (mean = %.2f).\n',...
         spotNum,min(spotSize),max(spotSize),mean(spotSize));

% Plot the spot size distribution
subplot(131);
histogram(spotSize,size(spotSize,1),'facecolor','k');
xlabel('Spot size [px]');
ylabel('Counts')
title('Spot size distribution','fontsize',14)
axis square;

% Plot the band-passed image
subplot(132);
imshow(bandPassedIm,[min(bandPassedIm(:)), max(bandPassedIm(:))]);
title('Band-pass filtered image','fontsize',14);

% Plot the binary copy
subplot(133);
imshow(bwIm,[min(bwIm(:)), max(bwIm(:))]);
title('Binary image','fontsize',14);

end
