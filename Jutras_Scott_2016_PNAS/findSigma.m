
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function findSigma(File1)
%@author: Brad Parry
%@date: March 30, 2015
%==========================================================================
%************************Output**********************:
%sigmas_heights:        Cell arrays containing (1) relative Z position of spots
%                       tracked through a Z series (heights) and (2) the
%                       Gaussian fitted sigmas of these spots. Adjust the
%                       file name by altering the final line of this code.
%************************Input**********************:
%file1:                 The first image files in a z series
%==========================================================================
%This function tracks spots through a z series and outputs cell arrays
%containing information about each spot, tracked through all frames. It
%requires access to the files: filterIM.m, g2d.m, and gaussfit2d.m. Make
%sure that your images are saved to a single folder and that they have a
%similar naming structure, such that all images in a Z series have the same
%name, but end in a number (0-9).tif.
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 

function findSigma(File1)

% this script requires access to filterIM, g2d and gaussfit2d. 

%the name of the first image in a series of images. so long as the names
%end in a number, i.e., (0-9).tif this script will be able to find all of
%the subsequent images

% the number of frames above and below the frame with greatest intensity
% that will analyzed
frameOffset = [7, 7];

%bounds for spot - must have at least spotSize(1) pixels, but fewer than
%spotSize(2) pixels
spotSize = [5, 25];

%just an anonymous function helpful for displaying data
isho = @(x) imshow(x,[min(x(:)) max(x(:))],'InitialMagnification','fit');

%regular expressions to figure out how the user decided to name all of the
%images
[fluo_file1,ff_ext] = strtok(File1,'.');
ind = regexp(fluo_file1,'\d+$');
first_ind = str2num(fluo_file1(ind:end));
fff=['%0',num2str(length(ind:length(fluo_file1))),'d'];

%the next loop goes through all of the images it can find and searches for
%the brightest image -- this *should* approximately correspond to the frame
%that is the most in-focus. 
maxIntensities = [];
Frame = first_ind;
disp('identifying the brightest frame')
while 1
%     read all sequential files with the same prefix
    filename = [fluo_file1(1:ind-1), num2str(Frame + first_ind - 1,fff), ff_ext];
    if exist(filename,'file')
        image = imread(filename);
        %build a vector of image intensities, max only
        maxIntensities(end+1) = max(image(:));
    else
        break
    end
    Frame = Frame + 1;
end

% take the brightest frame and try to identify spots. 
disp('looking for spots in the brightest frame')
[~,maxFrame] = max(maxIntensities);
% load the brightest image
fff=['%0',num2str(length(ind:length(fluo_file1))),'d'];
brightFrame = [fluo_file1(1:ind-1), num2str(maxFrame + first_ind - 1,fff), ff_ext];
brightImage = imread(brightFrame);

%use a band-pass gaussian filter to identify spots (filterIM). depending on the image
%set, the argument ithresh may require some tweaking. ithresh modulates
%variance based thresholding as implemented from Otsu's method. values less
%than 1 make the thresholding less stringent and more than one make the
%thresholding more stringent. The bandpass filter is pretty good at
%removing high-frequency (i.e., pixel-based) noise, so I'm using a weak threshold
%here. If you want to see how the filter works, right-click on filterIM to
%open it; if you ran the addpath line above, matlab will find and open
%filterIM
fim = filterIM(brightImage,'ithresh',.3);

%If you want to view the result of filterIM, uncomment and run the
%following line:
figure(1),isho(fim)

%put a hard threshold on the filtered image...
fim(fim>0) = 1;
%If you want to view the result of thresholding, uncomment and run the
%following line:
figure(2),isho(fim)

bwc = bwconncomp(fim); %identify blobs of signal
bwx = bwc.PixelIdxList; %extract the pixels in each blob
lll = cellfun(@length, bwx);
%remove all signal blobs that are outside of the pixel range determined by 
%spotSize
bwx((lll < spotSize(1)) | (lll > spotSize(2))) = []; 

%extrach each blob from the original image, measure its centroid with a
%center of mass based approach (g2d) and save these in the centers array.
%Later on, the centers array will be used to guess where each spot is in
%all of the other z-slices
for k = 1:length(bwx)
    [r,c] = ind2sub(bwc.ImageSize,bwx{k});
    xy = g2d(brightImage(min(r)-1:max(r)+1,min(c)-1:max(c)+1));
    centers(k,1:2) = xy + [min(c)-2, min(r)-2];
end

%construct regions of interest (ROI) around each spot in the centers array. this
%routine determines ROIs that are of size 2*sz, 2*sz
sz = 10;
addBox = [-sz, sz, -sz, sz]; %set the ROI offsets in row and column space
%apply ROI offsets (addBox) to the array of spot centers
dims = round([centers(:,2), centers(:,2), centers(:,1), centers(:,1)]) + repmat(addBox, [size(centers,1),1]);
%it is possible that some ROIs were created that are outside of the image
%(e.g., by spots that were identified very close to the edge of the image)
%this line finds bad ROIs quickly...sorry, it is an ugly expression
kill = ~((~(dims(:,1) < 1)) .* (~(dims(:,3) < 1)) .* (~(dims(:,2) > size(brightImage,1))) .* (~(dims(:,4) > size(brightImage,2))));
%and now we kill the bad ROIs
centers(kill,:) = [];

%2d gauss fn for re-evaluating the fits; it is used in the for loop below
fn = @(A,k,sigma,x,x0,y,y0) A+k.*exp(-(x-x0).^2/(2*sigma^2)-(y-y0).^2/(2*sigma^2));

%figure out the frame indices offset from the brightest frame
getFrames = (maxFrame - frameOffset(1)):(maxFrame + frameOffset(2));
%initialize an array to hold all of the [sigma, height] values for all frames
%height is spot amplitude and can be used to get rid of spots that are
%multiple beads
%sigmas will be a matrix arranged like so:
%     spot1,  spot2,  spot3,  ... spotN
% Z1 
% Z2
% Z3
% .
% .
% .
% ZN
sigmas = zeros(length(getFrames), size(centers,1));
heights = sigmas;
for Frame = getFrames
    %read in the image...
    fff=['%0',num2str(length(ind:length(fluo_file1))),'d'];
    filename = [fluo_file1(1:ind-1), num2str(Frame + first_ind - 1,fff), ff_ext];
    image = imread(filename);
    
    %iterate through all ROIs...
    for k = 1:size(centers,1)
        %dtermine the image rows and columns to extract. We have to round
        %because pixels are integers and values in center are floats
        b = round([centers(k,2), centers(k,2), centers(k,1), centers(k,1)]) + addBox;
        %extract one ROI
        croppedIM = double(image(b(1):b(2),b(3):b(4)));
        %perform 2d Gaussian fitting
        [sfit, gof] = gaussfit2d(croppedIM);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %these next 9 lines are dumb - they just save fitting data and evaluate the fit. i will
        %not describe them. they can be safely commented out.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fdat{k}{Frame-min(getFrames)+1}.f = sfit;
        fdat{k}{Frame-min(getFrames)+1}.im = croppedIM;
        A = fdat{k}{Frame-min(getFrames)+1}.f.A;
        kkk = fdat{k}{Frame-min(getFrames)+1}.f.k;
        SIG = fdat{k}{Frame-min(getFrames)+1}.f.sigma;
        x0 = fdat{k}{Frame-min(getFrames)+1}.f.x0;
        y0 = fdat{k}{Frame-min(getFrames)+1}.f.y0;
        [mx,my] = meshgrid(1:size(croppedIM,1), 1:size(croppedIM,2));
        fdat{k}{Frame-min(getFrames)+1}.fim = fn(A,kkk,SIG,mx,x0,my,y0);      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %back to real code again (below)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %extract the measured sigma and height values
        sigmas(Frame-min(getFrames)+1, k) = sfit.sigma;
        heights(Frame-min(getFrames)+1, k) = sfit.A;
        %give some feedback to the user....
        disp([num2str(Frame),':   ',num2str(sfit.A),' ',num2str(sfit.sigma)])
        
    end
end

% %uncomment these lines to save the data
save('sigmas_heights_fdat','sigmas','heights','fdat')
save('sigmas_heights','sigmas','heights')

end
