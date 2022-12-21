function cellList = findIrregularSpotsParallel(cellList, images, varargin)
%locate spots that may or may not be diffraction limited inside cells found
%in cellList in provided images. Finds local peaks and then constructs a
%shell of user determined thickness at radius away from bright pixels. 
%The intensity of the spot is measured as the ratio of the peak intensity 
%to the intensity of pixels within the shell.
%
%cellList   - cellList returned from Oufti
%images     - a cell array of image locations such that images{1} 
%             corresponds to the frame at cellList.meshData{1}
%peakRadius - a pixel will only be considered a local maxima if it is the
%             brightest peak within radius defined by peakRadius 
%edgeDist   - Radius from the peak where the outter edge of the shell will
%             be constructed
%shellThickness 
%           - thickness of the shell constructed at edgeDist to be used in 
%             intensity ratio calculation
%centerDist - Radius from the peak that will be used to calculate peak
%             intensity. Intensity of all pixels within centerDist from the
%             peak will be used in the intensity ratio calculation. Should 
%             typically be set to small values (i.e., 1)
%intensityRatioThreshold 
%           - if the ratio if peak to shell intensity is less
%             than intensityRatioThreshold, the spot will be discarded
%quantileThreshold 
%           - pixels must be greater than quantileThreshold in the
%             spot to be fit, else they will be ignored even if they are
%             within fitRadius. This permits a large fit radius to be
%             selected for irregular spots. Pixels to fit can then be
%             chosen based on intensity by modulating quantileThreshold
%fitRadius  - radius in pixels from the peak used to perform centroid 
%             calculation
%
%Brad Parry, 2015 October

%an anonymous function to convert an oufti cellmesh into a polygon
m2model = @(x) double(cat(1, x(1:end-1,1:2), flipud(x(:,3:4))));

%specify default settings in case the user specifies no additional values
% fitRadius = 3;
% edgeDist = 2.5;
% centerDist = 1;
% peakRadius = 11;
% shellThickness = 1;
% quantileThreshold = .1;
% intensityRatioThreshold = 1.195;

%SeqA-mCherry values EC M9gluCAAT
%  fitRadius = 3;
%  edgeDist = 2.5;
%  centerDist = 1;
%  peakRadius = 3;
%  shellThickness = 1;
%  quantileThreshold = .1;
%  intensityRatioThreshold = 1.08;

%SeqA-mCherry values EC M9gly
fitRadius = 2;
edgeDist = 2.5;
centerDist = 1;
peakRadius = 3;
shellThickness = 1;
quantileThreshold = .1;
intensityRatioThreshold = 1.18;

for k = 1:length(varargin)

    if strcmpi(varargin{k},'fitRadius')
        fitRadius = varargin{k+1};
    elseif strcmpi(varargin{k},'edgeDist')
        edgeDist = varargin{k+1};
    elseif strcmpi(varargin{k},'centerDist')
        centerDist = varargin{k+1};
    elseif strcmpi(varargin{k},'peakRadius')
        peakRadius = varargin{k+1};
    elseif strcmpi(varargin{k},'shellThickness')
        shellThickness = varargin{k+1};
    elseif strcmpi(varargin{k},'quantileThreshold')
        quantileThreshold = varargin{k+1};
    elseif strcmpi(varargin{k},'intensityRatioThreshold')
        intensityRatioThreshold = varargin{k+1};
    end
end

%make sure a couple parameters are integers...
edgeDist = ceil(edgeDist);
centerDist = ceil(centerDist);

%construct a mophological filtering element and prepare it for
%morphCompFilter
se = ballElement(ceil(peakRadius*2)+2, peakRadius);
center = floor(size(se,1)/2);
se = double(se);
se(se==0) = nan;
se(center+1,center+1) = nan;

%iterate through all frames present in meshData
total_meshData = cellList.meshData;
parfor F = 1:length(cellList.meshData)

%     This is for when F is a cell array of  the file locations of your
%     images
%     im = imread(images{F});

%   This is for when F is a 'stack' of your images that have already been
%   pre-loaded
    im = images(:,:,F);
    
%     in principle, this could be done outside the loop, assuming image size never changes
    imsize = size(im); 

    %perform initial estimation of peak locations
    peaks = morphCompFilter(im, se, '>');
    %initial guesses of peaks will be where peaks == 1, so get these indices
    [r,c] = find(peaks == 1);

    %find the relative locations of pixels that are distance edgeDist away from
    %a bright pixel and within a thin shell of thickness shellThickness
    %determine the pixels that are within the shell...
    e = ballElement(edgeDist*2+1,edgeDist);
    e1 = ballElement(edgeDist*2+1,edgeDist-shellThickness);
    edgePixels = e - e1;
    %indicies of pixels located within the shell 
    [pixelShellR,pixelShellC] = find(edgePixels == 1);
    if isempty(pixelShellR), error('shellThickness is too small, increase to a positive integer'), end
    %determine relative locations
    pixelShellR = pixelShellR - (edgeDist + 1);
    pixelShellC = pixelShellC - (edgeDist + 1);

    %determine the relative locations of pixels to be used in center
    %calculation: centerDist away from the central peak pixel
    e = ballElement(ceil(centerDist)*2 + 1, centerDist);
    [centralR,centralC] = find(e == 1);
    %determine relative locations
    centralR = centralR - (centerDist + 1);
    centralC = centralC - (centerDist + 1);

    %do an initial sort through the peaks and eliminate the ones that are
    %not at least as bright as the intensityRatioThreshold
    intensityRatio = nan(1,length(r));
    for k = 1:length(r)
        %convert subscripts to index for current peak after applying
        %offsets specified by pixelShellR and pixelShellC
        shellInds = (r(k) + pixelShellR) + imsize(1)*((c(k) + pixelShellC) - 1);
        centralInds = (r(k) + centralR) + imsize(1)*((c(k) + centralC) - 1);
    
        %compare the intensity in the center to the intensity in the outter
        %shell
        shellMean = mean(double(im(shellInds)));
        centralMean = mean(double(im(centralInds)));
        intensityRatio(k) = centralMean / shellMean;    
    end
    
    %spots with intensityRatios less than the threshold will be eliminated
    kill = intensityRatio < intensityRatioThreshold;
    intensityRatio(kill) = [];
    r(kill) = [];
    c(kill) = [];

    spots = [];
    pos = [];

    %This is the slowest part of the code. Now that spots have been
    %identified, they have to be assigned to cells in the cellList
    disp(['assigning spots to cells in frame ',num2str(F)])
    for C = 1:length(total_meshData{F})

        if isempty(total_meshData{F}{C}) || ~isfield(total_meshData{F}{C},'mesh') ||...
                length(total_meshData{F}{C}.mesh) < 6
            %you know the drill...
            continue
        end
        
        %convert meshes to polygons and polygons to masks and check if
        %spots are inside. This method may be a little slow, if too slow,
        %we can check if spots are inpolygon instead. I chose the mask
        %approach because converting sub-pixel polygons to masks has the
        %effect of slightly enlarging the cell area. This is beneficial for
        %the case when the cell mesh may be too tight or shifted between
        %phase and fluorescence acquisition and/or if spots are at the cell
        %boundary
        co = m2model(total_meshData{F}{C}.mesh);
        msk = poly2mask(co(:,1),co(:,2),imsize(1),imsize(2));

        %have to loop from the end of the vector to the begining because
        %items will be dropped from the vector as we go
        for k = length(r):-1:1

            %if the spot is outside of the current cell (outer loop C), skip
            %the spot
            if msk(r(k),c(k)) == 0,continue, end

            spot = double(im(r(k)-fitRadius:r(k)+fitRadius,c(k)-fitRadius:c(k)+fitRadius));
            spotOrig = spot;
            spots{end+1} = spot;
            spot = spot - quantile(spot(:),quantileThreshold);

            pos(end+1,:) = g2d(spot) +  [c(k), r(k)] - (fitRadius+1);

            if ~isfield(total_meshData{F}{C},'spotPosition')
                total_meshData{F}{C}.spotPosition = [];
                total_meshData{F}{C}.intensityRatio = [];
                total_meshData{F}{C}.spotIntensity = [];
            end

            total_meshData{F}{C}.spotPosition(end+1,:) = pos(end,:);
            total_meshData{F}{C}.intensityRatio(end+1) = intensityRatio(k);
            total_meshData{F}{C}.spotIntensity(end+1) = sum(spotOrig(:));

            %eliminate used spots
            r(k) = [];
            c(k) = [];

        end

    end
end
cellList.meshData=total_meshData;
end

function peaks = morphCompFilter(im, se, operation)
%use a morphological element to perform isotropic pixel comparison 
%according to the operation specified by operator.
%
%im = image to be filtered
%
%se = morphological element to perform comparisons, matrix of Nan and 1
%the pixel at the center. The central value (as determined by center, see
%below), must be NaN or the operator should be one of [<>]=
%all dimensions of se must be odd and have the same length
%
%center = specifies the central location of se. center*2 + 1
%, will be broadcast center = center + [0,0]
%operation for comparison: >,<,<=,>=
%
%
%Brad Parry, 2015 October

if mod(size(se,1),2) ~= 1 || mod(size(se,2),2) ~= 1
    error('All dimensions of se must be of odd length')
end
if size(se,1) ~= size(se,2)
    error('Dimensions of se must be the same length')
end

%center needs to be chosen so that size of se = 2*center+1
center = floor(size(se,1)/2);
[x,y] = meshgrid(-center:center,-center:center);
%convert the structure element to an array of relative indices
x = x.*se;
y = y.*se;
%remove indices which are NaN
x(isnan(x)) = [];
y(isnan(y)) = [];

%some housekeeping to make sure that we do not try to examine pixels that
%don't exist, that is, pixels that would be outside of the image bounds
ix0 = center;
ix1 = size(im,1)-center+1;
ix2 = size(im,2)-center+1;
%finally extract the 'central' part of the image -- the structuring
%element, se, can be placed anywhere on this image without asking for
%values that do not exist in the full image
IMcenter = im(ix0:ix1,ix0:ix2);

%initialze the array that will store the peaks
peaks = zeros(size(im));
peaks(ix0:ix1,ix0:ix2) = 1;


for k = 1:length(x)
    
    %The central part of the image (IMcenter) will
    %be held still and the full image will be shifted past according to the
    %shifts that were determined from the structuring element se. Each
    %shift is compared against IMcenter by the method specified by
    %operation. Successive shifted comparisons are multiplied together so
    %that all shifts of a given pixel must satisfy the condition determined
    %by operator for that main pixel to be kept as 1 in peaks. For a given
    %pixel r,c to satisfy the operation of all pixels in the region of r,c
    %(as determined by se) must be satisfied.
    
    %create relative indices
    R0 = ix0 + y(k);
    R1 = ix1 + y(k);
    C0 = ix0 + x(k);
    C1 = ix2 + x(k);
    
    if operation == '>'
        peaks(ix0:ix1,ix0:ix2) = peaks(ix0:ix1,ix0:ix2).* (IMcenter > im(R0:R1,C0:C1));
    elseif operation == '>='
        peaks(ix0:ix1,ix0:ix1) = peaks(ix0:ix1,ix0:ix1).* (IMcenter >= im(R0:R1,C0:C1));
    elseif operation == '<'
        peaks(ix0:ix1,ix0:ix1) = peaks(ix0:ix1,ix0:ix1).* (IMcenter < im(R0:R1,C0:C1));
    elseif operation == '<='
        peaks(ix0:ix1,ix0:ix1) = peaks(ix0:ix1,ix0:ix1).* (IMcenter <= im(R0:R1,C0:C1));
    end
    
end
end

function se = ballElement(sz,radius)
%a structuring element of size sz will be constructed so that all pixels within 
%radius of the center will be set to 1 and all pixels of distance > radius from 
%center will be 0
%
%Brad Parry 2015 October
sz = floor(sz/2)*2+1;
center = ceil(sz/2);
[x,y] = meshgrid(1:sz,1:sz);
se = ((x-center).^2 + (y-center).^2).^(1/2) <= radius;
end

function [xy] = g2d(image)
%calculate the centroid, weighted mean of a 2d matrix
image = double(image);
M0 = sum(image(:));
x = repmat(1:size(image,2),[size(image,1), 1]);
Mx = sum(sum(x.*image));
ay(:,1) = 1:size(image,1);
y = repmat(ay,[1,size(image,2)]);
My = sum(sum(y.*image));
xy = [Mx/M0, My/M0];
end