function [SpotDat] = devTracking(image, ROI, parameters)
%{
-About-
Performs 2D localization of diffraction limited spots in an image delimited
by a provided ROI over parameters given. If no ROI is given, particle
localization is performed over the entire image. Localization is performed
by fitting raw pixel intensity values to a 2D Gaussian

-Inputs-
image: the image to perform spot identification and localization in

ROI: region of interest for spot localization, should be an n-by-2 polygon
specifying valid image locations to search for spots

Parameters: a structure of parameters controlling 2d Gaussian fitting
    pixelsPerSpotRange: the lower and upper bounds of the number of pixels
        in a spot. If [6,90], a spot must have at least 6 pixels but fewer
        than 90
    spotlength: anticipated spot size, used for construction of bandpass
        filter for noise removal prior to spot identification
    frequency: frequency response for bandpass filter
    ithresh: modifies threshold value determined from Otsu's method: after
        the image is filtered, Otsu's method is used to naively segment it.
        it ithresh = 1, the value from Otsu's method will be used, if
        ithresh = 0.1, 0.1*otsu will be used
    nSigmaCutoff: a spot must be nSigmaCutoff above background noise to be
        labeled as a spot
    fitTol: upper and lower bounds of spot statistics to be kept as spots
        checks the error on fitted terms [A,k,sigma,x,y]

shifts: an n-by-2 matrix specifying image shifts for each frame in x
    (shifts(:,1)) and y (shifts(:,2)) directions. the nth row corresponds to
    the nth frame. if shifts are not provided, the images will be analzed
    as is

-varargin-
na

-Outputs-
SpotDat: nested structure of cell arrays {Frame}{Cell}.(tracking output)
with fields (the kth element of each field corresponds to the kth spot):
    Centroid: xy coordinates of spot peaks
    rmse: root mean-squared error of the fit
    sigma: sigma of the fitted gaussian function
    gtr: a measure of spot intensity over the background
    confint: confidence intervals of the
    fitdata: statistics of the fit

-Example-
   
-Supplementary-

-Keywords-

-Dependencies-
imcrop.m
devTracking.m

-References-

-Author-
Brad Parry, 2014 July 28
%}

%initialize the output
SpotDat.Centroid = [];
SpotDat.rmse = [];
SpotDat.sigma = [];
SpotDat.gtr = [];
SpotDat.confint = [];
SpotDat.fitdata = [];

%initialize default parameters
def_params.pixelsPerSpotRange = [6 50];
def_params.spotlength = 9;
def_params.frequency = .95;
def_params.ithresh = .1;
def_params.fitTol = [-inf inf;-inf inf;-inf inf;-inf inf;-inf inf];
def_params.nSigmaCutoff = 2;

%sort parameters and populate workspace with their corresponding variables
%this segment extracts fieldnames from the parameters structure and creates
%those variables explicitly in the workspace
fields = fieldnames(def_params);
if ~exist('parameters','var') || ~isstruct(parameters)
    for k = 1:length(fields)
        v = genvarname(fields{k});
        eval([v ' = def_params.(fields{k});']);
    end
else
    a = fieldnames(parameters);
    badfields = setdiff(a, fields);
    parameters = rmfield(parameters, badfields);
    if ~isempty(badfields)
        warning('Bad fieldnames were found in the parameters list')
        tx = ['The following field names were removed from parameters:'];
        disp(tx)
        for k = 1:length(badfields)
            disp(['    ',char(badfields{k})])
        end
    end
    fieldslogical = isfield(parameters,cellfun(@(x) x, fields,'UniformOutput',0));
    needfields = ~fieldslogical;
    ix = needfields.*(1:length(fields))'; ix(ix == 0) = [];
    for k = ix'
        parameters.(fields{k}) = def_params.(fields{k});
    end
    %sprouting variables...
    for k = 1:length(fields)
        v = genvarname(fields{k});
        eval([v ' = parameters.(fields{k});']);
    end
end


[r,c] = size(image);
image = double(image);
dxFromCenter = spotlength / 2 + 1;
if nargin == 1 || isempty(ROI)
    %if a Region of interest in not provided, compute a region of interest
    %equivalent to the entire frame
    y = [1:r,repmat(r,[1,c]),r:-1:1,ones([1,c])]';
    x = [ones([1,r]),1:c,repmat(c,[1,r]),c:-1:1]';
    ROI = cat(2,x,y);
end
ROI = double(ROI);
ROIpxls = poly2mask(ROI(:,1),ROI(:,2),r,c);
filteredImage = filterIM(image,'frequency',frequency,'sl',spotlength,'ithresh',ithresh);
thFiltImage = filteredImage;
thFiltImage(thFiltImage>0) = 1;

imgMasked = image.*double(ROIpxls);

% find the mean intensity in the ROI
meanCellInt = mean(imgMasked(imgMasked>0));
stdCellInt = std(imgMasked(imgMasked>0));

%find rough location of peaks, by guessing they should be contain the
%bright pixels in the image
[pixvals, sortorder] = sort(imgMasked(:),'descend');

%save some computation time later by ignoring relatively dim pixels
%initially
del = pixvals < (meanCellInt + nSigmaCutoff*stdCellInt);
sortorder(del) = [];

%build matrix of neighbors to ignore the bright neighbors of brightest
%pixels, but ignore pixels adjacent to edge pixels
del = 1:r:c*r;
for k = [1:r, del, del + r-1, ((c-1)*r + 1):c*r]
    sortorder(sortorder == k) = [];
end

%check if the initial peak estimates are local maxima, that is check if
%they are brighter than their neighbors.
repeatedvals = repmat(image(sortorder),[1,8]);
neighbors = [image(sortorder - r - 1),image(sortorder - r),image(sortorder - r + 1),image(sortorder + 1),...
    image(sortorder + r + 1),image(sortorder + r),image(sortorder + r - 1),image(sortorder - 1)];
di = repeatedvals - neighbors;
peaks = sortorder(min(di,[],2)>0);
[n, m] = meshgrid(1:c,1:r);

%eliminate the dimmest peak of two peaks that are close together
%index peaks as peaks(:) to ensure a column vector
rowlocation = (rem(peaks(:)-1,r)+1);
collocation = ceil(peaks(:)./r);

dr = repmat(rowlocation, [1, length(rowlocation)]) - repmat(rowlocation', [length(rowlocation), 1]);
dc = repmat(collocation, [1, length(collocation)]) - repmat(collocation', [length(collocation), 1]);
drc = dxFromCenter*tril(ones(size(dr)))+(dr.^2 + dc.^2).^(1/2);
% the rows contain the indices of the brightest of two peaks dxFromCenter apart
[~, cHat] = find(drc < dxFromCenter);
peaks(cHat) = [];

%iterate through the initial peak estimates and try to perform fitting
for ind = peaks(:)'
    rowlocation = (rem(ind-1,r)+1);
    collocation = ceil(ind./r);
    if thFiltImage(rowlocation,collocation) == 0,continue, end

%         choose which pixels to use in fit...
%         -choose pixels within dxFromCenter of the peak and of value '1' in
%         the thresholded, filtered image
    mask = ((m-rowlocation).^2 + (n-collocation).^2).^(1/2);
    mask(mask <= dxFromCenter) = -1;
    mask(mask>0) = 0;        
    mask(mask==-1) = 1;

%         Continue if there is no spot
    if sum(mask(:)) == 0, continue, end
%         extract the box of ones in mask
    cs = (sum(mask) > 0) .* n(1,:);
    cs = [min(cs(cs>0)),max(cs)];
    rs = (sum(mask,2) > 0) .* m(:,1);
    rs = [min(rs(rs>0)),max(rs)];

    subImage = image(rs(1):rs(2),cs(1):cs(2));
%         continue if the spot is too small
    if length(subImage(:)) < pixelsPerSpotRange(1), continue, end
    
    %try to fit the cropped image
    [fitdata, gof] = gaussfit2d(subImage);
    if isempty(fitdata),continue,end
%     construct errr that is the ratio of confidence intervals to the
%     fitted data, in other words a normalized confidence range.
    errr = diff(confint(fitdata)) ./ ([fitdata.A,fitdata.k,fitdata.sigma,fitdata.x0,fitdata.y0]);
    errr = errr(:);
    %verify that the errors are within the specified fit tolerance
    checkErr = (errr >= fitTol(:,1)) .* (errr <= fitTol(:,2));

%     transform to fit tolerances
    if sum(checkErr) == length(checkErr)

        SpotDat.confint{end+1} = confint(fitdata);
        SpotDat.gtr(end+1) = (image(ind)-meanCellInt)/stdCellInt;
        SpotDat.Centroid(end+1,1:2) = [cs(1)+fitdata.x0-1, rs(1)+fitdata.y0-1];
        SpotDat.sigma(end+1) = fitdata.sigma;
        SpotDat.rmse(end + 1) = gof.rmse;
        SpotDat.fitdata{end+1} = [fitdata.A,fitdata.k,fitdata.sigma,fitdata.x0,fitdata.y0];
    end
end