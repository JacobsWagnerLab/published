function [SpotDat, def_params] = particleTracking(File1, cellList, parameters, shifts)
%{
-About-

-Inputs-
File1: location of the first fluorescent image in a series of fluorescent
images

cellList: an Oufti created cellList whose first frame corresponds to File1

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

%initialize default parameters.
shift_cellList_meshes = 0;
def_params.pixelsPerSpotRange = [6 50];
def_params.spotlength = 9;
def_params.frequency = .95;
def_params.ithresh = .1;
def_params.nSigmaCutoff = 3;
def_params.fitTol = [-inf inf;-inf inf;-inf inf;-inf inf;-inf inf];

%check input parameters and re-assign them as defaults for current session
if nargin >= 3
    % convert all fields to lower case, this will help interpretation of
    % fieldnames below
    fn = fieldnames(parameters);
    fn_ = cellfun(@lower, fn, 'UniformOutput',0);
    for k = 1:length(fn)
        parameters.(fn_{k}) = parameters.(fn{k});
    end
    %interpret fieldnames and extract values
    if isfield(parameters,'pixelsperspotrange')
        def_params.pixelsPerSpotRange = parameters.pixelsperspotrange;
    end
    if isfield(parameters,'spotlength')
        def_params.spotlength = parameters.spotlength;
    end
    if isfield(parameters,'frequency')
        def_params.frequency = parameters.frequency;
    end
    if isfield(parameters,'ithresh')
        def_params.ithresh = parameters.ithresh;
    end
%     if isfield(parameters,'xytol')
%         def_params.xyTol = parameters.xytol;
%     end
    if isfield(parameters,'nsigmacutoff')
        def_params.nSigmaCutoff = parameters.nsigmacutoff;
    end
    if isfield(parameters,'fittol')
        def_params.fitTol = parameters.fittol;
    end
end

if nargin == 4
    frameRange = 1:size(shifts,1);
    shift_cellList_meshes = 1;
end

%regex stuff in preparation to read across all files
[fluo_file1,ff_ext] = strtok(File1,'.');
ind = regexp(fluo_file1,'\d+$');
first_ind = str2num(fluo_file1(ind:end));

if shift_cellList_meshes == 0
    SpotDat = cell(1,length(cellList.meshData));
    for k = 1:length(cellList.cellId)
        SpotDat{k} = cell(1,max(cellList.cellId{k}));
    end
% was parfor
    parfor Frame = 1:length(cellList.meshData)
        fff=['%0',num2str(length(ind:length(fluo_file1))),'d'];
        filename = [fluo_file1(1:ind-1), num2str(Frame + first_ind - 1,fff), ff_ext];
        image = imread(filename);

        for cell_ix = 1:length(cellList.meshData{Frame})
            %check for a valid cell. if one does not exist, continue
            if length(cellList.meshData{Frame}{cell_ix}.mesh) < 6
                continue
            end
            SpotDatINDEX = cellList.cellId{Frame}(cell_ix);

            %get the cell bounding box and polygonal outline shifted to the
            %box
            box = cellList.meshData{Frame}{cell_ix}.box;
            sub = repmat(box(1:2), [2*size(cellList.meshData{Frame}{cell_ix}.mesh,1), 1]) + 1;
            cellpoly = cat(1,cellList.meshData{Frame}{cell_ix}.mesh(:,1:2),flipud(cellList.meshData{Frame}{cell_ix}.mesh(:,3:4))) - sub;

            img = imcrop(image,box);
            
            %perform spot localization in devTracking
            [TMP] = devTracking(img, cellpoly, def_params);

            if isempty(TMP.Centroid)
                SpotDat{Frame}{SpotDatINDEX}.mesh = cellList.meshData{Frame}{cell_ix}.mesh;
                SpotDat{Frame}{SpotDatINDEX}.box = cellList.meshData{Frame}{cell_ix}.box;
                continue
            end
            
            %assign the data to SpotDat and shift back to full image
            %coordinates
            SpotDat{Frame}{SpotDatINDEX} = TMP;
            SpotDat{Frame}{SpotDatINDEX}.Centroid(:,1) = SpotDat{Frame}{SpotDatINDEX}.Centroid(:,1) + cellList.meshData{Frame}{cell_ix}.box(1) - 1;
            SpotDat{Frame}{SpotDatINDEX}.Centroid(:,2) = SpotDat{Frame}{SpotDatINDEX}.Centroid(:,2) + cellList.meshData{Frame}{cell_ix}.box(2) - 1;
            SpotDat{Frame}{SpotDatINDEX}.mesh = cellList.meshData{Frame}{cell_ix}.mesh;
            SpotDat{Frame}{SpotDatINDEX}.box = cellList.meshData{Frame}{cell_ix}.box;
        end 
        disp(num2str(Frame))
    end
elseif shift_cellList_meshes == 1
    
    %similar to above, used if the images are to be shifted
    
    SpotDat = cell(1,length(frameRange));
    makeNcells = max(cellList.cellId{1});
    for k = 1:length(SpotDat)
        SpotDat{k} = cell(1,makeNcells);
    end
    for k = 1:length(frameRange)
        Frame = frameRange(k);
        fff=['%0',num2str(length(ind:length(fluo_file1))),'d'];
        filename = [fluo_file1(1:ind-1), num2str(Frame + first_ind - 1,fff), ff_ext];
        image = imread(filename);
        
        for cell_ix = 1:length(cellList.meshData{1})
            if length(cellList.meshData{1}{cell_ix}.mesh) < 6
                continue
            end
            
            SpotDatINDEX = cellList.cellId{1}(cell_ix);

            %grab cell mesh data, applying image shifts
            %add box data to spotdat
            box = cellList.meshData{1}{cell_ix}.box;
            box(1:2) = round(box(1:2) - shifts(k,1:2));
            
            mesh = cellList.meshData{1}{cell_ix}.mesh;
            sub = repmat(shifts(k,1:2), [size(mesh,1), 2]);
            cellMesh = mesh - sub;
            
            sub = repmat(box(1:2), [2*size(cellMesh,1), 1]);
            cellpoly = cat(1,cellMesh(:,1:2),flipud(cellMesh(:,3:4))) - sub + 1;

            img = imcrop(image,box);
            [TMP] = devTracking(img, cellpoly, def_params);
            if isempty(TMP.Centroid),continue,end
            SpotDat{k}{SpotDatINDEX} = TMP;
            SpotDat{k}{SpotDatINDEX}.Centroid(:,1) = SpotDat{k}{SpotDatINDEX}.Centroid(:,1) + box(1) - 1;
            SpotDat{k}{SpotDatINDEX}.Centroid(:,2) = SpotDat{k}{SpotDatINDEX}.Centroid(:,2) + box(2) - 1;
            SpotDat{k}{SpotDatINDEX}.box = box;
            SpotDat{k}{SpotDatINDEX}.mesh = cellMesh;
            SpotDat{k}{SpotDatINDEX}.shift = shifts(k,1:2);
        end
        disp(['complete frame ',num2str(Frame)])
    end
end
