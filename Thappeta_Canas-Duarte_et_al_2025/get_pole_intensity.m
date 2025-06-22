%% Get the fluorescence intensity at the cell pole
%
% - Input
%   1. cellList: Oufti cellList
%   2. im: a specific frame of image
%   3. frame: the corresponding frame ID of the image provided
%   4. bound: a 1 by 4 matrix: [x_min, xmax, y_min, y_max] that is used to
%       define where the poles of the cell are. Also can be a single
%       integer that represents the number of pixels to take at the pole
%       (i.e. top 'bound' closest pixels to the cell edge at the pole)
%   5. varargin
%       - 'nucleoid_field':
%           keyword followed by a char array that represents the name of
%           the field in the cellStruct that stores the information about
%           the objects detected in the cell.
%           Default value: 'objects'
%       - 'skip_nucleoidless'
%           keyword followed by a logical (true or false) flag that decides
%           if to skip cells without any nucleoid detected. If false,
%           values related to the nucleoid will be set to -1.
%           Default value: true
%       - 'output_table'
%           keyword followed by a logical (true or false) flag that
%           indicates if the user wants the final output in the format of a
%           MATLAB table with named column names. If not the raw numeric
%           array (N by 16) will be output.
%           Default value: true
%
% - Output
%   A N by 16 numeric array or table, where N is the number of cells
%   Going from left to right, each column represents:
%   (1) frame ID
%   (2) cell ID
%   (3) cell area [px^2]
%   (4) nucleoid area [px^2]
%   (5) average intensity in the pole further away from the nucleoid
%   (6) average intensity in the pole closer to the nucleoid
%   (7) distance (px) between the further-away pole to the nucleoid
%   (8) distance (px) between the closer pole to the nucleoid
%   (9) average intensity in the entire cell
%   (10) relative nucleoid centroid position (x)
%   (11) relative nucleoid centroid position (y)
%   (12) degree of cell constriction
%   (13) position of cell constriction
%   (14) cell length [px]
%   (15) cell width [px]
%   (16) cell volume [px^3]
%
% - Dependency
%   getextradata.m, projectToMesh.m, constDegreeSingleCell.m
%
% - Author
%   Yingjie Xiang, CJW Lab, Stanford University

function res = get_pole_intensity(cellList,im,frame,bound,varargin)

% Default name of the field storing object information
nucleoid_field = 'objects';
% Default to skip cells without any nucleoid detected
skip_nucleoidless = true;
% Default to output table
output_table = true;

% Parser for optional input arguments
for ii = 1:2:length(varargin)
    switch lower(varargin{ii})
        case 'nucleoid_field'
            nucleoid_field = varargin{ii+1};
        case 'skip_nucleoidless'
            skip_nucleoidless = varargin{ii+1};
        case 'output_table'
            output_table = varargin{ii+1};
        otherwise
            error('Invalid optional arguments');
    end
end

% Lambda expression to calculate the cell width
% We calculate the lengths of individual 'ribs' of the cell mesh
% We then sort the lengths from large to small and average the top 1/3 of
% them. The average is taken as the cell width.
get_cell_width = @(cc) mean(maxk(sqrt((cc.mesh(:,1)-cc.mesh(:,3)).^2+...
    (cc.mesh(:,2)-cc.mesh(:,4)).^2),floor(size(cc.mesh,1)/3)));

% Result container
res = zeros(length(cellList.meshData{frame}),11);

% First colum is the frame number
res(:,1) = frame;

% Go through each cell in this frame
for cc = 1:length(cellList.meshData{frame})
    
    % Load in a cell mesh, which must be a N by 4 matrix
    mesh = cellList.meshData{frame}{cc}.mesh;
    % Validate the cell mesh
    if size(mesh,2) ~= 4
        warning('Frame %d, Cell %d has an invalid cell mesh, skipping now.',frame,cc);
        continue;
    end
    
    % Boolean flag to indicate if this cell contains any detected objects
    has_nucleoid = true;
    % If the 'nucleoid_field' does not exist OR there was no 'outlines'
    % there was no object (i.e. nucleoid)
    if ~isfield(cellList.meshData{frame}{cc},nucleoid_field) ||...
            isempty(cellList.meshData{frame}{cc}.(nucleoid_field).outlines)
        has_nucleoid = false;
    end
    
    % WARNING: if no there was no nucleoid detected and 'skip_nucleoidless'
    % was set to be true, the cell will be skipped SILENTLY.
    if ~has_nucleoid && skip_nucleoidless
        continue;
    end
    
    % Load the box of the cell
    box = cellList.meshData{frame}{cc}.box;
    c0 = box(1);
    r0 = box(2);
    dc = box(3);
    dr = box(4);
    
    % Here we offset the coordinates of the cell mesh,
    % so that we do not need to mask the entire image, but just image
    % within the box of the cell
    mesh(:,[1,3]) = mesh(:,[1,3]) - c0 + 1;
    mesh(:,[2,4]) = mesh(:,[2,4]) - r0 + 1;
    
    % For efficiency, we only deal with the image in the box of the cell
    sub_im = im(r0:r0+dr,c0:c0+dc);
    
    % Convert the N by 4 matrix to N by 2 matrix for the mask construction
    m = double([mesh(:,1:2);flip(mesh(:,3:4),1)]);
    
    % Construct a mask of the cell
    cell_mask = poly2mask(m(:,1),m(:,2),size(sub_im,1),size(sub_im,2));
    
    % Find pixel coordinates within the mask
    [ys,xs] = find(cell_mask);
    
    % Project to cell coordinates
    [L,W] = projectToMesh(xs,ys,mesh);
    
    % Check if need to get the cell length
    if ~isfield(cellList.meshData{frame}{cc},'length')
        cellList.meshData{frame}{cc} = getextradata(cellList.meshData{frame}{cc});
    end
    
    % Normalize by the cell length
    len = cellList.meshData{frame}{cc}.length;
    L = L ./ len - 0.5;
    
    % Normalize by the cell width
    wid = get_cell_width(cellList.meshData{frame}{cc});
    W = W ./ wid;
    
    % Below are nucleoid-related calculations
    nuc_xs_rel = -1;
    nuc_ys_rel = -1;
    nucleoid_area_sum = 0;
    
    if has_nucleoid
        % Get the nucleoid outline
        nuc = cellList.meshData{frame}{cc}.(nucleoid_field).outlines;
        % Nucleoid centroids
        nuc_xs = zeros(length(nuc),1);
        nuc_ys = zeros(length(nuc),1);
        % Go through each nucleoid object
        for ii = 1:length(nuc)
            % Get the centroid of the nucleoid
            nuc_xs(ii) = mean(nuc{ii}(:,1) - c0 + 1);
            nuc_ys(ii) = mean(nuc{ii}(:,2) - r0 + 1);
            % Accumulate the nucleoid area
            nucleoid_area_sum = nucleoid_area_sum + ...
                cellList.meshData{frame}{cc}.(nucleoid_field).area{ii};
        end
        % Get the relative position of the nucleoid centroids
        [nuc_xs_rel, nuc_ys_rel] = projectToMesh(nuc_xs,nuc_ys,mesh);
        nuc_xs_rel = nuc_xs_rel ./ len - 0.5;
        nuc_ys_rel = nuc_ys_rel ./ wid;
    end
    
    % Pixels at the poles
    switch numel(bound)
        case 1
            [sorted_L,ix] = sort(L,'ascend');
            ix(sorted_L < -0.5 | sorted_L > 0.5) = [];
            assert(bound < numel(ix),'Too many pixels for a pole.');
            ixl = ix(1:bound);
            ixr = ix(end-bound:end);
        case 4
            ix = abs(L) > bound(1) & abs(L) < bound(2) &...
                abs(W) > bound(3) & abs(W) < bound(4);
            % Left pole pixels
            ixl = L < 0 & ix;
            % Right pole pixels
            ixr = L > 0 & ix;
        otherwise
            error('Could not determine the poles using the bound.')
    end
    
    if has_nucleoid
        % Get the minimal distance between the nucleoid and the pole
        xl0 = mean(xs(ixl));
        yl0 = mean(ys(ixl));
        dl = min(sqrt((xl0-nuc_xs).^2+(yl0-nuc_ys).^2));
        xr0 = mean(xs(ixr));
        yr0 = mean(ys(ixr));
        dr = min(sqrt((xr0-nuc_xs).^2+(yr0-nuc_ys).^2));
    else
        % Flag to indicate that no pole-nucleoid distance was calculated
        dl = -1;
        dr = -1;
    end
    
    % Pixels only in within the cell outline
    cell_im = cell_mask .* sub_im;
    % Pixels of the entire cell
    pixels_entire_cell = cell_im(:);
    % Pixels at the left pole
    pixels_left_pole = cell_im(sub2ind(size(sub_im),ys(ixl),xs(ixl)));
    % Pixels at the right pole
    pixels_right_pole = cell_im(sub2ind(size(sub_im),ys(ixr),xs(ixr)));
    
    % Remove 0 pixels
    pixels_entire_cell(pixels_entire_cell == 0) = [];
    pixels_left_pole = pixels_left_pole(pixels_left_pole ~= 0);
    pixels_right_pole = pixels_right_pole(pixels_right_pole ~= 0);
    
    % Calculation of constriction
    [constriction_degree,constriction_pos_rel,~] = constDegreeSingleCell(cellList.meshData{frame}{cc},0.0642);
    
    % Collect the results
    res(cc,2) = cc;
    res(cc,3) = cellList.meshData{frame}{cc}.area;
    res(cc,4) = nucleoid_area_sum;
    res(cc,5) = mean(pixels_left_pole);
    res(cc,6) = mean(pixels_right_pole);
    res(cc,7) = dl;
    res(cc,8) = dr;
    res(cc,9) = mean(pixels_entire_cell);
    res(cc,10) = mean(nuc_xs_rel);
    res(cc,11) = mean(nuc_ys_rel);
    res(cc,12) = constriction_degree;
    res(cc,13) = constriction_pos_rel;
    res(cc,14) = len;
    res(cc,15) = wid;
    res(cc,16) = cellList.meshData{frame}{cc}.volume;
end

% Clean up before the output
% Remove any rows with a NaN value
% res(any(isnan(res),2), :) = [];
% Remove any rows with cc == 0, that is, that cell was skipped
res(res(:,2) == 0, :) = [];
% Swap the 7th and 8th column, so that the pole further away from the
% nucleoid is always in the 7th column
if_swap = res(:,7) < res(:,8);
tmp = res(:,7);
res(if_swap,7) = res(if_swap,8);
res(if_swap,8) = tmp(if_swap);
% Similarly, swap the intensity columns
tmp = res(:,5);
res(if_swap,5) = res(if_swap,6);
res(if_swap,6) = tmp(if_swap);

% Check if need to ouput the table
if output_table
    if isempty(res)
        res = table();
    else
        col_names = {'Frame','CellID','CA','NA','iL','iR','dL','dR',...
            'iA','NucX','NucY','Constriction','ConsPos','CLen','CWid','CVol'};
        res = array2table(res,'VariableNames',col_names);
    end
end

end
