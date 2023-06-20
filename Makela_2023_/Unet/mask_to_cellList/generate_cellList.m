%% generate_cellList
%
% -Purpose-
%   Generate cellList struct for every cell within the same frame
%
% -Input-
%   - mask (logical array): binary mask indicating the cells for this frame
%   - frame (numeric, int): the frame number 
%   - algo (numeric, int):  choice for contour fitting algorithm
%                           1 for pixel based, 4 for subpixel based
%   - is_time_lapse (bool): boolean flag indicating if time lapse
%
% -Output-
%   - res (cell array): all cellList structs for each cell within the frame
%
% -Varargin-
%   - STEP_SIZE (numeric, double): default to 1.2, unclear, used in model2mesh()
%   - TOLERANCE (numeric, double): default to 1e-5, unclear, used in model2mesh()
%   - MESH_WIDTH (numeric, double): default to 500, unclear, used in model2mesh()
%
% -Dependency-:
%   1.  intxy2.m
%   2.  intxy2C.mex (compiled by intxy2C.c)
%   3.  intxyMulti.m
%   4.  intxyMultiC.mex (compiled by intxyMultiC.c)
%   5.  spsmooth.m
%   6.  model2mesh.m
%   7.  mask2mesh.m
%
% -Author-
%   Yingjie Xiang, 2019-01-17
%
% -Patch Notes-
%   2019-01-17: created the function


function [res, unique_bits] = generate_cellList(mask,frame,algo,is_time_lapse,cell_info)

% Turn off the polyshape warning
warning('off','MATLAB:polyshape:repairedBySimplify');

% Initiate default parameters
switch nargin
    case 2
        algo = 1;
        STEP_SIZE = 1.2;
        TOLERANCE = 1e-5;
        MESH_WIDTH = 500;
    case 3
        STEP_SIZE = 1.2;
        TOLERANCE = 1e-5;
        MESH_WIDTH = 500;
    case 4
        STEP_SIZE = 1.2;
        TOLERANCE = 1e-5;
        MESH_WIDTH = 500;
    case 5
        STEP_SIZE = 1.2;
        TOLERANCE = 1e-5;
        MESH_WIDTH = 500;
    case 6
        STEP_SIZE = 1.2;
        TOLERANCE = 1e-5;
        MESH_WIDTH = 500;
    otherwise
        error('Incorrect number of arguments.');
end

% Core function, convert the mask to meshes and models (cell outlines)
% The variables, STEP_SIZE, TOLERANCE, MESH_WIDTH are parameters for the
% function model2mesh(). The exact algorithm of that function is unclear.
% Care should be taken when selecting these parameters. Failure to
% correctly select the parameters cause abnormal mesh formation or a plain
% output '0' in some cases. The function mask2mesh(), which is a wrapper
% for model2mesh() attempts to resolve these issues.

% whole frame is processed at once (no time lapse information here)
[meshes,models,unique_bits] = mask2mesh(mask,STEP_SIZE,TOLERANCE,MESH_WIDTH);

warning('off')
% Convert all models into Polyshape objects
polygons = cellfun(@(model) polyshape(model), models,'UniformOutput',false);
% Get the bounding boxes for each Polyshape object
[xlims,ylims] = cellfun(@(polygon) polygon.boundingbox,polygons,'UniformOutput',false);
% Convert the xlims and ylims
% Incorrect conversion was seen to cause incorrect display size in Oufti
boxes = cellfun(@(xlim,ylim) [xlim(1),ylim(1),xlim(2)-xlim(1),ylim(2)-ylim(1)],xlims,ylims,'UniformOutput',false);
warning('on')

% Memory allocated to store each cellList struct
res = cell(1,length(meshes));
% Loop through each meshes
for ii = 1:length(meshes)
    % Filling all the fields based on the default Oufti output
    oneCell = struct;
    % Pass in the fitting algorithm choice here
    oneCell.algorithm = algo;
    % Pass in the models and meshes here
    oneCell.model = models{ii};
    oneCell.mesh = meshes{ii};
    oneCell.polarity = 0;
    oneCell.stage = 1;
    oneCell.divisions = [];
    
    if is_time_lapse == 1
        % find correct cell from cell_info
        index = find(unique_bits(ii) == [cell_info.ID]);
        % birth frame
        oneCell.birthframe = cell_info(index).birth;
        % only direct mother cell
        if ~isnan(cell_info(index).parent)
            oneCell.ancestors = cell_info(index).parent;
        else
            oneCell.ancestors = [];
        end
        % only daughter cells
        if ~isempty(cell_info(index).daughters) && length(cell_info(index).daughters) == 2
            oneCell.descendants = cell_info(index).daughters;
        else
            oneCell.descendants = [];
        end
        oneCell.timelapse = 1;
        
    else % individual frames
        oneCell.birthframe = frame;
        oneCell.ancestors = [];
        oneCell.descendants = [];
        oneCell.timelapse = 0;
    end
    
    % Pass in the bounding box for the cell here
    oneCell.box = boxes{ii};
    % Store the result from a single cell
    res{ii} = oneCell;
end

end

