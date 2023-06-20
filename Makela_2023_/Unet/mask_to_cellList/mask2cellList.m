%% mask2cellList
%
% -Purpose-
%   Generate .mat files containing cellList structure and parameters to be
%   loaded in Oufti. This implementation enables the isolation between the
%   segmentation step and active contour fitting step in Oufti. Better
%   segmentations, such as convolutional neural networks established in
%   Tensorflow could be integrated to Oufti.
%   
%   A general overview of work flow:
%   - read_params() & generate_params() create 'p' parameter struct
%   - generate_cellList calls mask2mesh() to find all meshes given a mask
%   - mask2mesh() find the model from the mask and converts each mask to
%     mesh using model2mesh(),
%   - intxy2,intxy2C,intxyMulti,intxyMultiC,spsmooth are required for
%     model2mesh()
%   - mask2cellList (this function) uses all of above to generate .mat file
%     that can be loaded in Oufti
%
% -Input-
%   - masks (cell array): image masks indicating the cells for each frame
%   - param_filename (char array): path to the .set Oufti parameter file
%   - output_filename (char array): filename for the final output
%   - varargin (bool): boolean flag to indicate if for time lapse
%
% -Output-
%   WARNING: this function writes to the current working directory with
%   the filename specified as output_filename. The resulted .mat file can
%   be loaded in Oufti for further proecessing, e.g. refinement
%
% -Dependency-:
%   1.  intxy2.m
%   2.  intxy2C.mex (compiled by intxy2C.c)
%   3.  intxyMulti.m
%   4.  intxyMultiC.mex (compiled by intxyMultiC.c)
%   5.  spsmooth.m
%   6.  model2mesh.m
%   7.  read_params.m
%   8.  generate_params.m
%   9.  mask2mesh.m
%   10. generate_cellList.m
%
% -Author-
%   Yingjie Xiang, 2019-01-17
%
% -Patch Notes-
%   2019-01-17: created the function
%   2020-11-24: added support for time lapse

function mask2cellList(masks,param_filename,output_filename,is_time_lapse,cell_info)

% Default not time lapse
% is_time_lapse = false;
% if nargin == 4
%     is_time_lapse = varargin{1};
% end

% Create the 2 global variables that Oufti generates natively
% Not declaring them as global may have undefined behavior
% Need a better understanding on the Oufti algorithm to declare these as
% local variables
global p cellList;

% Parse the parameters and generate the structure
p = generate_params(param_filename);
% Some other variables that Oufti natively generates
coefPCA = [];
mCell = [];
paramString = read_params(param_filename);
rawPhaseFolder = [];
shiftfluo = [0,0;0,0];
shiftframes = [];
weights = [];

% Number of images 
im_num = length(masks);

% Memory allocation for arrays contained in cellList
meshData = cell(1,im_num);
cellId = cell(1,im_num);

% Number of cells within each frame 
cellListN = zeros(1,im_num);

% Choice of algorithm
% Create a copy here for the following parfor loop
% algo == 1: pixel based contour fitting
% algo == 4: subpixel based contour fitting
algo = 4;%p.algorithm;

% Loop through each frame
% parfor frame = 1:im_num
for frame = 1:im_num
    % if time lapse
    if is_time_lapse == 1
        % Convert mask to cellList structs for each mask within this frame
        [meshes, unique_bits] = generate_cellList(masks{frame},frame,algo,is_time_lapse,cell_info);
        % Generate the corresponding cellID
        cellId{frame} = unique_bits; 
    else % independent frames
        % Convert mask to cellList structs for each mask within this frame
        meshes = generate_cellList(masks{frame},frame,algo,is_time_lapse);
        % Generate the corresponding cellID
        cellId{frame} = 1:length(meshes);
    end
    % Store the meshes and number of meshes created
    meshData{frame} = meshes;
	cellListN(frame) = length(meshes);
end

%%%%generate_cellList(mask,frame,algo,is_time_lapse,cell_info)

% Finishing touch of the cellList structure
cellList.meshData = meshData;
cellList.cellId = cellId;

% Save varialbes to be loaded by Oufti 
save(output_filename,'cellList','cellListN','coefPCA','mCell','p',...
    'paramString','rawPhaseFolder','shiftfluo','shiftframes','weights');

% Notify the user of the result
fprintf('Saved Oufti cellList to %s.mat \n',output_filename);
end