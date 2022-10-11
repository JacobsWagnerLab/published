%% Add spots detected in fluorescent images to existing cellList with meshes
% -Input-
%   - cellList: Oufti cellList
%   - intensityRatioThreshold: threshold for Modified_Find_Irregular_Spots
%   - varargin:
%       Provide either 'che' or 'gfp' followed by path to images, case does
%       not matter.
%       - 'che': keyword, followed by path to mCherry images
%       - 'gfp': keyword, followed by path to GFP images
%
% -Output-
%   - cellList with additional fields within the channel (either 'che' or 'gfp') field:
%       1. spotPosition
%           Output from Modified_Find_Irregular_Spots
%           Spot position in image coordinates
%       2. spotIntensity
%           Output from Modified_Find_Irregular_Spots
%           Sum of pixel intensity arond the intensity peak of the spot
%           The spot has a default "radius" of 5 px, see Modified_Find_Irregular_Spots
%       3. intensityRatio
%           Output from Modified_Find_Irregular_Spots
%           Ratio between spotIntensity to a shell surrounding the spot, see Modified_Find_Irregular_Spots
%       4. spotPosition_rel
%           Spot positions relative to the cell mesh
%       5. pixels_in_cell
%           Pixel values within the cell mesh in the fluoresence channel
%       6. pixels_outside_spots
%           Sames Item 5 above, except all the pixels that are within 5
%           pixels to the spot centers (based on spotPosition) are removed.
%       7. cell_length_in_um
%           Cell length in micrometer (assumed px = PIXEL_SIZE)
%       8. average_fluor
%           Average fluoresence intensity within the cell mesh
%       9. number_of_spots
%           Number of spot detected
%       10. spot_density_per_10um
%           Number of spots per 10 micrometers
%       11. mean_spot_intensity
%           Average fluoresence among all spots detected
%       12. cv_spot_intensity
%           Standard deviation of fluoresence among all spots divided by
%           item 11
%       13. spot_fluor_density
%           Spots fluoresence divided by 121, which is the default area of
%           each spot (in pixels), assuming the "radius" is 5
%       14. mean_spot_fluor_density
%           Average value of Item 13
%       15. mean_nonspot_fluor_density
%           Average intensity of pixels outside the spots in a cell
%       16. mean_spot_nonspot_fluor_density_diff
%           Difference between Items 14 and 15
%       17. interpuncta_distance
%           The smallest distance between two nearby puncta in the same
%           channel. The unit is micrometer
%       18. edge_distance
%           The minimal distances between the ends of the cell to the
%           nearest puncta detected. The unit is micrometer
%
%       Also added fields added by getextradata.
%
% -Author-
%   Yingjie Xiang, 2020-04-02, Jacobs-Wagner Lab, Stanford Universiy


function cellList_with_spots = add_spots_to_cellList(cellList,intensityRatioThreshold,varargin)

% Global variable here, this is the pixel number of a spot
% By default, fitRadius = 5 in Modified_Find_Irregular_Spots.m
% With the center pixel at the peak, the total number of pixels for a peak
% is 11x11 = 121
SPOT_AREA = 121;

% Default pixel size (micron / pixel)
PIXEL_SIZE = 0.0643;

% Boolean flag to indicate if the path to GFP images and/or mCherry images
% were provided
has_che = false;
has_gfp = false;
channel = '';

% Parser for optional input arguments
for ii = 1:2:length(varargin)
    switch lower(varargin{ii})
        case 'che'
            path_to_images = varargin{ii+1};
            has_che = true;
            channel = 'che';
        case 'gfp'
            path_to_images = varargin{ii+1};
            has_gfp = true;
            channel = 'gfp';
        otherwise
            error(['Invalid optional arguments. Options: ''che''',...
                ' or ''gfp'' followed by the path to images']);
    end
end

% Check if either GFP or mCherry images were provided
assert(has_che | has_gfp, 'Neither mCherry nor GFP images were provided.');
assert(xor(has_che,has_gfp), 'Provide either mCherry or GFP images, but not both.');

% Add image path to the cellList if there is any che or gfp image
if has_che
    cellList.che_path = path_to_images;
elseif has_gfp
    cellList.gfp_path = path_to_images;
end

% Add spots to the cellList
% Load the corresponding images
images = loadimseries(path_to_images);
% Copy the cellList to avoid changing the original cellList
cellList_with_spots = cellList;
% Add spots data
cellList_with_spots = add_spots(cellList_with_spots,images,intensityRatioThreshold,channel);
% Add pixels within the cell mesh in the fluoresence channel into cellList
cellList_with_spots = add_pixels(cellList_with_spots,images,channel);
% Add spot statistics
cellList_with_spots = add_stats(cellList_with_spots,channel);
% Add distance calculations
cellList_with_spots = add_distances(cellList_with_spots,channel);

% Modified field name and add some fields
    function cl = tidy_cellList(cl,channel)
        assert(strcmp(channel,'che') || strcmp(channel,'gfp'),'False channel name.');
        for frame = 1:length(cl.meshData)
            for cc = 1:length(cl.meshData{frame})
                % Check if spots were detected for this cell
                if isfield(cl.meshData{frame}{cc},'spotPosition')
                    % Move the spotPosition into the channel field
                    cl.meshData{frame}{cc}.(channel).spotPosition = cl.meshData{frame}{cc}.spotPosition;
                    % Calculate the relative coordinates
                    [lA,dA] = projectToMesh(cl.meshData{frame}{cc}.spotPosition(:,1),...
                        cl.meshData{frame}{cc}.spotPosition(:,2),cl.meshData{frame}{cc}.mesh);
                    cl.meshData{frame}{cc}.(channel).spotPosition_rel = [lA,dA];
                    % Remove the spotPosition field
                    cl.meshData{frame}{cc} = rmfield(cl.meshData{frame}{cc},'spotPosition');
                end
                if isfield(cl.meshData{frame}{cc},'spotIntensity')
                    % Move the spotIntensity into the channel field
                    cl.meshData{frame}{cc}.(channel).spotIntensity = cl.meshData{frame}{cc}.spotIntensity;
                    % Remove the spotIntensity field
                    cl.meshData{frame}{cc} = rmfield(cl.meshData{frame}{cc},'spotIntensity');
                end
                if isfield(cl.meshData{frame}{cc},'intensityRatio')
                    % Move the intensityRatio into the channel field
                    cl.meshData{frame}{cc}.(channel).intensityRatio = cl.meshData{frame}{cc}.intensityRatio;
                    % Remove the intensityRatio field
                    cl.meshData{frame}{cc} = rmfield(cl.meshData{frame}{cc},'intensityRatio');
                end
            end
        end
    end

% Add spots detected by the "Modified_Find_Irregular_Spots" to the cellList
    function cl = add_spots(cl,images,ithreshold,channel)
        fprintf('Adding %s spots to the cellList. \n',channel);
        cl = Modified_Find_Irregular_Spots(cl, images, length(cl.meshData), ...
            'intensityRatioThreshold',ithreshold);
        cl = tidy_cellList(cl,channel);
    end

% Get pixel intensity within the cell mesh
    function cl = add_pixels(cl,images,channel)
        fprintf('Adding %s pixels to the cellList.\n',channel);
        assert(strcmp(channel,'che') || strcmp(channel,'gfp'),'False channel name.');
        for frame = 1:length(cl.meshData)
            % Read a frame
            % im0 and im1 are the same, except that in im1, the "spot
            % pixels" (5-pixel from the spot center) were zeroed out
            % The length difference between im0 and im1 should therefore be
            % around N x 121, where N is the total number of spots detected
            im0 = double(images(:,:,frame));
            im1 = double(images(:,:,frame));
            has_spots = false;
            for cc = 1:length(cl.meshData{frame})
                if isfield(cl.meshData{frame}{cc},channel) && isfield(cl.meshData{frame}{cc}.(channel),'spotPosition')
                    has_spots = true;
                    % Zero out pixels around the spot centers
                    spots = round(cl.meshData{frame}{cc}.(channel).spotPosition);
                    for kk = 1:size(spots,1)
                        % 5 pixels were zeroed out in each direction
                        im1(spots(kk,2)-5:spots(kk,2)+5,spots(kk,1)-5:spots(kk,1)+5) = 0;
                    end
                end
                mesh = cl.meshData{frame}{cc}.mesh;
                % Check if the mesh is valid
                if size(mesh,2) == 4
                    mesh = double(flip([mesh(:,1:2);flip(mesh(:,3:4),1)],1));
                    % Create the cell mask
                    cell_mask = poly2mask(mesh(:,1),mesh(:,2),size(im0,1),size(im0,2));
                    % Find pixels within the cell mesh
                    cell_im0 = im0.*cell_mask;
                    cell_im0 = cell_im0(:);
                    cell_im0 = cell_im0(cell_im0 ~= 0);
                    % Store the pixels in the channel field
                    cl.meshData{frame}{cc}.(channel).pixels_in_cell = cell_im0;
                    % Do the same for the image with "spot pixels" removed
                    if has_spots
                        cell_im1 = im1.*cell_mask;
                        cell_im1 = cell_im1(:);
                        cell_im1 = cell_im1(cell_im1 ~= 0);
                        cl.meshData{frame}{cc}.(channel).pixels_outside_spot = cell_im1;
                    end
                end
            end
        end
    end

% Add statistics calculated
    function cl = add_stats(cl,channel)
        fprintf('Calculating statistics for the %s spots.\n',channel);
        assert(strcmp(channel,'che') || strcmp(channel,'gfp'),'False channel name.');
        for frame = 1:length(cl.meshData)
            for cc = 1:length(cl.meshData{frame})
                % Check if the mesh is valid
                if size(cl.meshData{frame}{cc}.mesh,2) ~= 4
                    continue;
                end
                % Check if the length field exists
                if ~isfield(cl.meshData{frame}{cc},'length')
                    cl.meshData{frame}{cc} = getextradata(cl.meshData{frame}{cc});
                end
                % Convert cell length to micrometers
                cl.meshData{frame}{cc}.cell_length_in_um = cl.meshData{frame}{cc}.length*PIXEL_SIZE;
                % Average pixel intensity within the cell mesh
                cl.meshData{frame}{cc}.(channel).average_fluor = mean(cl.meshData{frame}{cc}.(channel).pixels_in_cell);
                % Skip the cell if no spot was detected
                if ~isfield(cl.meshData{frame}{cc}.(channel),'spotPosition')
                    cl.meshData{frame}{cc}.(channel).number_of_spots = 0;
                    cl.meshData{frame}{cc}.(channel).spot_density_per_10um = 0;
                    continue;
                end
                % Number of spots detected
                cl.meshData{frame}{cc}.(channel).number_of_spots = size(cl.meshData{frame}{cc}.(channel).spotPosition,1);
                % Number of spots per 10 micrometers
                cl.meshData{frame}{cc}.(channel).spot_density_per_10um = 10 * cl.meshData{frame}{cc}.(channel).number_of_spots / cl.meshData{frame}{cc}.cell_length_in_um;
                % Average pixel intensity within all spots
                cl.meshData{frame}{cc}.(channel).mean_spot_intensity = mean(cl.meshData{frame}{cc}.(channel).spotIntensity);
                % Std dev of pixel intensity / Average pixel intensity within all spots
                cl.meshData{frame}{cc}.(channel).cv_spot_intensity = std(cl.meshData{frame}{cc}.(channel).spotIntensity) / mean(cl.meshData{frame}{cc}.(channel).spotIntensity);
                % Spot fluoresence density
                cl.meshData{frame}{cc}.(channel).spot_fluor_density = cl.meshData{frame}{cc}.(channel).spotIntensity ./ SPOT_AREA;
                % Average spot fluoresence density
                cl.meshData{frame}{cc}.(channel).mean_spot_fluor_density = mean(cl.meshData{frame}{cc}.(channel).spot_fluor_density);
                % Average fluoresence density of non-spot pixels
                cl.meshData{frame}{cc}.(channel).mean_nonspot_fluor_density = mean(cl.meshData{frame}{cc}.(channel).pixels_outside_spot);
                % Difference between the average fluoresence density between spot and non-spot pixels
                cl.meshData{frame}{cc}.(channel).mean_spot_nonspot_fluor_density_diff = cl.meshData{frame}{cc}.(channel).mean_spot_fluor_density - cl.meshData{frame}{cc}.(channel).mean_nonspot_fluor_density;
            end
        end
    end

% Add distance related calculations
    function cl = add_distances(cl,channel)
        fprintf('Calculating distances for the %s spots.\n',channel);
        assert(strcmp(channel,'che') || strcmp(channel,'gfp'),'False channel name.');
        for frame = 1:length(cl.meshData)
            for cc = 1:length(cl.meshData{frame})
                if isfield(cl.meshData{frame}{cc},channel) && isfield(cl.meshData{frame}{cc}.(channel),'spotPosition_rel')
                    % Get the interpuncta distances
                    cl.meshData{frame}{cc}.(channel).interpuncta_distance = ...
                        helper(cl.meshData{frame}{cc}.(channel).spotPosition_rel).*PIXEL_SIZE;
                    % Get the edge distances
                    max_pos = max(cl.meshData{frame}{cc}.(channel).spotPosition_rel(:,1));
                    min_pos = min(cl.meshData{frame}{cc}.(channel).spotPosition_rel(:,1));
                    cl.meshData{frame}{cc}.(channel).edge_distance = ...
                        [cl.meshData{frame}{cc}.length-max_pos;min_pos].*PIXEL_SIZE;
                end
            end
        end
    end

% Helper function for calculating the distance
    function res = helper(pos)
        % Sort the position based on the first column
        % Borrelia cells are long and narrow, the first column here is the relative
        % particle locations (cell coordinates), x. The second column is almost
        % always close to 0, because the particles are almost always sandwiched at
        % the center by two inner membranes.
        [~,sort_ix] = sort(pos(:,1));
        pos = pos(sort_ix,:);
        % Calculate the euclidean distances between two neighboring points
        res = sqrt(sum((pos(2:end,:) - pos(1:end-1,:)).^2,2));
    end

end