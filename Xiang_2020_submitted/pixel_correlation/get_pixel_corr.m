%% Get the correlation between the pixel signals
% -About-
%   This function calculates the correlation between fluoresence signals in
%   two channels (GFP and DAPI).
%
% -Input-
%   - cellList: Oufti cellList
%   - flst: file lists containing paths to image files
%   - type: type of correlation to calculate
%   - nucleoid_mesh_expansion: negative value for shrinking the nucleoid
%
% -Output-
%   - res: a N by 5 numeric array, from left to right each column represents:
%     1. frame
%     2. cell ID
%     3. correlation results calculated for this cell
%     4. correlation score
%     5. area ratio between the shrunk and orignal nucleoid
%
% -Author-
%   Yingjie Xiang, Stanford University, May 2020
function res = get_pixel_corr(cellList,flst,type,nuc_mesh_expansion)

% Turn of polygon warning
warning('off','MATLAB:polyshape:repairedBySimplify');

% Extract file paths
gfp_flst = flst{1};
dapi_flst = flst{2};

res = cell(length(cellList.meshData),1);

% Go through each frame
for frame = 1:length(cellList.meshData)

    % Read in the frame
    gfp = double(imread(fullfile(gfp_flst(frame).folder,gfp_flst(frame).name)));
    dapi = double(imread(fullfile(dapi_flst(frame).folder,dapi_flst(frame).name)));

    % Normalize the image
    gfp = (gfp - min(gfp(:))) ./ (max(gfp(:)) - min(gfp(:)));
    dapi = (dapi - min(dapi(:))) ./ (max(dapi(:)) - min(dapi(:)));

    res_per_frame = ones(length(cellList.meshData{frame}),5).*-2;

    % Go through each cell in this frame
    for cc = 1:length(cellList.meshData{frame})
        % Load in a cell mesh, which must be a N by 4 matrix
        mesh = cellList.meshData{frame}{cc}.mesh;
        % Validate the cell mesh
        if size(mesh,2) ~= 4
            warning('Frame %d, Cell %d has an invalid cell mesh, skipping now.',frame,cc);
            continue;
        end

        % Check if the nucleoid exists
        if ~isfield(cellList.meshData{frame}{cc},'objects') ||...
                isempty(cellList.meshData{frame}{cc}.objects.outlines)
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
        gfp_sub = gfp(r0:r0+dr,c0:c0+dc);
        dapi_sub = dapi(r0:r0+dr,c0:c0+dc);

        % Convert the N by 4 matrix to N by 2 matrix for the mask construction
        m = double([mesh(:,1:2);flip(mesh(:,3:4),1)]);

        % Construct a mask of the cell
        cell_mask = poly2mask(m(:,1),m(:,2),size(gfp_sub,1),size(gfp_sub,2));

        % Get the nucleoid outline
        nuc = cellList.meshData{frame}{cc}.objects.outlines;

        % Construct a mask for the nucleoid
        nuc_mask = false(size(gfp_sub));
        % Go through each nucleoid object
        for ii = 1:length(nuc)
            % Offset the nucleoid position
            nuc{ii}(:,1) = nuc{ii}(:,1) - c0 + 1;
            nuc{ii}(:,2) = nuc{ii}(:,2) - r0 + 1;

            % Construct a polyshape object using the nucleoid outline
            nuc_poly = polyshape(nuc{ii}(:,1),nuc{ii}(:,2));

            % Nucleoid mesh expansion or shrinking
            nuc_poly = polybuffer(nuc_poly,nuc_mesh_expansion);

            % If the nucleoid is highly constricted, shrinking could cause a
            % single object to be broken into multiple (usually 2)
            if nuc_poly.NumRegions > 1
                % Find all the division sites of the object
                division = find(any(isnan(nuc_poly.Vertices),2) == 1);
                division_sites = [division;length(nuc_poly.Vertices)+1];
                assert(nuc_poly.NumRegions == length(division_sites),...
                      'Mismatch between numbers of polygons');
                last_div = 0;
                nuc_chunck_xs = zeros(nuc_poly.NumRegions,1);
                nuc_chunck_ys = zeros(nuc_poly.NumRegions,1);
                % Go through each 'nucleoid chuck'
                for jj = 1:length(division_sites)
                    div = division_sites(jj);
                    chunck = nuc_poly.Vertices(last_div+1:div-1,:);
                    nuc_chunck = polyshape(chunck(:,1),chunck(:,2));
                    [nuc_chunck_xs(jj),nuc_chunck_ys(jj)] = nuc_chunck.centroid;
                    nuc_chunck_mask = poly2mask(chunck(:,1),chunck(:,2),...
                                                size(gfp_sub,1),size(gfp_sub,2));
                    nuc_mask = nuc_mask | nuc_chunck_mask;
                    last_div = div;
                end
            else
                nuc_mask = nuc_mask | poly2mask(nuc_poly.Vertices(:,1),...
                       nuc_poly.Vertices(:,2),size(gfp_sub,1),size(gfp_sub,2));
            end
        end

        % Get only pixels within the nucleoid
        gfp_sub = gfp_sub .* (cell_mask & nuc_mask);
        dapi_sub = dapi_sub .* (cell_mask & nuc_mask);

        g = gfp_sub(:);
        d = dapi_sub(:);
        
        g(g == 0) = [];
        d(d == 0) = [];

        % Collect results
        res_per_frame(cc,1) = frame;
        res_per_frame(cc,2) = cc;
        res_per_frame(cc,3) = length(nuc);
        res_per_frame(cc,4) = corr(g,d,'type',type);
        res_per_frame(cc,5) = sum(nuc_mask(:)) / sum([cellList.meshData{frame}{cc}.objects.area{:}]);
    end
    res_per_frame(res_per_frame(:,3) == -2, :) = [];
    res{frame} = res_per_frame;
end

res = vertcat(res{:});
end
