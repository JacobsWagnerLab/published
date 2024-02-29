function trackCells
% Tracks cells based on overlap from masks files where each cell has it's
% separate ID. For drift correction, consecutive frames are aligned using
% phase correlation (only translation, no rotation or scaling). The 
% aligment is only used for tracking and not saved to actual image data. 
% If multiple masks are paired with the same mask, this is counted as 
% division event. Single frame joining/splitting are fixed automatically.
% Re-writes masks IDs with new unique cell ID and saves cell info to a
% cell_info.mat file. Also, estimates when each XY should be truncated and
% saves the last frame to XY_limits vector

% input: a folder with mask images with unique number for each mask
%       and background as zeros. 
%       filename format([name]t**xy**c1_mask.tif)(e.g. name_t01xy01c1_mask.tif)

% output: overwrites the images in the folder with unique cell IDs that are
%       consistent across time series

%       cell_info - a struct with each row corresponding to an unique  cell 
%       with following fields:
%           ID - cell ID corresponding to ID in Oufti (for all cell frames)
%           birth - frame when the cell was born
%           parent - ID of parent cell
%           daughters - IDs of daughter cells 

% 7/16/2021 Jarno Makela

size_min = 150;                 % minimum size for a cell in px  

% for truncation of XY
growth_fold_change = 0.5;       % growth treshold below average for truncation
window_size = 3;                % how many frames below the treshold are required    

% choose folder with masks
pathname = pwd;
% pathname = uigetdir(pwd, 'Select a folder for Masks');
images_mask_all = dir(fullfile(pathname, '*c1_mask.*'));

% find image mask files and extract XY number
for ii = 1:length(images_mask_all)
    newStr = split(images_mask_all(ii).name,["xy","c1_mask"]);
    images_mask_all(ii).XY = str2num(newStr{2});
end

% struct to store cell info on parent, daughters, birth, division
cell_info = struct;
cell_info.XY = [];
cell_info.ID = [];
cell_info.birth = [];
cell_info.division = [];
cell_info.parent = [];
cell_info.daughters = [];
cell_info.area = [];
cell_info.coord = [];
n = 1;

uniq_XY = unique([images_mask_all.XY]);
XY_areas = NaN(length(uniq_XY),length(images_mask_all)./length(uniq_XY));
for hh = 1:length(uniq_XY)  % XY fields
    % image_mask files for single XY
    images_mask = images_mask_all([images_mask_all.XY] == uniq_XY(hh));
    
    % linking matrix (rows frames, columns masks)
    % parent mask number in previous frame
    % no parent - 0, no mask for the number - NaN
    linking_masks = NaN(length(images_mask),10000);
    % reverse direction linking for problem divisions
    linking_masks_rev = NaN(length(images_mask),10000);
    
	% cell areas for future use
    cell_areas = NaN(length(images_mask),10000);
    
    % cell coordinates for future use
    x_coord = NaN(length(images_mask),10000);
    y_coord = NaN(length(images_mask),10000);
    
    disp('Calculating cell mask overlap between frames')
    for kk = length(images_mask):-1:2 % time (from end to beginning)
        % mask1 is later (t) than mask2 (t-1) in time
        bw_mask1 = imread(fullfile(pathname,images_mask(kk).name));
        bw_mask2 = imread(fullfile(pathname,images_mask(kk-1).name));
        
        % total mask areas for growth estimates
        XY_areas(hh,kk) = sum(bw_mask1(:)>0);
        if kk == 2
            XY_areas(hh,1) = sum(bw_mask2(:)>0);
        end
        
        % remove masks touching image border
        border_mask = true(size(bw_mask1));
        border_mask(2:(size(bw_mask1,1)-1),2:(size(bw_mask1,2)-1)) = 0;
        ind1 = ismember(bw_mask1, unique(bw_mask1(border_mask)));
        bw_mask1(ind1) = 0;
        ind2 = ismember(bw_mask2, unique(bw_mask2(border_mask)));
        bw_mask2(ind2) = 0;

        % remove masks smaller than X (4-connected) 
        % by multiplying the mask on the image
        bw_mask1 = bw_mask1.*uint16(bwareaopen(bw_mask1,size_min,4));
        bw_mask2 = bw_mask2.*uint16(bwareaopen(bw_mask2,size_min,4));
        
        % estimate translation to align two masks using phase
        % correlation and apply it to bw_mask2 (no interpolation)
        % Note! this only affects linking, not saved to masks
        tformEstimate = imregcorr(bw_mask2>0,bw_mask1>0,'translation');
        Rfixed = imref2d(size(bw_mask1));
        bw_mask2 = imwarp(bw_mask2,tformEstimate,'OutputView',Rfixed,'interp','nearest');

        % cellIDs
        cellIDs1 = unique(bw_mask1(:));
        cellIDs1 = double(cellIDs1(cellIDs1 ~= 0));
        cellIDs2 = unique(bw_mask2(:));
        cellIDs2 = double(cellIDs2(cellIDs2 ~= 0));

        % for first frame save cell areas
        if kk == 2
            % go through masks
            for jj = 1:length(cellIDs2)  
                mask_cell_first = bw_mask2 == cellIDs2(jj);
                cell_areas(kk-1,cellIDs2(jj)) = sum(mask_cell_first(:));
                % coords
                [y,x] = find(mask_cell_first);
                x_coord(kk-1,cellIDs2(jj)) = mean(x);
                y_coord(kk-1,cellIDs2(jj)) = mean(y);
            end
        end

        % overlap from later time point to earlier one (t -> t-1)
        for ii = 1:length(cellIDs1) % masks
            % cell mask from first time point
            mask_cell = bw_mask1 == cellIDs1(ii);
            
            % cell areas
            cell_areas(kk,cellIDs1(ii)) = sum(mask_cell(:));
            
            % coords
            [y,x] = find(mask_cell);
            x_coord(kk,cellIDs1(ii)) = mean(x);
            y_coord(kk,cellIDs1(ii)) = mean(y);

            % mask the later time point using cell mask
            overlap_mask = mask_cell.*double(bw_mask2);
            overlap_mask = overlap_mask(overlap_mask ~= 0);

            % find the largest overlap by number of mask indeces
            edges = unique(overlap_mask);
            counts = histc(overlap_mask,edges);
            [~, ind] = max(counts);
            maxOverlapId = edges(ind);

            if ~isempty(maxOverlapId)
                % save to linking matrix, position defines mask number and
                % frame, value defines mask number from previous frame
                linking_masks(kk,cellIDs1(ii)) = maxOverlapId;
            else
                % if no parent found, put 0
                linking_masks(kk,cellIDs1(ii)) = 0;
            end
        end

        % opposite direction, overlap from earlier to later time (t-1 -> t)
        for ii = 1:length(cellIDs2) % masks
            % cell mask from first time point
            mask_cell = bw_mask2 == cellIDs2(ii);

            % mask the earlier time point using cell mask
            overlap_mask = mask_cell.*double(bw_mask1);
            overlap_mask = overlap_mask(overlap_mask ~= 0);

            % find the largest overlap by number of mask indeces
            edges = unique(overlap_mask);
            counts = histc(overlap_mask,edges);
            [~, ind] = max(counts);
            maxOverlapId = edges(ind);

            if ~isempty(maxOverlapId)
                % save to linking matrix, position defines mask number and
                % frame, value defines mask number from previous frame
                linking_masks_rev(kk-1,cellIDs2(ii)) = maxOverlapId;
            end
        end
    end
    % first frame values
    linking_masks(1,cellIDs2) = cellIDs2;
    % remove extra columns with only NaNs
    % limit is the maximum number of mask in all frames
    limit = find((sum(~isnan(linking_masks),1) > 0) == 1,1,'last');
    linking_masks = linking_masks(:,1:limit);
    linking_masks_rev = linking_masks_rev(:,1:limit);
    cell_areas = cell_areas(:,1:limit);
    x_coord = x_coord(:,1:limit);
    y_coord = y_coord(:,1:limit);

    % unique node numbers for each mask and frame
    % rows and columns are flipped to get indexing by time
    linking_nodes = linking_masks';
    mask = ~isnan(linking_nodes);
    indeces = find(~isnan(linking_nodes));
    linking_nodes(indeces) = 1:sum(mask(:));
    linking_nodes = linking_nodes';

    % generate directed graph from this data using both masks (front and reverse)
    % edges: s -> t
    disp('Generating a directed network of cell masks')
    s = NaN(size(linking_masks,1),2*size(linking_masks,2));
    t = NaN(size(linking_masks,1),2*size(linking_masks,2));
    for ii = 1:size(linking_masks,2) % masks
        % mask from time -> time-1, linking_masks
        for jj = 2:size(linking_masks,1) % times
            % ignore Nan and 0 that indicate new cells
            if ~isnan(linking_masks(jj,ii)) && linking_masks(jj,ii) ~= 0
                s(jj,ii) = linking_nodes(jj-1,linking_masks(jj,ii));
                t(jj,ii) = linking_nodes(jj,ii);
            end
        end
        % mask from time-1 -> time, linking_masks_rev
        for jj = 1:(size(linking_masks_rev,1)-1) % times
            % ignore Nan
            if ~isnan(linking_masks_rev(jj,ii))
                s(jj,ii + size(linking_masks,2)) = linking_nodes(jj,ii);
                t(jj,ii + size(linking_masks,2)) = linking_nodes(jj+1,linking_masks_rev(jj,ii));
            end
        end
    end
    % concatenate and remove NaNs
    s = s(:);
    t = t(:);
    s = s(~isnan(s));
    t = t(~isnan(t));

    % digraph and remove duplicated edges
    G = digraph(s,t,[],sum(mask(:)));
    G = simplify(G);

    % create a struct to hold information of all nodes:
    nodes = struct;
    nodes.cellID = [];          % cell ID
    nodes.parentNodeID = [];    % node ID of parent
    nodes.daughterNodeID = [];  % node ID of daughter
    nodes.inedges = [];         % inedges
    nodes.outedges = [];        % outedges
    nodes.splittedNodeIDs = []; % nodeID was split between existing node (first) and new nodeID (second)
    nodes.joinedNodeIDs = [];   % nodeIDs were joined by removing edges to nodeID to second (does not have connection anymore)
    nodes.maskCoord = [];       % mask coordinates
    nodes.area = [];            % mask area
    nodes.frame = [];           % frame number
    nodes.maskID = [];          % mask ID

    % pre-cache the correct size
    nodes(height(G.Nodes)).cellID = [];

    % fix 1 frame joining by finding nodes with inedges and outedges == 2 and creating a new
    % node and redirecting edges
    for ii = 1:height(G.Nodes) % loop over nodes
        [~, inid] = inedges(G,ii);
        [~, onid] = outedges(G,ii);
        nodes(ii).inedges = inid;
        nodes(ii).outedges = onid;

        % nodes with 2 of both in and out edges
        if length(inid) == 2 && length(onid) == 2 
            % create new node that gets the extra connections
            G = addnode(G,1);
            newIndex = height(G.Nodes);
            G = addedge(G, [inid(2) newIndex], [newIndex onid(2)]);
            G = rmedge(G, [inid(2) ii], [ii onid(2)]);
            % mark this information to ii and sister cell
            nodes(ii).splittedNodeIDs = [ii newIndex];   
            nodes(newIndex).splittedNodeIDs = [ii newIndex];  
        end
    end

    % fix 1 frame splitting cells (only 2 cells version)
    for ii = 1:height(G.Nodes) % loop over nodes
        [~, inid] = inedges(G,ii);
        [~, onid] = outedges(G,ii);
        nodes(ii).inedges = inid;
        nodes(ii).outedges = onid;

        % nodes with 1 in edge and 2 out edges
        if length(inid) == 1 && length(onid) == 2 
            % make sure that joined nodes are not newly added from above loop
            % node ID must be within original network size otherwise causes
            % problems with masks
            if onid(1) <= sum(mask(:)) && onid(2) <= sum(mask(:))
                % make also sure that the same nodes are not used in
                % joining fixing above
                if isempty(find(onid(1) == [nodes.splittedNodeIDs])) && ...
                        isempty(find(onid(2) == [nodes.splittedNodeIDs]))
                    [~, onidDes1] = outedges(G,onid(1));
                    [~, onidDes2] = outedges(G,onid(2));
                    % 2 nodes converge back to a single node AND
                    % neither of cells have been joined before
                    if length(onidDes1) == 1 && length(onidDes2) == 1 && ...
                            onidDes1 == onidDes2 && ...
                            isempty(find(onid(1) == [nodes.joinedNodeIDs])) && ...
                            isempty(find(onid(2) == [nodes.joinedNodeIDs]))
                        
                        % remove edges to the larger node index 
                        G = rmedge(G, [ii max(onid)], [max(onid) onidDes2]);
                        % mark this information to daughter cells
                        nodes(onid(1)).joinedNodeIDs = [min(onid) max(onid)];
                        nodes(onid(2)).joinedNodeIDs = [min(onid) max(onid)];
                    end
                end
            end
        end
    end

    % save network before splitting to individual cells
    G_whole = G;

    % split whole lineage trees into single cells by removing edges whenever 
    % a node has two or more out edges/in edges
    for ii = 1:height(G.Nodes) % loop over nodes
        [~, inid] = inedges(G,ii);
        [~, onid] = outedges(G,ii);

        % 1 node becomes 2 nodes (cell division)
        % remove outedges and save parent node
        if ~isempty(onid) && length(onid) > 1
            % save daughter nodes
            nodes(ii).daughterNodeID = onid;

            for jj = 1:length(onid) % daughter nodes
                % remove edge and save parent
                G = rmedge(G, ii, onid(jj));
                nodes(onid(jj)).parentNodeID = ii;
                % update edges
                [~, inidDaugh] = inedges(G,onid(jj));
                [~, onidDaugh] = outedges(G,onid(jj));
                nodes(onid(jj)).inedges = inidDaugh;
                nodes(onid(jj)).outedges = onidDaugh;
            end 
        end 
        [~, inid] = inedges(G,ii);
        [~, onid] = outedges(G,ii);

        % 2 nodes become 1 node (incorrect segmentation)
        % remove inedges - don't save parent for these 
        if ~isempty(inid) && length(inid) > 1
            for jj = 1:length(inid) % parent nodes
                % remove edge
                G = rmedge(G, inid(jj), ii);
                % update edges
                [~, inidDaugh] = inedges(G,inid(jj));
                [~, onidDaugh] = outedges(G,inid(jj));            
                nodes(inid(jj)).inedges = inidDaugh;
                nodes(inid(jj)).outedges = onidDaugh;
            end 
        end
        [~, inid] = inedges(G,ii);
        [~, onid] = outedges(G,ii);
        
        % update edges
        [~, inid] = inedges(G,ii);
        [~, onid] = outedges(G,ii);
        nodes(ii).inedges = inid;
        nodes(ii).outedges = onid;
    end

    % find nodes corresponding to the same cell using graph distance matrix
    % all distances from node i to node j
    d = distances(G);
    % linking matrix with unique cell IDs
    linking_ID = NaN(size(linking_masks));
    ID = 1;
    for ii = 1:height(G.Nodes) % loop over nodes
        % find nodes with no parents
        if isempty(nodes(ii).inedges)
            vector = d(ii,:);
            nodeInd = find(vector < Inf);
            % cell ID
            nodes(ii).cellID = ID;
            % assign ID to right nodes
            for jj = 1:length(nodeInd)
                index = find(linking_nodes == nodeInd(jj));
                linking_ID(index) = ID;
                nodes(nodeInd(jj)).cellID = ID;
            end
            ID = ID + 1;
        end
        % save cell coordinates and mask size
        [row, col] = find(ii == linking_nodes);
        nodes(ii).maskCoord = [x_coord(row, col) y_coord(row, col)];
        nodes(ii).area = cell_areas(row, col);
        nodes(ii).frame = row;
        nodes(ii).maskID = col; 
        
        % newly created node doesn't have frame or coord,  copy from
        % duplicated one
        if isempty(row)
            % find sister index and copy the info
            sisterNodeID = min([nodes(ii).splittedNodeIDs]);
            nodes(ii).area = nodes(sisterNodeID).area;
            nodes(ii).maskCoord = nodes(sisterNodeID).maskCoord;
            nodes(ii).frame = nodes(sisterNodeID).frame;
        end
    end
    
    % save cellIDs to network
    G.Nodes.ID = [nodes.cellID]';
    G_whole.Nodes.ID = [nodes.cellID]';

    % rewrite mask files with new unique cell IDs
    disp('Re-writing mask files')
    for ii = 1:length(images_mask) % time        
        % load mask
        bw_mask = imread(fullfile(pathname,images_mask(ii).name));
        % ignore cells touching border as before
        bw_mask = imclearborder(bw_mask,4);
        % remove masks smaller than X as before
        bw_mask = bw_mask.*uint16(bwareaopen(bw_mask,size_min,4));
        % find unique mask IDs
        maskIDs = unique(bw_mask(:));
        maskIDs = double(maskIDs(maskIDs ~= 0));
        
        % make new empty bw_mask where new IDs will be written
        bw_mask_rev = bw_mask;
        bw_mask_rev(:) = 0;

        % loop over mask IDs in this frame
        for jj = 1:length(maskIDs) % cells
            % only nodes that exist in linking_nodes
            if ~isnan(linking_nodes(ii,maskIDs(jj)))     
                % check for 1 frame splitting problem
                if ~isempty(nodes(linking_nodes(ii,maskIDs(jj))).joinedNodeIDs)
                    % cells are joined by filling the mask between the cells
                    % node IDs
                    node_cell1 = nodes(linking_nodes(ii,maskIDs(jj))).joinedNodeIDs(1);
                    node_cell2 = nodes(linking_nodes(ii,maskIDs(jj))).joinedNodeIDs(2);

                    % mask numbers in bw_mask
                    mask_index1 = find(node_cell1 == linking_nodes(ii,:));
                    mask_index2 = find(node_cell2 == linking_nodes(ii,:));

                    % replace mask number with a new one from linking_ID
                    indeces = bw_mask == mask_index1 | bw_mask == mask_index2;
                    bw_mask_rev(indeces) = linking_ID(ii,maskIDs(jj));
                    
                    % distance transforms to find shortest path
                    D1 = -bwdist(bw_mask == mask_index1);
                    D2 = -bwdist(bw_mask == mask_index2);
                    D = D1+D2;
                    
                    % shortest distance value
                    D_max = max(D(:));
                    
                    % add shortest distance px to the mask
                    bw_mask_rev(D == D_max) = linking_ID(ii,maskIDs(jj));

                    % mask coordinates
                    [y, x] = find(indeces);
                    cell1_x = mean(x);
                    cell1_y = mean(y);

                    % correct area and coordinates
                    nodes(node_cell1).maskCoord = [cell1_x cell1_y];
                    nodes(node_cell2).maskCoord = NaN;
                    nodes(node_cell1).area = nodes(node_cell1).area + nodes(node_cell2).area;
                    nodes(node_cell2).area = NaN;

                    % remove info on the other mask
                    linking_nodes(ii,mask_index2) = NaN;
                    nodes(node_cell2).cellID = NaN;
                    nodes(node_cell1).joinedNodeIDs = [];
                    
                % check for 1 frame joining problems
                elseif ~isempty(nodes(linking_nodes(ii,maskIDs(jj))).splittedNodeIDs)
                    % node IDs
                    node_cell1 = nodes(linking_nodes(ii,maskIDs(jj))).splittedNodeIDs(1);
                    node_cell2 = nodes(linking_nodes(ii,maskIDs(jj))).splittedNodeIDs(2);
                    % mask number and indeces
                    mask_index = find(min([node_cell1 node_cell2]) == linking_nodes(ii,:));
                    mask_cell = bw_mask == mask_index;
                    
                    % parent node IDs are taken from G_whole that still has
                    % all edges intact before splitting to individual cells
                    [~, inid1] = inedges(G_whole,node_cell1);
                    [~, inid2] = inedges(G_whole,node_cell2);
                    
                    % previous frame masks_IDs
                    mask_ID_prev1 = nodes(inid1).cellID;
                    mask_ID_prev2 = nodes(inid2).cellID;
                    mask_cell_prev1 = bw_mask_previous == mask_ID_prev1;
                    mask_cell_prev2 = bw_mask_previous == mask_ID_prev2;
                    mask_cell_prev = bw_mask_previous == mask_ID_prev1 ...
                        | bw_mask_previous == mask_ID_prev2;
   
                    % estimate translation to align two masks using phase
                    % correlation and apply it to previous mask (no interpolation)
                    % Note! this only affects linking, not saved to masks
                    tformEstimate = imregcorr(mask_cell_prev,mask_cell,'translation');
                    Rfixed = imref2d(size(mask_cell));
                    mask_cell_prev = imwarp(mask_cell_prev,tformEstimate,'OutputView',Rfixed,'interp','nearest');
                    mask_cell_prev1 = imwarp(mask_cell_prev1,tformEstimate,'OutputView',Rfixed,'interp','nearest');
                    
                    % reduce mask size by 1 px
                    mask_cell_prev = imerode(mask_cell_prev, strel('disk', 2));
                    
                    % distance transform
                    D = -bwdist(~mask_cell);
                    % make previous mask areas deep wells for watershed
                    D(mask_cell_prev) = -100;
                    % minima transform to avoid other local minima
                    Dmin = 4;
                    wd_mask = imextendedmin(D,Dmin);
                    D2 = imimposemin(D,wd_mask);
                    % watershed
                    Ld2 = watershed(D2);
                    % split mask by watershed lines and label them
                    mask_cell(Ld2 == 0) = 0;
                    [labels, no_areas] = bwlabel(mask_cell);
                    
                    % only use 2 largest areas
                    if no_areas > 2
                        % sizes of nucleoid areas
                        areas_label = zeros(1,no_areas);
                        for ff = 1:no_areas
                            areas_label(ff) = sum(sum(labels == ff));
                        end

                        [~,indLabel] = sort(areas_label,'descend');

                        for ff = 3:length(indLabel)
                            labels(labels == indLabel(ff)) = 0;
                        end
                    elseif no_areas == 1
                        % if only one area use skel to find endings of
                        % region and use maxima to split the region
                        D = -bwdist(~mask_cell);
                        skel = bwskel(mask_cell);
                        [y, x] = find(skel);
                        % define maxima pixels of skeleton to -200 in
                        % distance transform to force 2 areas
                        if (max(x)-min(x)) > (max(y)-min(y))
                            [~,indMax] = max(x);
                            [~,indMin] = min(x);
                            D(y(indMax(1)),x(indMax(1))) = -200;
                            D(y(indMin(1)),x(indMin(1))) = -200;
                        else
                            [~,indMax] = max(y);
                            [~,indMin] = min(y);
                            D(y(indMax(1)),x(indMax(1))) = -200;
                            D(y(indMin(1)),x(indMin(1))) = -200;
                        end
                        wd_mask = imextendedmin(D,Dmin);
                        D2 = imimposemin(D,wd_mask);
                        % watershed
                        Ld2 = watershed(D2);
                        % split mask by watershed lines and label them
                        mask_cell(Ld2 == 0) = 0;
                        [labels, no_areas] = bwlabel(mask_cell);
                    end
                    
                    % find max overlap with labels for node 1
                    overlap_mask = mask_cell_prev1.*double(labels);
                    overlap_mask = overlap_mask(overlap_mask ~= 0);

                    % find label number with highest overlap for node 1
                    edges = unique(overlap_mask);
                    counts = histc(overlap_mask,edges);
                    [~, ind] = max(counts);
                    maxId = edges(ind);
                    if isempty(maxId)
                        uniq = unique(labels);
                        maxId = uniq(end);
                    end
                    
                    % find label id for node 2
                    uniq_labels = unique(labels);
                    node2label = uniq_labels(uniq_labels ~= 0 & uniq_labels ~= maxId);
                    
                    % find new cell masks
                    node1_mask = labels == maxId;
                    node2_mask = labels == node2label(1);
    
                    % rewrite the masks IDs
                    bw_mask_rev(node1_mask) = nodes(node_cell1).cellID;
                    bw_mask_rev(node2_mask) = nodes(node_cell2).cellID;

                    % add mask area and coordinate information into nodes
                    [y, x] = find(node1_mask);
                    cell1_x = mean(x);
                    cell1_y = mean(y);
                    [y, x] = find(node2_mask);
                    cell2_x = mean(x);
                    cell2_y = mean(y);
                    nodes(node_cell1).maskCoord = [cell1_x cell1_y];
                    nodes(node_cell2).maskCoord = [cell2_x cell2_y];
                    nodes(node_cell1).area = sum(node1_mask(:));
                    nodes(node_cell2).area = sum(node2_mask(:));
                else
                    % replace mask number with a new one from linking_ID
                    indeces = bw_mask == maskIDs(jj);
                    bw_mask_rev(indeces) = linking_ID(ii,maskIDs(jj)); 
                end
            end
        end
        % keep the mask for the next round of mask revision
        bw_mask_previous = bw_mask_rev;

        % save mask image with the same name
        imwrite(bw_mask_rev,fullfile(pathname,images_mask(ii).name))
    end

    % create cell ID based information from node based
    % unique cellIDs
    unique_cellIDs = unique([nodes.cellID]);
    unique_cellIDs = unique_cellIDs(~isnan(unique_cellIDs));
    ID_vector = [nodes.cellID];
    for ii = 1:length(unique_cellIDs) % cell IDs
        indeces = find(ID_vector == unique_cellIDs(ii));
        cell_info(n).XY = uniq_XY(hh);
        cell_info(n).ID = unique_cellIDs(ii);

        % sort frames
        frames = [nodes(indeces).frame]; 
        [~,indSort] = sort(frames);
        cell_info(n).birth = min(frames);
        cell_info(n).division = max(frames);

        % sort area and coords using frame indeces
        area = [nodes(indeces).area];
        cell_info(n).area = area(indSort);
        coords = [nodes(indeces).maskCoord];
        coords = [coords(1:2:end)', coords(2:2:end)'];
        cell_info(n).coord = coords(indSort,:);

        % parents
        parentNode = [nodes(indeces).parentNodeID];
        if ~isempty(parentNode) && length(parentNode) == 1
            cell_info(n).parent = nodes(parentNode).cellID;
        end

        % daughter
        daughterNodes = [nodes(indeces).daughterNodeID];
        if ~isempty(daughterNodes) && length(daughterNodes) == 2
            cell_info(n).daughters = [nodes(daughterNodes).cellID];
        end

        n = n + 1;
    end

    disp(['XY' num2str(uniq_XY(hh)) ' done'])
end

% estimate when each XY should be truncated by comparing overall growth rate
% before and after each time point
XY_limits = NaN(size(XY_areas,1),1);
for jj = 1:size(XY_areas,1)
    % calculate mean area before and after time point
    before = NaN(1,size(XY_areas,2));
    after = NaN(1,size(XY_areas,2));
    for ii = 1:size(XY_areas,2)
        before(ii) = mean(XY_areas(jj,1:ii));
        after(ii) = mean(XY_areas(jj,ii:end));
    end
    % growth rate at each point
    diff_before = diff(before);
    diff_after = diff(after);
    
    % see if there are at least window_size number of points where
    % diff_after is X fold smaller than average growth (diff_after(1))
    if sum(diff_after/diff_after(1) < growth_fold_change) >= window_size 
        % filter out time points until after is below threshold
        points = find(diff_after/diff_after(1) < growth_fold_change);
        diff_before(1:points(1)) = 0;
        [~,XY_limits(jj)] = max(diff_before);
    else
        XY_limits(jj) = size(XY_areas,2);
    end
end

% save cell_info for future use
save(fullfile(pathname,'cell_info.mat'),'cell_info','XY_limits','XY_areas')

% save network and linking information
save(fullfile(pathname,'extra_info.mat'),'G','G_whole','linking_ID',...
    'linking_masks','linking_masks_rev','linking_nodes','nodes')

% Extra: how to visualize a graph G:
% figure;h = plot(G,'Layout','force')
% highlight(h,[2240],'NodeColor','g','MarkerSize',10)

end

