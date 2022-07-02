%% Calculate the distance ratios
% -Input-
%   cellList: Oufti cellList with 'che' and/or 'gfp' fields
% -Output-
%   ratios: a struct with 4 different fields:
%       1. che_include_edge: calculation based on the che puncta, distance
%       ratios calculated based on the entire cell length
%       2. che_exclude_edge: calculation based on the che puncta, distance
%       ratios calculated between the 1st and last punctum
%       3. gfp_include_edge: calculation based on the gfp puncta, distance
%       ratios calculated based on the entire cell length
%       4. gfp_exclude_edge: calculation based on the gfp puncta, distance
%       ratios calculated between the 1st and last punctum
%       5. che_deviation_stats_include_edge:
%       stats on how much the distance ratios deviate from uniform
%       distribution, calculated based on the entire cell length
%       6. che_deviation_stats_exclude_edge:
%       stats on how much the distance ratios deviate from uniform
%       distribution, calculated between the 1st and last punctum
%       7. gfp_deviation_stats_include_edge:
%       stats on how much the distance ratios deviate from uniform
%       distribution, calculated based on the entire cell length
%       8. gfp_deviation_stats_exclude_edge:
%       stats on how much the distance ratios deviate from uniform
%       distribution, calculated between the 1st and last punctum
function ratios = calculate_distance_ratio(cellList)

all_ratio_che_include_edge = cell(length(cellList.meshData),1);
all_ratio_che_exclude_edge = cell(length(cellList.meshData),1);
all_ratio_gfp_include_edge = cell(length(cellList.meshData),1);
all_ratio_gfp_exclude_edge = cell(length(cellList.meshData),1);
all_che_stats_include_edge = cell(length(cellList.meshData),1);
all_che_stats_exclude_edge = cell(length(cellList.meshData),1);
all_gfp_stats_include_edge = cell(length(cellList.meshData),1);
all_gfp_stats_exclude_edge = cell(length(cellList.meshData),1);
for frame = 1:length(cellList.meshData)
    ratio_che_include_edge = cell(length(cellList.meshData{frame}),1);
    ratio_che_exclude_edge = cell(length(cellList.meshData{frame}),1);
    ratio_gfp_include_edge = cell(length(cellList.meshData{frame}),1);
    ratio_gfp_exclude_edge = cell(length(cellList.meshData{frame}),1);
    che_stats_include_edge = cell(length(cellList.meshData{frame}),1);
    che_stats_exclude_edge = cell(length(cellList.meshData{frame}),1);
    gfp_stats_include_edge = cell(length(cellList.meshData{frame}),1);
    gfp_stats_exclude_edge = cell(length(cellList.meshData{frame}),1);
    for cc = 1:length(cellList.meshData{frame})
        if isfield(cellList.meshData{frame}{cc},'che') && ...
                isfield(cellList.meshData{frame}{cc}.che,'spotPosition_rel') && ...
                size(cellList.meshData{frame}{cc}.che.spotPosition_rel,1) > 1
            ratio_che_include_edge{cc} = get_distance_ratio(cellList.meshData{frame}{cc},'che',true);
            ratio_che_exclude_edge{cc} = get_distance_ratio(cellList.meshData{frame}{cc},'che',false);
            che_stats_include_edge{cc} = get_deviation_stats(ratio_che_include_edge{cc});
            che_stats_exclude_edge{cc} = get_deviation_stats(ratio_che_exclude_edge{cc});
        end
        if isfield(cellList.meshData{frame}{cc},'gfp') && ...
                isfield(cellList.meshData{frame}{cc}.gfp,'spotPosition_rel') && ...
                size(cellList.meshData{frame}{cc}.gfp.spotPosition_rel,1) > 1
            ratio_gfp_include_edge{cc} = get_distance_ratio(cellList.meshData{frame}{cc},'gfp',true);
            ratio_gfp_exclude_edge{cc} = get_distance_ratio(cellList.meshData{frame}{cc},'gfp',false);
            gfp_stats_include_edge{cc} = get_deviation_stats(ratio_gfp_include_edge{cc});
            gfp_stats_exclude_edge{cc} = get_deviation_stats(ratio_gfp_exclude_edge{cc});
        end
    end
    all_ratio_che_include_edge{frame} = vertcat(ratio_che_include_edge{:});
    all_ratio_che_exclude_edge{frame} = vertcat(ratio_che_exclude_edge{:});
    all_ratio_gfp_include_edge{frame} = vertcat(ratio_gfp_include_edge{:});
    all_ratio_gfp_exclude_edge{frame} = vertcat(ratio_gfp_exclude_edge{:});
    all_che_stats_include_edge{frame} = vertcat(che_stats_include_edge{:});
    all_che_stats_exclude_edge{frame} = vertcat(che_stats_exclude_edge{:});
    all_gfp_stats_include_edge{frame} = vertcat(gfp_stats_include_edge{:});
    all_gfp_stats_exclude_edge{frame} = vertcat(gfp_stats_exclude_edge{:});
end

ratios.che_include_edge = vertcat(all_ratio_che_include_edge{:});
ratios.che_exclude_edge = vertcat(all_ratio_che_exclude_edge{:});
ratios.gfp_include_edge = vertcat(all_ratio_gfp_include_edge{:});
ratios.gfp_exclude_edge = vertcat(all_ratio_gfp_exclude_edge{:});

ratios.che_deviation_stats_include_edge = vertcat(all_che_stats_include_edge{:});
ratios.che_deviation_stats_exclude_edge = vertcat(all_che_stats_exclude_edge{:});
ratios.gfp_deviation_stats_include_edge = vertcat(all_gfp_stats_include_edge{:});
ratios.gfp_deviation_stats_exclude_edge = vertcat(all_gfp_stats_exclude_edge{:});

% Get the stats of deviations of distance ratios from 1
    function res = get_deviation_stats(ratios)
        deviations = abs(ratios - 1);
        res = zeros(1, 5);
        res(1) = min(deviations);
        res(2) = max(deviations);
        res(3) = median(deviations);
        res(4) = mean(deviations);
        res(5) = std(deviations);
    end

% Helper function for calculating the distance
    function res = calculate_interpuncta_distance(pos)
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

    function ratio = get_distance_ratio(cellStruct,puncta_color,include_edge)
        separation = calculate_interpuncta_distance(cellStruct.(puncta_color).spotPosition_rel);
        number_of_puncta = size(cellStruct.(puncta_color).spotPosition_rel,1);
        if include_edge
            ratio = separation./(cellStruct.length/(number_of_puncta+1));
        else
            % Calculate length of the sum of edge length
            max_pos = max(cellStruct.(puncta_color).spotPosition_rel(:,1));
            min_pos = min(cellStruct.(puncta_color).spotPosition_rel(:,1));
            edgeSum = cellStruct.length - max_pos + min_pos;
            ratio = separation./((cellStruct.length-edgeSum)/(number_of_puncta-1));
        end
    end
end