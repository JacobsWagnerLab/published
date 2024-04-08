
% This code map cell cordinate of a given spot.
%cell_id = 85;
%data = spot_data2.cell_cord{cell_id};
%spot_id = data{1}.spot_index;
%spot_ensemble = spot_data2.ensemble;

function [ind_flag_percentile] = map_spot_on_cell_cord(cell_id, spot_id, tile, spot_ensemble)

% From cell id, get the cell coordinate (relative)
cell_cord_data = tile.cell_cord{cell_id};

cell_corner = cell_cord_data.corner;    % marking absolute position, {r,c}
cell_midline = cell_cord_data.midline;  % relative position wrt to corner, {r,c} 

spot_location = spot_ensemble.location;
temp = spot_location(spot_id).Centroid;  % {c,r} coordinate

spot_relative_location = [temp(2) temp(1)] - [cell_corner(1) cell_corner(2)] + [1 1];  %relative {r,c}

% Finding the nearest cell midline point

dist_flag = norm(spot_relative_location - cell_midline(1,1:2));
ind_flag = 1;

for j = 1:size(cell_midline, 1)
    
    distJ = norm(spot_relative_location - cell_midline(j,1:2));
    
    if (distJ < dist_flag)
        
        dist_flag = distJ;
        ind_flag = j;
    end
    
end

% Calculate the percentile of nearest position
ind_flag_percentile = ind_flag/size(cell_midline,1);

%[dist_flag ind_flag]

%figure;
%plot(cell_midline(:,2), cell_midline(:,1), 'b-');  hold on;
%plot(spot_relative_location(2), spot_relative_location(1), 'ro');


