
% Extract given quantity of cell cycle array

function [output] = cell_cycle_get_Rsc_traj(cc_array)

% Example:
%tagA = 'Rsc_stat';
%tagB = 'mean';

output = {};

for j = 1:length(cc_array)
    
    output{j} = cc_array{j}.Rsc_traj;
    
end
