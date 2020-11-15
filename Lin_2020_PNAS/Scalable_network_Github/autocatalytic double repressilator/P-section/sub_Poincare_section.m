
% This script generate the Poincare section.
% (1) The 7-dimensional trajectory Y(t) is first projected on Y2-Y3 space
%     to obain an projected attractor.
% (2) The line segment {Y2=0.2} is chosen as the Poincare section for the 
%     projected attractor. 

% Input: 
% - param: the value of parameter, here stands for alpha
% - y: the trajectory data (already converged to attractor)


function [rec] = sub_Poincare_section(param, y)

section = 2;        % Coordinate for Poincare section
section_pos = 0.2;  % Section position. Here, choose as {Y2=0.2}
coordinate = 3;     % Another coordinate for the projected atractor, here as Y3
        
rec = [];           
% 'rec': record of coordinate of the points on Poincare section
% it is an k-by-2 array, where k is number of points on the section. 
% first column: value of alpha (constant)
% second column: projection of points on Poincare section.         

% ============== Finding points on the section ================

num = size(y,1);

for j = 1:(num-1)
    
    % Test if the trajectory segment between y(j) and y(j+1) cross the section
    
    test1 = ( y(j, section) <= section_pos ) && (  y(j+1, section) > section_pos  );     
    
    if (test1 == 1)  
        
        % Trajectory cross the section. 
        % Intrapolate the data from Y2-Y3 coordinate
        
        % Y3 component
        coordinate_temp = [y(j, coordinate) y(j+1, coordinate)];        
        
        % Y2 component
        section_temp = [y(j, section) y(j+1, section)];
        
        % intrapolation, finding the location where the segment {y(j)-y(j+1)} intersect
        % with the Poincare section
        slope = (coordinate_temp(2) - coordinate_temp(1)) / (section_temp(2) - section_temp(1));        
        y_coordinate_pos = coordinate_temp(1) + slope* (section_pos - section_temp(1) );
        
        % Mark one pair of coordinate (parameter, Y3) for this point
        mark  = [param y_coordinate_pos];
        % Appending this point to the recording array
        rec = [rec; mark];
        
    end    

end

% End of script

