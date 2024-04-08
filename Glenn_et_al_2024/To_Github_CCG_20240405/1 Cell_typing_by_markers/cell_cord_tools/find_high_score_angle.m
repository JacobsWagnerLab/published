
% Calculate the angle with higest score

function [opt_theta_degree, alert_flag] = find_high_score_angle(score_mat, focus_pos, move_dist, degree_range, minimal_level)

%score_mat = exp(score_mat);
%focus_point = ini_pixel = max_pixel(1:2);
%move_dist = 3;
%degree_range = [10 200];

alert_flag = 0;  % Alert for hitting the boundary. when hitting boundary, flag = 1
opt_theta_degree = NaN; % optimal theta (in degree)

theta_degree = degree_range(1):degree_range(2);
theta_rad = (2*pi/360)*theta_degree';
circle_point = [cos(theta_rad) sin(theta_rad)];  % row, col

tn = size(theta_rad,1);
candidate = NaN(tn, 2);
candidate_score = NaN(tn,1);

for j = 1:tn
    candidate(j,:) = focus_pos + (move_dist * circle_point(j,:));
end

% Test if the candidate point sets hit the boundary

flag = sum(min(candidate) > [1 1]) + sum(max(candidate) < size(score_mat)); 

if (flag == 4)

    alert_flag = 0;
    
    for j = 1:tn       
        candidate_score(j) = interp2(score_mat, candidate(j,2), candidate(j,1), 'linear');        
    end

    [val_j, ind_j] = max(candidate_score);
    opt_theta_degree = theta_degree(ind_j);
    
    if (val_j < minimal_level)
        alert_flag = 1;
    end

elseif (flag < 4)
    
    alert_flag = 1;
    
end

%figure; plot(theta_degree, candidate_score, 'o-');

