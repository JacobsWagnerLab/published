
% Finding mid-axis of cells

function [output_midline, mid_axis_output] = find_middle_axis(level_mat)
%level_mat = obj.level;

move_dist = 2;
increment_dist = 1;
max_bending_angle = 15;
minimal_level = 1;
max_step = 150;
crossing_dist = 3;
mid_line_scale = 0.1;

%%% ========================= Step 1 =========================== %%%

% Finding pixels with maximal intensity as initial point of mid_axis.

GS_param = [5 1.5];
score_mat = Img_GS_conv(level_mat, GS_param);

max_pixel = [1 1 0];  % row, col, value

for r = 1:size(level_mat,1)
    for c = 1:size(level_mat,2)
        
        if (score_mat(r,c) > max_pixel(3))  
            
            max_pixel = [r c score_mat(r,c)];
            
        end
        
    end
end

% Enhance the difference of score matrix (optional)
score_matE = score_mat;


%%% ========================= Step 2 =========================== %%%

mid_axis1 = [];
opt_theta_rec1 = [];
mid_axis2 = [];
opt_theta_rec2 = [];

% find two initial directions 
focus_pos = max_pixel(1:2);

degree_range1 = [1 360];
[opt_theta_degree1] = find_high_score_angle(score_matE, focus_pos, move_dist, degree_range1, minimal_level);

degree_range2 = opt_theta_degree1 + 180 + 90*[-1 1] ;
[opt_theta_degree2] = find_high_score_angle(score_matE, focus_pos, move_dist, degree_range2, minimal_level);


% Start from opt_theta_degree1

alert = 0;
focus_pos = max_pixel(1:2);
mid_axis1 = [mid_axis1; focus_pos];
opt_theta_rec1 = [opt_theta_rec1; NaN];
opt_theta_degree = opt_theta_degree1;

for h = 1:max_step
    
    if (alert == 0)
        opt_theta_rad = (2*pi/360) *opt_theta_degree;
        focus_pos = focus_pos + increment_dist *[cos(opt_theta_rad) sin(opt_theta_rad)];
        degree_range = opt_theta_degree + max_bending_angle*[-1 1];
        [opt_theta_degree, alert] = find_high_score_angle(score_matE, focus_pos, move_dist, degree_range, minimal_level);

        mid_axis1 = [mid_axis1; focus_pos];
        opt_theta_rec1 = [opt_theta_rec1; opt_theta_degree];
    
    elseif (alert == 1)
        break;
    end
        
end

% Start from opt_theta_degree2

alert = 0;
focus_pos = max_pixel(1:2);
mid_axis2 = [mid_axis2; focus_pos];
opt_theta_rec2 = [opt_theta_rec2; NaN];
opt_theta_degree = opt_theta_degree2;

for h = 1:max_step

    if (alert == 0)
        opt_theta_rad = (2*pi/360) *opt_theta_degree;
        focus_pos = focus_pos + increment_dist *[cos(opt_theta_rad) sin(opt_theta_rad)];
        degree_range = opt_theta_degree + max_bending_angle*[-1 1];
        [opt_theta_degree, alert] = find_high_score_angle(score_matE, focus_pos, move_dist, degree_range, minimal_level);

        mid_axis2 = [mid_axis2; focus_pos];
        opt_theta_rec2 = [opt_theta_rec2; opt_theta_degree];
    
    elseif (alert == 1)
        break;
    end
        
end

%%% ========================= Step 3 =========================== %%%

% Combine two branches of mid-axis

mid_axis = [];

L1 = size(mid_axis1, 1);
L2 = size(mid_axis2, 1);

% Reverse data in mid_axis1 
for j = 1:L1
    
    temp = mid_axis1(L1-j+1,:);
    mid_axis = [mid_axis; temp];

end

% Add data in mid_axis2
% Note that the initial point is repetitive and start from index j=2.

for j = 2:L2
    
    temp = mid_axis2(j,:);
    mid_axis = [mid_axis; temp];

end

%%% ========================= Step 4 =========================== %%%

% Remove self-crossing points

crossing_flag = NaN(L1+L2-1,1);

% (1) Caculate if crossing happens for remoted mid-axis points
%     for first branch (start from terminal)

for j = 1:L1
    
    ind = j;
    
    dist_vec = mid_axis - kron(mid_axis(ind,:), ones(size(mid_axis,1),1) );
    dist_ind = sqrt(dist_vec(:,1).^2 + dist_vec(:,2).^2);

    near_flag = (dist_ind < crossing_dist);   % points which are closed to point j
    
    for k = 1:(L1+L2-1)
       
        if ( abs(ind-k) > 2*crossing_dist ) && (near_flag(k) == 1)            
            crossing_flag(ind) = 1; 
        else
            crossing_flag(ind) = 0;
        end        
    end    
    
end

% (2) Caculate if crossing happens for remoted mid-axis points
%     for second branch (start from terminal)

for j = 1:L2
    
    ind = (L1+L2-1)-j+1;
    
    dist_vec = mid_axis - kron(mid_axis(ind,:), ones(size(mid_axis,1),1) );
    dist_ind = sqrt(dist_vec(:,1).^2 + dist_vec(:,2).^2);

    near_flag = (dist_ind < crossing_dist);   % points which are closed to point j
    
    for k = 1:(L1+L2-1)
       
        if ( abs(ind-k) > 2*crossing_dist ) && (near_flag(k) == 1)            
            crossing_flag(ind) = 1;  
        else
            crossing_flag(ind) = 0;
        end        
    end    
    
end

mid_axis_output = [];

for r = 1: (L1+L2-1)
    
    if (crossing_flag(r) == 0)        
        mid_axis_output = [mid_axis_output; mid_axis(r,:)];
    end
    
end

%%% ========================= Step 5 =========================== %%%

% Spline fit the row and col coordinate of mid_axis_output 
% to polynomial functions

L = size(mid_axis_output,1);

mid_line = [];
mid_line_s = [(-2):mid_line_scale:(L+2)]';

row_parameter = polyfit((1:L)', mid_axis_output(:,1), 3); 
mid_line(:,1) = polyval(row_parameter, mid_line_s);

col_parameter = polyfit((1:L)', mid_axis_output(:,2), 3); 
mid_line(:,2) = polyval(col_parameter, mid_line_s);

output_midline = [mid_line mid_line_s];

%%% ========================= Step 6 =========================== %%%

%{
figure;
subplot(211);   
plot(mid_axis(:,1), 'ro'); hold on;
plot(mid_line_s, mid_line(:,1), '.-');

subplot(212);  
plot(mid_axis(:,2), 'ro'); hold on;
plot(mid_line_s, mid_line(:,2), '.-');

%%%


figure; 

imagesc(score_mat);  hold on;
plot(mid_axis(:,2), mid_axis(:,1), 'ko-'); hold on;
plot(mid_line(:,2), mid_line(:,1), 'r-'); hold on;
hold off;

colormap(gray);
%}
% subplot(121);  imagesc(data);
% subplot(122);  imagesc(data2); 




