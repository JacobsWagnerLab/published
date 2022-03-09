
%{
-About-
This is a code for shift correction of a time-series image stack. 
The image stack is usually exported from time-lapse experiments.
The code provide a estimated shift of image for each frame.

2022 Feb 17th: Final version.

-Inputs-

input_folder: Folder containing the image files. Each image need to be
              in the format of "(img_prefix)-(image index)-(img_postfix)"
              where 'image index' is the number of frame (e.g. 001) 
              and 'img_prefix' and 'img_postfix' are strings.

output_folder: Folder for exporting the output file 'data_interpolate.mat'. 

img_prefix: a string for the prefix of the file name
img_postfix: = a string for the postfix of the file name

stack_rg: a row vector specify the beginning and ending indecies 
          of the image stack

max_digit:  maximal digit of the number of image stack
            used for AddZeros.m

frame_gap: number of frame between images using for shift detection

mark_X_loc:  a vector specify the range (in fraction between 0 and 1)
                in x-direction using for marker of shift detection

mark_Y_loc = a vector specify the range (in fraction between 0 and 1)
               in y-direction using for marker of shift detection

-Outputs-

(1) Shift correction result is saved into 'data_interpolate.mat' 
under output_folder. 

data_interpolate is an n-by-3 vector. 
Each row correspond to data of one frame.
Column 1 = frame number
(Column2 Column3) = shifted distance in (Row, Column) direction

(2) Input Parameters are saved into 'param_Stack_shift_Correction.mat'  
under output_folder

-Example-
input_folder = '/Users/WHLin/Desktop/Code_sharereport/20170823/Demo_dataset/';
output_folder = input_folder;
img_prefix = 'waf3_003t';
img_postfix = 'c1.tif';
stack_rg = [1 100];
max_digit = 3; 
frame_gap = 20;
mark_X_range = [1000 1500];   
mark_Y_rage = [500 800];

%}

% ====================================================================== %
%  First, need to open EnterGlobalVar and define the global variables.
% ====================================================================== %

%  Modify following variables:

% (1) Folder containing the image stack 

input_folder = main_folder;
output_folder = input_folder;

% (2) Enter the information about image stack:
%     Each image file name follows the pattern  "(img_prefix)-(image index)-(img_postfix)"

img_prefix = 'waf7_001xy3t';
img_postfix = 'c1.tif';

% (3) Enter the range (first and last index) of the image stack
%     The max_digit specify the maximal digit of image index
%     For example, a stack with 450 images has max_digit=3.

stack_rg = [1 2500]; 

% (4) Enter the frame gap between image for detecting shift

frame_gap = 1;

% (5) Specify the region of shift correction reference marker.
% input coordinate range for x (column) and y (row).

mark_X_loc = [500 650];   
mark_Y_loc = [760 865];

% (6) Specify the fraction of high/low pixel used for shift correction
pixel_frac = 0.0005;


% ====================================================================== %
%   Script starts here. Do not modify variables bellow. 
% ====================================================================== %

% (a) Define internal variables:

% Determine how many intrapolation points are used for the shift correction.
% The first and last frames are always included, so at least two points
% is used for the interpolation.
interpolate_num = ceil( (stack_rg(2) - stack_rg(1) + 1 ) / frame_gap ) + 1;  

% Array containing frame indices using for interpolation
frame_index = zeros(interpolate_num, 1);

% Get frame_index list

for k = 1: interpolate_num
    
    if (k == 1)
        
        tempJ  = 1;
    
    elseif ( ( (k>1) * (k<interpolate_num) ) == 1 )
        
        tempJ = (k-1)*frame_gap;
        
    elseif (k == interpolate_num)
        
        tempJ = (stack_rg(2) - stack_rg(1) + 1 );
    
    end
    
    frame_index(k) = stack_rg(1) + tempJ - 1;

end

% ========================================================================
% (b) Inferiing the shift correction by marker centroid

% (1) Preview of the cropped area. This is just verification of the marker
%     is within the cropped range

preview_indA = stack_rg(1);
preview_indB = stack_rg(2);

% Get image file names
format_string = strcat('%0', num2str(max_digit), 'd');
index_stringA = sprintf(format_string, preview_indA);    
index_stringB = sprintf(format_string, preview_indB);

file_name_A = strcat( img_prefix, index_stringA, img_postfix );    
file_name_B = strcat( img_prefix, index_stringB, img_postfix );
    
% Load two image files
cd(input_folder);
imgA0 = imread(file_name_A);    
imgB0 = imread(file_name_B);
        
% Crop the region with marker
imgA = imgA0( mark_Y_loc(1):mark_Y_loc(2), mark_X_loc(1):mark_X_loc(2) );    
imgB = imgB0( mark_Y_loc(1):mark_Y_loc(2), mark_X_loc(1):mark_X_loc(2) );    
    
figure;
subplot(211); imagesc(imgA);
subplot(212); imagesc(imgB);


% (2) Using the selected frames, estimating the drift by marker

% For each selected frame, find the centroid of markers

centroid_rec = zeros(interpolate_num, 2);

parfor j = 1:interpolate_num
   
    k = frame_index(j);
    k
    
    % Generate file name (string)     
    format_string = strcat('%0', num2str(max_digit), 'd');
    index_string = sprintf(format_string, k);
    file_name = strcat( img_prefix, index_string, img_postfix );
    
    % Load and crop the kth image
    cd(input_folder);    
    img_k0 = imread(file_name);
    img_k = img_k0( mark_Y_loc(1):mark_Y_loc(2), mark_X_loc(1):mark_X_loc(2) ); 
    
    [centroid] = MarkerCentroid(img_k, 1, pixel_frac);  % param set to be 1, for bright markers
    centroid_rec(j,:) = centroid;
    
end

% Calculate the shift dynamics according to the centroid
shift_array = centroid_rec;
shift_array(:,1) = shift_array(:,1) - centroid_rec(1,1);
shift_array(:,2) = shift_array(:,2) - centroid_rec(1,2);

% Interpolate data for all frames
[data_interpolate] = Interpolate(frame_index, shift_array);


% ========================================================================

% (c) Save the intrapolate shift data in ouput folder

cd (output_folder);

% Saving shift correction data
shift_data = 'data_interpolate.mat';            
save(shift_data, 'data_interpolate'); 

% Saving parameter data
param_data = 'param_Stack_shift_Correction.mat';
save(param_data, 'input_folder', 'output_folder', 'img_prefix', 'img_postfix', ...
     'stack_rg', 'max_digit', 'frame_gap', 'mark_X_loc', 'mark_Y_loc', 'pixel_frac');

 
% ====================================================================== %
%  Subscrpts
% ====================================================================== %

function [centroid] = MarkerCentroid(img, param, pixel_frac)

%{
This script find the centroid coordinate of the brightest/darkest region
of the image.

- Input -
img: a grayscale image
param: specify to find the centroid of bright region (set to 1) 
       or darkest region (set to 2)  

- Output -
centroid: the (row, col) coordinate for centroid
%}


% Columnize the image and sort the brightest/darkest pixels

vec = reshape(img, [], 1);

img_filter = img;

if (param == 1) % threshold above
    
    qH = quantile(vec, 1-pixel_frac);
    img_filter(img < qH) = 0;

elseif (param == 2)
    
    qL = quantile(vec, pixel_frac);
    img_filter(img > qL) = 0;

end

% figure;imagesc(img_filter);

% Finding the centoid of the sorted pixels

moment = [0 0 0];  % pixel count, 1st moment in row, 1st moment in column
centroid = [0 0];  

for r = 1:size(img_filter, 1)
    
    for c = 1:size(img_filter, 2)
       
        temp = double(img_filter(r,c));
        
        moment(1) = moment(1) + temp;
        moment(2) = moment(2) + r*temp;
        moment(3) = moment(3) + c*temp;                  
        
    end
    
end

centroid(1) = moment(2)/moment(1);
centroid(2) = moment(3)/moment(1);

end

% ====================================================================== %

function [interpolated_data] = Interpolate(t, x)

%{

The function perform linear interpolate the data vector 
for integer time points between first and last time point. 
This function is specilized for image stacks, and all time point, 
position vector data must have integers values.

-Inputs-

t:  n-by-1 time-series vector. Need to be ascending intergers.
    n is the number of data points availible.

x:  n-by-m position vector. Need to be intergers.
    m is the dimension of the position vector (m=2 for monochrome images)

-Outputs-
interpolated_data:  p-by-(m+1) vector, where p is the total time points 
                    between t(1) and t(end).

Each frame is one integer time point.
Column 1 = interpolated time data
Column 2 to Column (m+1) = interpolated position data

-Example-
t = [1 10 20 30 35]';
x = [55 59 52 51 55; 13 5 19 4 6]';
   
%}

% Let p be the total time points between t(1) and t(end).

% Define a p-by-1 vector for intrpolated time data
t_intra = ( t(1):t(end) )';

% Define a p-by-m vector for intrpolated position data
x_intra = zeros(t(end)-t(1)+1, size(x,2)) ;

% Define number of node of input data
node_num = size(t,1);

% Loop through every node pairs to do interpolation

for j = 1: (node_num-1)

    % Using data of node k and k+1, do intrapolation for all integer time
    % points in between. 
       
    indA = t(j);
    indB = t(j+1);
    
    posA = x(j,:);
    posB = x(j+1, :);
    
    % Do interpolation for all integer time points
    
    for q = indA:indB
       
        x_intra(q,:) = posA + round( (posB-posA)* ( (q-indA)/(indB-indA) ) );        
        
    end
    
end
   
interpolated_data = [];
interpolated_data(:,1) = t_intra;
interpolated_data(:,2:3) = x_intra;


% Plot the figures if plotflag = 1.
% (Only plot for first and second dimension)

plotflag = 1;

if (plotflag == 1)

    figure;
    
    plot(t, x(:,1), 'bo', t, x(:,2), 'ro');
%    legend('x','y');
    hold on;
    plot(t_intra, x_intra(:,1), 'b-.', t_intra, x_intra(:,2), 'r-.' );     
    hold off;
    
    legend('row', 'column', 'r-intra', 'c-intra');
    
end

end

% ====================================================================== %
