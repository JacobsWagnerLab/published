
function [seg_output] = Mask_1D(img0, qt_cut)

%{
-About-

The function Img1DSeg is used for segmentation of cells in microfluidic
chamber with only 1 row (sometimes referred to as the "mother machine"). 

The input image must be a verticle chamber with cells. Input can be a
fluorescent image or a debacked and inverted phase image. As long as the
cells has brighter intensity, and the background is uniformly dimmer, 
this script should work. 

The script first conver the input image to 1D image by summing up the 
second dimension. Then, the 1D image is thresholded by a quantile cutoff,
denoted by qt_cut, resulting an 1D binary image with 0 or 1 segments.  

Then, the 0 and 1 segments were converted into segments labeled from 0 to N,
where N is the total number of cells. The segment labeled by number k>0 
represents a range containing the kth cells. The segments labeled by 0 
represent the range with no cells.

The output of the script seg_output is a matrix with original image size, 
labeled with number from 0 to N. Similarly, The region labeled by number 
k>0 represents a range containing the kth cells. The regions labeled by 0 
represent the range with no cells.


-Inputs-

img0: original image for 1D segmentation. Input image must be a 
      verticle chamber with cells, with cell has higher signal than the 
      background.  

qt_cut: The quantile cutoff threshold for the 1D segment, ranged from 0 to 1.  


-varargin-
N.A.

-Outputs-
seg_output: a data matrix labeled with number from 0 to N. Region with 0 
            represents the background, and regions with 1 to N represent
            ranges containing cell 1 to N.


-Example-

img0: a fluorescent image with cell in one row
qt_cut = 0.6

Also see the demo_main.m for exmaple

-Supplementary-
Supplemental file location

-Keywords-
1D segmentation

-Dependencies-
N.A.

-References-
WHLin lab presentation: WHLin 20180219.pptx

-Author-
Wei-Hsiang Lin, 2018 Feb 19

%}

% ================================================================ %
% Script start. Do not modifiy the content below 
% ================================================================ %

% Define the quantile fravtion for normalization.
% The lower quantile will be normalized to 0,
% and the higher quntile will be normalized to 1.

% In the following, we choose lower quantile= 0.05 and 
% higher quantile to be 0.9.

low_qt_frac = 0.05;
high_qt_frac = 0.9;

% (0) Converting the original image to 1D binary image
img1D = sum(img0, 2);
L = size(img1D,1);


% (1) Normalize the 1D intensity profile by low/high quantile
qt_para = [quantile(img1D,low_qt_frac) quantile(img1D,high_qt_frac)];
img1D_norm = (img1D - qt_para(1)) / (qt_para(2)-qt_para(1));

% Segment the 1D image by qt_cut threshold
% note that thresh_1D = 0 for boundaries, 
%           thresh_1D = 1 for cells.

thresh_1D = (img1D_norm > qt_cut);

% (2) Let ind_mat be the index-recording matrix 
%     where first column are starting points of 1-segments
%     and the second column are ending points of 1-segments

ind_mat = [];
rowflag = 1;

for j = 1:L
    
    % Check the boundary j=1

    if (j == 1)        
        if ( thresh_1D(j) ==0 )  % start with a boundary
            ind_mat(rowflag, 1) = 1;
        end
    end
    
    % Check the interior 1<j<L
        
    if ( (j > 1) && (j < L) )
    
        if ( thresh_1D(j) == 0 ) &&  ( thresh_1D(j+1) == 1 )
            ind_mat(rowflag, 2) = j;
            rowflag = rowflag + 1;   % finish one segment            
        end
        
        if ( thresh_1D(j) == 1 ) &&  ( thresh_1D(j+1) == 0 )
            ind_mat(rowflag, 1) = j;
        end
        
    end
        
     % Check the boundary j=L
     
    if (j == L)
        
        if ( thresh_1D(j) == 0 )   % finish with boundary
            ind_mat(rowflag, 2) = L;
        end
        
        if ( thresh_1D(j) == 1 )   % an incompleted cell
                                   % mark the last entry as boundary                                   
            ind_mat(rowflag, 1) = L;
            ind_mat(rowflag, 2) = L;
            
        end
    end
    
end

% Define the partition point c_k as midpoints of each boundary segments
c_ind = ( ind_mat(:,1) + ind_mat(:,2) ) / 2;
c_ind = floor(c_ind);

% Now, use the partition points vector {c_k} to assign number from 0 to c_n 
% to the segmentation index-vector 

seg_1D_vec = 0 *img1D;

for m = 1: (size(c_ind,1) - 1)

    seg_1D_vec(c_ind(m):c_ind(m+1),1) = m*ones(c_ind(m+1)-c_ind(m)+1,1);
    
end

% Generate a datamatrix with original size, with index 0 to c_n. 
seg_output = kron( seg_1D_vec, ones(1 ,size(img0,2)) );


% This is the end of the script

%{
figure;
subplot(311); plot(img1D_norm');
subplot(312); plot(thresh_1D');
subplot(313); plot(seg_1D_vec');
%}
