
function [opt_deback, opt_data_mark, simple_deback] ...
       = Optimal_DebackV2(data_img, ref_img, data_mark, ref_mark, L, row_range, col_range)

%{
-About-

This is a code for background subtraction, useful in microfluidic images.
User assign two marks on dataerence image and ref image, respectively. 
The algorithm try to search for optimal background subtraction by 
perturbing the dataerence mark locally. The optimal result is defined by 
minimizing the absolute value of total image gradient of the debacked image.

This algorithm does not return the entire deback image as output. 
It only return a subset of image adjacent to the ref_img mark. 
The user need to specify proper a range adjacent to the mark that 
does not exceed the size of the image. Otherwise, the script return
an error flag.
       
The script also return a simple deback result (simple_deback) using initial
ref and dataerence marks without image gradient minimization.
       

-Inputs-

data_img: dataerence image containing the background feature appears in ref image.

ref_img: ref image. 

data_mark: a 1x2 vector specify the initial dataerence mark (in (row,col) coordinate)
          Need to be positive integers.

ref_mark: a 1x2 vector specify the ref mark (in (row,col) coordinate)
          Need to be positive integers.

L: the maximal distance (in pixel) for local perturbation of dataerence mark.
   Need to be positive integer.

row_range: a 1x2 vector specify the range for creating ROI of deback image.
           ROI of rows will be set by (mark-row_range(1) to mark+row_range(2))

col_range: a 1x2 vector specify the range for creating ROI of deback image.
           ROI of columns will be set by (mark-col_range(1) to mark+col_range(2))

-varargin-

N.A.

-Outputs-

opt_deback: the optimal debacked image, with size 
            (row_range(2)+row_range(1)+1, col_range(2)+col_range(1)+1)

opt_data_mark: a 1x2 vector specify the optimal dataerence mark
              which minimize the absolute value of gradient of debacked
              image.

simple_deback: the simple deback result using initial mark positions
               without image gradient minimization.


-Example-

data_img, ref_img: ref matricies
data_mark = [10 15]; 
ref_mark = [20 25];
L = 4;
row_range = [5 5];
col_range = [5 5];

Also see the demo_main.m for exmaple

-Supplementary-
Supplemental file location

-Keywords-

background subtraction
deback

-Dependencies-

N.A.

-dataerences-

WHLin lab presentation: WHLin 20170823.pptx

-Author-
Wei-Hsiang Lin, 2017 Aug 23

%}

       
% ====================================================================== %
%   Script starts here. 
% ====================================================================== %
       
% (a) Test if the image ROI cropping exceeds the size of ref or dataerence images 
%     if so, return an error flag and terminate the script.
%     if not, proceed to the next step.

% Test if the ROI of ref image can be set properly
error_flag_ref = [];
error_flag_ref(1) = ( ref_mark(1)-row_range(1) <= 0 ) ;
error_flag_ref(2) = ( ref_mark(1)+row_range(2) > size(ref_img, 1) ) ;
error_flag_ref(3) = ( ref_mark(2)-col_range(1) <= 0) ;
error_flag_ref(4) = ( ref_mark(2)+col_range(2) > size(ref_img, 2) ) ;

% Test if the ROI of dataerence image can be set properly
error_flag_data = [];
error_flag_data(1) = ( data_mark(1)-row_range(1)-L <= 0 ) ;
error_flag_data(2) = ( data_mark(1)+row_range(2)+L > size(data_img, 1) ) ;
error_flag_data(3) = ( data_mark(2)-col_range(1)-L <= 0) ;
error_flag_data(4) = ( data_mark(2)+col_range(2)+L > size(data_img, 2) ) ;


if ( sum(error_flag_ref) > 0 )
    
    %  ROI of ref image exceed the image boundary
    disp('Error_flag: the ROI region specified exceed the size of ref image');
    opt_deback = [];
    opt_data_mark = [];
    simple_deback = [];
    
    error_flag_ref
    
elseif ( sum(error_flag_data) > 0 )
    
    %  ROI of dataerence image exceed the image boundary
    disp('Error_flag: the ROI region specified exceed the size of dataerence image');
    opt_deback = [];
    opt_data_mark = []; 
    simple_deback = [];
    
    error_flag_data
    
else
    
% (b) After testing the ROI cropping, start generate ROI for ref image.

ref_ROI = ref_img(ref_mark(1)-row_range(1) : ref_mark(1)+row_range(2), ... 
                  ref_mark(2)-col_range(1) : ref_mark(2)+col_range(2) ) ;

% Define the matrix for recording the absolute image gradient for different
% candidate deback images.
gradient_record = NaN*ones((2*L)+1, (2*L)+1);  

% Test different deback images by shifting the dataerence mark 
% along row and column with respect to its initial position

for col_shift = (-L) : L
    
    for row_shift = (-L) : L       
        
        % temporal dataerence mark after local shift
        data_mark_temp = data_mark + [row_shift col_shift];        
        
        % temporal ROI of dataerence image
        data_ROI_temp = data_img( (data_mark_temp(1)-row_range(1)) : (data_mark_temp(1)+row_range(2)), ...
                                (data_mark_temp(2)-col_range(1)) : (data_mark_temp(2)+col_range(2)) );
        
        % temporal debacked image
        deback_img =  ref_ROI - data_ROI_temp; 
        
        % calculate the image gradient of the temporal debacked image
        Gmag = ImgBD(deback_img, 5, 0);
        
        % save the value of absolute image gradient in the recording matrix
        gradient_record(row_shift +(L+1), col_shift +(L+1) ) = sum(sum(abs(Gmag)));         

    end
    
end

% figure; imagesc(gradient_record);

% Finding the 'optimal shift' that minimizes the total gradient of debacked
% image. 

[opt_row_shift, opt_col_shift] = find(gradient_record == min(min(gradient_record)) );

% Find the optimal dataerence mark position correspond to the optimal shift
opt_data_mark = data_mark + [opt_row_shift-(L+1) opt_col_shift-(L+1)];

% Find the optimal ROI of dataerence image according to the optimal dataerence mark position
opt_data_ROI = data_img( (opt_data_mark(1)-row_range(1)) : (opt_data_mark(1)+row_range(2)), ...
                         (opt_data_mark(2)-col_range(1)) : (opt_data_mark(2)+col_range(2)) );

% Calculating the optimal debacked image                   
opt_deback = ref_ROI - opt_data_ROI;


% (c) The next step is calculating the debacked image with original marks
% without image gradient minimization. This is for comparison purpose.

% Calculating the ROI of dataerence image according to the original
% dataerence mark

simple_data_ROI = data_img( (data_mark(1)-row_range(1)) : (data_mark(1)+row_range(2)), ...
                            (data_mark(2)-col_range(1)) : (data_mark(2)+col_range(2)) );  
                      
% Calculating the debacked image by original dataerence mark
simple_deback = ref_ROI - simple_data_ROI;


end




end


% ===================================================================== 

function BdMat = ImgBD(InpMat, Dsize, flag)

% May 19, 2012

% Dsize must be odd integer larger than 3. 
% flag = 0 : Abs of gradient. 
% flag = 1 : x-component of gradient.
% flag = 2 : y-component of gradient.

Lx = size(InpMat,1);  Ly = size(InpMat,2);
HDsize = floor(Dsize/2);  % half diameter size

% Generate gradient matrices, used for convolution in the next step
DMx = kron(ones(Dsize,1), [-HDsize:HDsize]);
DMy = DMx';

NF = Dsize*Dsize*(HDsize+1)*(HDsize)/3;
% Normalized factor of matrices: 
% Let n = HDsize, this equal to 2*(2n+1)*[1^2 + 2^2 + ...+n^2]
% simnplify to (2n+1)(2n+1)(n+1)(n/3)

Gx0 = (1/NF)*conv2(DMx,InpMat); 
Gy0 = (1/NF)*conv2(DMy,InpMat);

% Trimmed gradient matrices that consistent to InpMat
Gx = Gx0( HDsize+1 : Lx+HDsize, HDsize+1 : Ly+HDsize);
Gy = Gy0( HDsize+1 : Lx+HDsize, HDsize+1 : Ly+HDsize);  

if (flag == 0)
    BdMat = sqrt(Gx.^2 + Gy.^2);

elseif (flag == 1)    
    BdMat = Gx;

elseif (flag == 2)
    BdMat = Gy;

end    

end

% ===================================================================== 
