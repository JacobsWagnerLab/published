
% Measure intensity of cell objects (with background subtracted)

function fluo_rec = measure_intensity(segmatN, path, paramA, paramB)

fn_range = paramA.Nrange;

c2_path = path.c2;
c3_path = path.c3;
name = paramB.name;
frame_interval = paramB.frame_interval; 
max_digit = paramB.max_digit;

% Estimate background fluorescence and perform fluorescence statistics

NUbg_width = 8;  % Pixel number at boundary used to estimate background

fluo_rec = {};

for fn = fn_range(1):fn_range(2)

    % Check if the c2, c3 fluorescence channel exist 

    remflag = rem(fn, frame_interval);
    
    data = {};
    
    if (remflag == 1) || (frame_interval == 1)
        
        % Load c2, c3 data        
        c2_filename = strcat(name.c2.prefix, AddZeros(max_digit, fn), name.c2.postfix);
        c3_filename = strcat(name.c3.prefix, AddZeros(max_digit, fn), name.c3.postfix);        
        cd(c2_path); data.c2 = double(imread(c2_filename));        
        cd(c3_path); data.c3 = double(imread(c3_filename));
        
        % Remove background         
        data.c2I = NUBg(data.c2, NUbg_width);        
        data.c3I = NUBg(data.c3, NUbg_width);
                
        % Perform fluorescence statistics         
        mask = segmatN(:,:,fn);
        fluo_rec{fn}.c2 = ObjIntV3(mask, data.c2I, fn);
        fluo_rec{fn}.c3 = ObjIntV3(mask, data.c3I, fn);
    
    end


end


end


%========================================================================

function [ObjRec] = ObjIntV3(FilteredMat, DataMat, Index)

ObjNum = max(max(FilteredMat));
ObjRec= zeros(ObjNum, 7); 

% 1st column: object area
% 2nd column: first moment of intensity distribution
% 3rd column: second moment of intensity distribution
% 4th column: average intensity of this object
% 5th column: SD of internsity of this object 
% 6th column: index j
% 7th column: frame index

mat_size = size(FilteredMat);

for row = 1: mat_size(1)
    for col = 1: mat_size(2)
        
        if (FilteredMat(row,col) > 0) 
            ObjRec(FilteredMat(row,col),1) = ObjRec(FilteredMat(row,col),1) + 1;
            ObjRec(FilteredMat(row,col),2) = ObjRec(FilteredMat(row,col),2) + DataMat(row,col);
            ObjRec(FilteredMat(row,col),3) = ObjRec(FilteredMat(row,col),3) + (DataMat(row,col)^2);
        end
        
    end
end

for obj = 1:ObjNum
    
    ObjRec(obj,4) = ObjRec(obj,2)/ObjRec(obj,1);
    ObjRec(obj,5) = sqrt( ( ObjRec(obj,3)/ObjRec(obj,1) ) - ( ObjRec(obj,2)/ObjRec(obj,1) )^2 ); 
    ObjRec(obj,6) = obj;
    ObjRec(obj,7) = Index;    
end

end

% =======================================================================

function [img_corrected] = NUBg(img, bd_width)

% Using left and right rims to estimate the y-dependent background

left_bd = img(:,1:bd_width);
right_bd = img(:, end-bd_width+1:end);
y_background = mean([left_bd right_bd], 2); % estimated background vector (has y-dependence)

NU_background = kron(y_background, ones(1,size(img,2)));  % estimated background matrix   
img_corrected = max(img - NU_background, 1);              % background-subtracted image

end

% =======================================================================

function [StrOutput] = AddZeros(digit, j)
%Adding zero before number and generate string
%Example: digit=3, j=132;

StrOutput = '';

% test how many digit for current string
StrTemp = num2str(j);
StrSize = size(StrTemp);    
ZeroCount = digit - StrSize(2);

if (ZeroCount>0)
    
    StrOutput = StrTemp;
    
    for m = 1:ZeroCount        
        StrOutput = strcat('0',StrOutput);        
    end
    
else
    
    StrOutput = StrTemp;
    
end

end
