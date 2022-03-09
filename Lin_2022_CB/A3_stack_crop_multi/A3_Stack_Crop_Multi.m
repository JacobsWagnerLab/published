
% 2022 Feb 17th: Final version
% Cropping multiple ROI from single image stack 
% according to the marks in parallel

% ====================================================================== %
% (1) First, need to open EnterGlobalVar and define the global variables
% ====================================================================== %

% (2) Enter following variables:
InputSubFolder = 'waf7_001xy3';
InputNameTag = InputSubFolder;  % do not enter 't' in the end of prefix

% Input ROI (1) Width, (2) L1 lengths, (3) L2 lengths and (4) type pameters.
% type=1: point-up.  type=2: point-down
ROIParam = [100 80 570 1];
StackInd = [1 2000];

% Enter a list of XY_mark (n-by-2) which specify the tip points of 
% microfluidic chambers

XY_Mark = [];

% Example:
% XY_Mark = [145 198; 450 199; 602 201; 755 201];


% ====================================================================== %
%  Script starts, do not modify the script below
% ====================================================================== %

% Define ROI for cropping. 
for mj = 1:size(XY_Mark, 1)
        
   XYRg(mj,:) = ROI_generation(XY_Mark(mj,:), ROIParam);
        
end

% For each microfluidic chamber, crop into indepedent folders

cd(script_dir);

PosNum = size(XY_Mark,1);

parfor pos = 1:PosNum

    ExportSubFolder = strcat(InputSubFolder, '-', num2str(pos));
    XYRg_pos = XYRg(pos,:);
        
    [Ret] = sub_crop_one_chamber ...
    (script_dir, data_folder, InputSubFolder, ExportSubFolder, InputNameTag, ...
     StackInd, XYRg_pos, ChNum, ChannelFlag, ...
     InpFdM, InpChM, ExpChM, max_digit);

end
 
cd (main_folder);
LocSaveName = 'CropLocations.mat';            
save(LocSaveName, 'XYRg');
 
% =======================================================================
% Cropping one ROI region from image stacks

function [Ret] = sub_crop_one_chamber ...
    (ScriptDir, MainFolder, InputSubFolder, ExportSubFolder, InputNameTag, ...
     StackInd, XYRg, ChNum, ChannelFlag, ...
     InpFdM, InpChM, ExpChM, max_digit)

% --------------------------------------------------------------------
% Exmaple for input parameters: 

%CropStep = 2;
%PosTag = 'Pos0';
%ChNum = 2;
%Prefix = InputNameTag;
%InpChM = cellstr(['c1';'c2';'c3']);
%ExpChM = cellstr(['c1';'c2';'c3']);
%max_digit = 3;
%ChannelFlag = [1 1 3];  % image interval for each channel

% --------------------------------------------------------------------

cd(MainFolder);
mkdir(ExportSubFolder);

DirInp = strcat(MainFolder, InputSubFolder, '\');
DirExp = strcat(MainFolder, ExportSubFolder, '\');

Stack_head = StackInd(1);
Stack_tail = StackInd(2);
Xcrop1 = XYRg(1);
Xcrop2 = XYRg(2);
Ycrop1 = XYRg(3);
Ycrop2 = XYRg(4);

     cd(ScriptDir);
    
     for ch = 1:ChNum  % number of fluorescence chennels
        
        cd(DirExp); 
        mkdir(char(ExpChM(ch))); 
         
        for j = Stack_head : Stack_tail
                        
            % Only do crop for Cropflag = 1 
            
            CropFlag = 0;
            if (ChannelFlag(ch) == 1)
                CropFlag = 1;
            elseif (ChannelFlag(ch) > 1)
                CropFlag = (rem(j,ChannelFlag(ch)) == 1);
            end
            
            %
            
            if ( CropFlag == 1 )
                
                format_string = strcat('%0', num2str(max_digit), 'd');
                ImgInd = sprintf(format_string, j);               
                Input = strcat(DirInp, char(InpFdM(ch)), '\', InputNameTag, 't', ImgInd, char(InpChM(ch)),'.tif');
                
                A = imread(Input);

                mXcrop1 = max(Xcrop1, 1);
                mXcrop2 = min(Xcrop2, size(A, 2));
                mYcrop1 = max(Ycrop1, 1);
                mYcrop2 = min(Ycrop2, size(A, 1));
                
                A1 = A(mYcrop1:mYcrop2, mXcrop1:mXcrop2);

                NameOut = strcat(InputNameTag, 't', ImgInd, char(ExpChM(ch)), '.tif');
                ExpPath = strcat(DirExp, char(ExpChM(ch)), '\', NameOut);
    
                imwrite(A1, ExpPath, 'tiff');
                

            end
        end
     end
     
     cd(DirExp);
     
Coordinate = [Xcrop1 Xcrop2 Ycrop1 Ycrop2];
LocSaveName = 'CropLocation.mat';            
save(LocSaveName, 'Coordinate');

A = [];
A1=  [];

cd(ScriptDir);

Ret = 1;

end

% =======================================================================

% Assign marker manually and create ROI of crop. 
% Using for vertical microfluidic chambers.

function [xy_pos] = ROI_generation(mark_pos, param)

% Note that mark_pos = (x,y) in ImageJ
% which is equivalent to (col, row) in matlab.

xy_pos = [];  % [x1 x2 y1 y2];

% Constants for the channel size
w = param(1);
L1 = param(2);
L2 = param(3);
type = param(4);

% point-up channel:    use type = 1
% point-down channel:  use type = 2

%%================================================%% 

if (type == 1)  
    
    % point-up channel
    
    xy_pos(1) = mark_pos(1) - floor(w/2);
    xy_pos(2) = mark_pos(1) + floor(w/2);
    xy_pos(3) = mark_pos(2) - L1;
    xy_pos(4) = mark_pos(2) + L2;
    
elseif (type == 2)
    
    % point-down channel    
        
    xy_pos(1) = mark_pos(1) - floor(w/2);
    xy_pos(2) = mark_pos(1) + floor(w/2);  
    xy_pos(3) = mark_pos(2) - L2;
    xy_pos(4) = mark_pos(2) + L1;  
    
end

%%================================================%%    
    
end


    
