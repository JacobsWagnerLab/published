
% 2022 Feb 17th: Final version. This code is doing: 
% (a) remove the microfluidic chamber feature in phase images and inverting
% the signal (b) crop the fluorescence images accordingly with step (a)

% ======================================================================= %
% Enter following variables 
% ======================================================================= %

ref_dir = strcat(data_folder);
SubFolderName = 'waf7_001xy3';
ref_path = strcat(ref_dir, SubFolderName, '\ref.tif');

% (i) Reference marker on ref.tif
RefMark_XY = [297 199];
% (ii) Channel type for image crop
ChannelType = 1;  % Microfluidic chamber. type=1: point up; type=2: point down.
                       
% (iii) Data marker:  A: first frame,  B: last frame
DataNamePre = strcat(SubFolderName, 't');
StackRg = [1 2000];
DataMarkIni_XY = [51 81];  % Marker for each chamber. Use (x,y) in ImageJ

% Specified the subfolders
SubFolderPrefix = strcat(SubFolderName, '-');
SubFolderList = [1:10]'; % [1 2 3 4 5 6 7 8 9]';

Range_row = [20 560];  % short side, long side 
Range_col = [30 30];   % left, right

% ====================================================================== %
%  Scrip start, do not modify the script below
% ====================================================================== %

%  Combine above parameters into ParamA and ParamB lists.

if (ChannelType == 1)
    RgR = [Range_row(1) Range_row(2)];  % up, down directions
    RgC = Range_col;   % left, right directions
elseif (ChannelType == 2)
    RgR = [Range_row(2) Range_row(1)];  % up, down directions
    RgC = Range_col;   % left, right directions
end

RefMark = [RefMark_XY(2) RefMark_XY(1)];  %convert to (Row, Col) in Matlab
DataMarkIni = [DataMarkIni_XY(2) DataMarkIni_XY(1)];
LowerLimit = 5000;

ParamA = [ChannelFlag max_digit];
ParamB = [RgR RgC DataMarkIni StackRg RefMark LowerLimit];

% ----------------------------------- %
% Call multiple matlab using parfor
% ----------------------------------- %

parfor j = 1:length(SubFolderList)

    TagJ = num2str( SubFolderList(j));
    PathJ = strcat(data_folder, SubFolderPrefix, TagJ, '\' );
    
    cd(script_dir);
    
    [RetJ] = CallStack_debackV3B(script_dir, PathJ, ref_path, DataNamePre, ParamA, ParamB, ChNum);
        
end


