
% 2022 Feb 17th. Final version. 
% This code crops a region of a image stack (tif files).
% Supporting shift correction (need to specifiy the shift correction data)
% which is generated by "Stack_shift_correction" algorithm

% ======================================================================= %
% (1) First, need to open EnterGlobalVar and define the global variables
% ======================================================================= %

% (2) Enter following variables:

% (a) Enter nametag for croping. This will generate a new folder with this postfix 
ExportSubFolder = strcat(input_subfolder, 'xy3');

% (b) Prefix of the file (the common part of the file name)
Prefix = 'waf7_001xy3';  % do not add 't' before index

% (c) Index range of the image stack
Stack_Rg = [1 2500];

% (d) region to be cropped
% Note that when using shift correction, if the image has maximal shift Dmax,
% need to reserve a distaince of Dmax pixels to the boundary
Xcrop = [25 1950];
Ycrop = [875 1800];

% (e) Using shift cirrection (set to 1) or not (set to 0)
shift_correction = 1;

% (f) Crop type: 
% From single folder to single folder, set =1
% From single folder to multiple folders (different fluos), set =2 
% From multiple folders to multiple folders (different fluoes), set =3
CropType = 2; 


% ======================================================================= %
% (3) Script starts. Do not modify script below
% ======================================================================= %

cd(data_folder);
mkdir(ExportSubFolder);
DirInp = strcat(data_folder, input_subfolder, '\');
DirExp = strcat(data_folder, ExportSubFolder, '\');

% ================================================= %
%  (3-1) Preporcessing: shift correction for image    
% ================================================= %

% Import shift coorrection data 

if ( shift_correction == 1 )
    
    cd(DirInp);
    temp0 = load('data_interpolate.mat');
    data_intra = temp0.data_interpolate;
    
elseif(shift_correction == 0)
    
    data_intra = zeros(Stack_Rg(2), 3); 
    
end

% Test if the linear correction is within the range of entire figure
ImgIndIni = AddZeros(max_digit, Stack_Rg(1));            

if (CropType == 2)
    Input = strcat(DirInp, Prefix, 't', ImgIndIni, char(InpChM(1)),'.tif');
end

if (CropType == 3)
    Input = strcat(DirInp, 'c1/', Prefix, 't', ImgIndIni, char(InpChM(1)),'.tif');
end

AIni = imread(Input);
Lmin = [1 1];
Lmax = size(AIni);

shift_min = min(data_intra( Stack_Rg(1):Stack_Rg(2), 2:3 ) );  % (row,col) = (y,x)
shift_max = max(data_intra( Stack_Rg(1):Stack_Rg(2), 2:3 ) );

BdCheck = zeros(2,2);
BdCheck(1,1) = ( ( Xcrop(1) + shift_min(2) ) >= Lmin(2) );
BdCheck(1,2) = ( ( Xcrop(2) + shift_max(2) ) <= Lmax(2) );
BdCheck(2,1) = ( ( Ycrop(1) + shift_min(1) ) >= Lmin(1) );
BdCheck(2,2) = ( ( Ycrop(2) + shift_max(1) ) <= Lmax(1) );

% CheckSum: flag which report if the cropping violates the boundary
% If CheckSum = 4, all cropping are normal.

CheckSum = sum(sum(BdCheck));

if ( CheckSum < 4 )

    [CheckSum]  % Report error
        
elseif (CheckSum == 4)   % Start stack cropping

% ================================================= %
% (3-2) Cropping the image stacks
% ================================================= %

stackN = (Stack_Rg(2)-Stack_Rg(1)+1) ;

CorrX = data_intra(:,3);   % x <-> col
CorrY = data_intra(:,2);   % y <-> row

if (CropType == 1)   %% Crop from single folder to single folder

    for ch = 1:ChNum  % number of chennels
    
        parfor j = Stack_Rg(1) : Stack_Rg(2)
            
            % Only do crop for Cropflag = 1 
            
            CropFlag = 0;
            if (ChannelFlag(ch) == 1)
                CropFlag = 1;
            elseif (ChannelFlag(ch) > 1)
                CropFlag = (rem(j,ChannelFlag(ch)) == 1);
            end
            
            %
            
            if ( CropFlag == 1 )
              
                ImgInd = AddZeros(max_digit, j);                
                Input = strcat(DirInp, Prefix, char(InpChM(ch)),'t', ImgInd, '.tif');
                A = imread(Input);
                
                X1j = Xcrop(1) + CorrX(j);
                Y1j = Ycrop(1) + CorrY(j);
                X2j = Xcrop(2) + CorrX(j);
                Y2j = Ycrop(2) + CorrY(j);
                
                Aj = A(Y1j:Y2j, X1j:X2j);
                
                NameOut = strcat(Prefix, 't', ImgInd, char(InpChM(ch)), '.tif');
                ExpPath = strcat(DirExp, '/', NameOut);
                
                imwrite(Aj, ExpPath, 'tiff');
            
            end
        end
    end
    
elseif (CropType == 2)  %% Crop from single folders to multiple folders
            
    cd(script_dir);
    
     for ch = 1:ChNum  % number of chennels
        
        cd(DirExp); 
        mkdir(char(ExpChM(ch))); 
         
        parfor j = Stack_Rg(1) : Stack_Rg(2)
                        
            % Only do crop for Cropflag = 1 
            
            CropFlag = 0;
            if (ChannelFlag(ch) == 1)
                CropFlag = 1;
            elseif (ChannelFlag(ch) > 1)
                CropFlag = (rem(j,ChannelFlag(ch)) == 1);
            end
            
            %
            
            if ( CropFlag == 1 )
                
                ImgInd = AddZeros(max_digit, j);                
                Input = strcat(DirInp, Prefix, 't', ImgInd, char(InpChM(ch)),'.tif');
                
                A = imread(Input);
            
                X1j = Xcrop(1) + CorrX(j);
                Y1j = Ycrop(1) + CorrY(j);
                X2j = Xcrop(2) + CorrX(j);
                Y2j = Ycrop(2) + CorrY(j);
                
                Aj = A(Y1j:Y2j, X1j:X2j);

                NameOut = strcat(Prefix, 't', ImgInd, char(ExpChM(ch)), '.tif');
                ExpPath = strcat(DirExp, char(ExpChM(ch)), '/', NameOut);
    
                imwrite(Aj, ExpPath, 'tiff');
                

            end
        end
     end
     
     
elseif (CropType == 3)  %% Crop from multiple folders to multiple folders

     cd(script_dir);
    
     for ch = 1:ChNum  % number of chennels
        
        cd(DirExp); 
        mkdir(char(ExpChM(ch))); 
         
        parfor j = Stack_Rg(1) : Stack_Rg(2)
                     
            % Only do crop for Cropflag = 1 
            
            CropFlag = 0;
            if (ChannelFlag(ch) == 1)
                CropFlag = 1;
            elseif (ChannelFlag(ch) > 1)
                CropFlag = (rem(j,ChannelFlag(ch)) == 1);
            end
            
            %
            
            if ( CropFlag == 1 )
                
                ImgInd = AddZeros(max_digit, j);                
                Input = strcat(DirInp, char(InpChM(ch)), '\', Prefix, 't', ImgInd, char(InpChM(ch)),'.tif');
                
                A = imread(Input);

                X1j = Xcrop(1) + CorrX(j);
                Y1j = Ycrop(1) + CorrY(j);
                X2j = Xcrop(2) + CorrX(j);
                Y2j = Ycrop(2) + CorrY(j);
                
                Aj = A(Y1j:Y2j, X1j:X2j);

                NameOut = strcat(Prefix, 't', ImgInd, char(ExpChM(ch)), '.tif');
                ExpPath = strcat(DirExp, char(ExpChM(ch)), '/', NameOut);
    
                imwrite(Aj, ExpPath, 'tiff');
                

            end
        end
     end
     
     cd(DirExp);
    
end
   
cd(DirExp);

Coordinate = [Xcrop(1) Xcrop(2) Ycrop(1) Ycrop(2) Stack_Rg];
LocSaveName = 'CropLocation.mat';            
save(LocSaveName, 'Coordinate');

A = [];
Aj=  [];

cd(script_dir);

end   % checksum end


% ====================================================================== %
%   Subscript
% ====================================================================== %

function [StrOutput] = AddZeros(digit, j)

% Adding zero before number and generate string
% EXAMPLE: digit = 5; j = 12;

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

