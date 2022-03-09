
function [ret1] = CallStack_debackV3B(ScriptDir, MainDir, RefPath, DataNamePre, ParamA, ParamB, ChNum)

% =============================================================== %

ExpName1 = 'c1Exp'; 
ExpName2 = 'c2Exp';
ExpName3 = 'c3Exp';

cd(MainDir);  mkdir(ExpName1); mkdir(ExpName2); mkdir(ExpName3);

ImpDir1 = strcat(MainDir, 'c1\');
ImpDir2 = strcat(MainDir, 'c2\');
ImpDir3 = strcat(MainDir, 'c3\');

ExpDir1 = strcat(MainDir, ExpName1);
ExpDir2 = strcat(MainDir, ExpName2);
ExpDir3 = strcat(MainDir, ExpName3);

% =============================================================== %

%StepFlag = 1;
SampleInterval = ParamA(1:3); %[1 3 3];
MaxDigit = ParamA(4);      

% Range of background image crop 
RgR = ParamB(1:2);  %[200 20];  % up, down directions
RgC = ParamB(3:4);  %[20 20];   % left, right directions

% (1) Get the initial data mark (which will be moved and optimized)
DataMarkIni = ParamB(5:6);  %[512 740];  % Use (y,x) in ImageJ
% Note: (Vertical, Horizontal) in ImageJ = (Row, Col) in Matlab

% (2) Get fata initial marker
StackRg = ParamB(7:8);        % [1 800];
RefMark = ParamB(9:10);       %[481 40];  % Use (y,x) in ImageJ
LowerLimit = ParamB(11);      %-4000;

RefImg = double( imread( RefPath ) );

% =============================================================== %

% Find DataMarkOpt by local perturbation and minimize the image gradient.

DataMarkOpt = DataMarkIni;

n = StackRg(2) - StackRg(1) + 1;
OptMarkRecRD2 = zeros(n, 2);

L = 4;  
% The +/- range in X,Y (in pixel) that will search for optimal deback
% For example, in L=4 case, the script will test all conditions with 
% x = (-4) to +4 and y = (-4) to 4 and find optimal deback among these conditions.

rem_flag = 1;

for j = StackRg(1) : StackRg(2)
    
    cd(ScriptDir);
    ImgNameJ = strcat(DataNamePre, AddZeros(MaxDigit, j), 'c1', '.tif') ; 
    
    cd (ImpDir1);
    DataJ = double( imread (ImgNameJ) );
    
    cd(ScriptDir);
    
    [OptDebackData, DataMarkOptJ, Simple_Deback] = ...
        Optimal_DebackV2(DataJ, RefImg, DataMarkOpt, RefMark, L, RgR, RgC);            
    
    DataMarkOpt = DataMarkOptJ;
    
    %
    
    rec_ind = j - StackRg(1) + 1;
    OptMarkRecRD2(rec_ind, :) = DataMarkOptJ;    
    
    
    % (1) Phase image
    
    OptDebackData = OptDebackData + LowerLimit;     
    OptDebackData = OptDebackData.*(OptDebackData>0);   % remove the value below the lower limit

    cd(ScriptDir);
    ExportNameJ2 = strcat(DataNamePre, AddZeros(MaxDigit, j), 'c1DB', '.tif');    
    
    cd (ExpDir1);
    imwrite( uint16(OptDebackData), ExportNameJ2, 'tif');    
    
    
    % (2) 405 image
    
    if ( ( rem(j, SampleInterval(2)) == rem_flag  )||(SampleInterval(2) == 1) )
        
        if (ChNum >= 2) 
        
        cd (ImpDir2);    
        ImgNameJ2 = strcat(DataNamePre, AddZeros(MaxDigit, j), 'c2', '.tif') ;         
        DataJ2 = double( imread (ImgNameJ2) );
                 
        DataRegionJ2 = DataJ2( (DataMarkOpt(1)-RgR(1)) : (DataMarkOpt(1)+RgR(2)), ...
                               (DataMarkOpt(2)-RgC(1)) : (DataMarkOpt(2)+RgC(2)) );
                       
          
        cd(ScriptDir);
        ExportNameJ2 = strcat(DataNamePre, AddZeros(MaxDigit, j), 'c2s', '.tif');    
        
        cd(ExpDir2); 
        imwrite( uint16(DataRegionJ2), ExportNameJ2, 'tif'); 
                      
        end
        
    end
    
    % (3) 488 image
       
    if ( ( rem(j, SampleInterval(3)) == rem_flag  )||(SampleInterval(3) == 1) )
        
        if (ChNum >= 3) 
            
        cd (ImpDir3);    
        ImgNameJ3 = strcat(DataNamePre, AddZeros(MaxDigit, j), 'c3', '.tif') ;         
        DataJ3 = double( imread (ImgNameJ3) );
                 
        DataRegionJ3 = DataJ3( (DataMarkOpt(1)-RgR(1)) : (DataMarkOpt(1)+RgR(2)), ...
                               (DataMarkOpt(2)-RgC(1)) : (DataMarkOpt(2)+RgC(2)) );
                       
           
        cd(ScriptDir);
        ExportNameJ3 = strcat(DataNamePre, AddZeros(MaxDigit, j), 'c3s', '.tif');    
        
        cd(ExpDir3);
        imwrite( uint16(DataRegionJ3), ExportNameJ3, 'tif'); 
        
        end
                       
    end
    
     
    
end

cd(MainDir);  save('OptMarkRecRD2.mat', 'OptMarkRecRD2');

ret1 = 1;

end


% =======================================================================

function [StrOutput] = AddZeros(digit, j)
% Adding zero before number and generate string

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


