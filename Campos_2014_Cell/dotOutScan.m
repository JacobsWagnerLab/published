function input = dotOutScan(pathFile)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function input = dotOutScan(pathFile)
%@author:  Ahmad J Paintdakhi
%@date:    September 12, 2013
%@copyright 2013-2014 Yale University
%==========================================================================
%**********output********:
%input:     Structure with cell informations organized into fields that are
%           sorted similarly
%**********input********:
%pathFile:  Path to the *.out file containing the cell informations from
%           the high-throughput experiment of interest.
%==========================================================================
%The function selectivey transform the data from the csv type *.out file
%into a Matlab structure that can be used to sift through the cells and
%construct a database of well-detected cells with autoSelect.m (v4 to date)
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 

%% 1.  Find the byte location where the header information is stored in the .out file.
     
% fname = '\\128.36.22.176\Users 1\Manuel Campos\FM\Cell Size\2014-02-28 ssb_microF\tl1.out';
fname = pathFile;
tempText = memmapfile(fname);
tempChar = char(tempText.Data(1:100000))';
frameNumberLocation = strfind(tempChar,'cellId,');
clear tempChar tempText

%% 2.  Open the file and read the first 1-million bytes after skipping comments and headerline.
      
fid = fopen(fname);

fseek(fid,frameNumberLocation,'bof');
%               frameNumber,ancestors,birthframe,box,descendants,divisions,mesh,length,area,polarity,signal0,cellId,
tic;
data = textscan(fid,'%s         %s       %f32    %*s    %s           %s    %s   %*f32 %*f32 %f32     %*s     %f32',...
    3200000,'delimiter',',','HeaderLines',0);
toc;
% The stars indicate the fields (or parameters) that we do not desire to
% extract


% Clean out empty meshes
meshLength=cellfun(@length,data{1,6});
ixZM = meshLength>1;
for i=1:8 
    data{1,i}=data{1,i}(ixZM);
end

     
% Notice that only cellIdArray is in single format, whereas the frameArray
% and meshArray are in string format.
% Also, the frameArray has '#' in front of every character.
% To remove the hash-key, you can use 
% This cleaning can also be done later after you clean up your data.
data{1,1}  = regexprep(data{1,1},'#','');  

% Combine into a structure array for further usage (see autoSlectv4.m)
input.frameArray=data{1,1};
input.ancArray = data{1,2};
input.birthArray = data{1,3};
input.progArray = data{1,4};
input.divArray = data{1,5};
input.meshArray = data{1,6};
input.polarityArray = data{1,7};
input.cellIdArray = data{1,8};

end
% To convert a mesh vector from string to double I would use str2num. 
% For example, mesh1 = str2num(meshArray(1));

%% 3.  To test a given cell
% Find all instances in which a given cell is present:


%% 4.  To load another 1-million bytes you would step back a few bytes to create the overlap.  To do that you simply do

% currentByteLocation = ftell(fid);

% I believe to move back 10 frames or so it is about 25000 bytes.
% fseek(fid,frameNumberLocation+1e9,'bof');   %frewind(fid) will go to position 1 and we do not want that in this case.

% After we move back a given # of bytes then you can apply this line of code again.
% tic;
% data2 = textscan(fid,'%s %s %*f32 %*s %s %*s %s %*f32 %*f32 %*f32 %*s %f32',1000000,'delimiter',',','HeaderLines',1);
% toc;
