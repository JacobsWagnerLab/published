function img = loadimseries(pathToFolder,varargin)
%--------------------------------------------------------------------------
%function img = loadimseries(pathToFolder,varargin)
%@author:  Manuel Campos
%@date:    February, 2014
%@copyright 2013-2014 Yale University
%==========================================================================
%**********output********:
%img:       3D array of class unit16 containing the image data as 
%           nb line x nb column x nb frames
%           
%**********input********:
%pathToFolder: String specifying the path to the folder where images are
%              stored
%Non-mandatory input arguments
% 1-        Vector of length 2 specifying the staring and ending file
%           number on the list of *.tif files that are to be loaded
% 2-        Step size between the *.tif files to be loaded. 
%           i.e. set this argument to 2 if you want to load one *.tif file
%           out of two 
%           
%==========================================================================
% This function bypasses the constrains associated with loadimageseries.m
% that is formatted to fit the internal requirements of Oufti.
% It adds the possibility to load a specific subset of the images stored in
% the desired folder.
%--------------------------------------------------------------------------

if isdir(pathToFolder)
    % Select the file names with a *.tif extension
    files = dir([pathToFolder '/*.tif*']);
    % Set to boundaries of the files to be loaded to the maximum
    lb = 1;
    ub = length(files);
    loopArray = lb:ub;
    if nargin==2
        % if a first non-mandatory input is provided, get the values of the
        % first and last *.tif file to be loaded (lower and upper
        % boundarie)
        fbounds = varargin{1};
        lb = fbounds(1);
        ub = fbounds(2);
        loopArray=lb:1:ub;
    elseif nargin==3
        % if a second non-mandatory input is provided, get the value of the
        % step
        loopArray = lb:varargin{2}:ub;
    end
    % Load the first file and use it to allocate the memory for the output
    % img variable
    firstFile = imread([pathToFolder '/' files(1).name]);
    img = uint16(zeros(size(firstFile,1),size(firstFile,2),length(loopArray)));
    % Loop through the list o f*.tif file names to load the requested
    % images
    w = waitbar(0, 'Loading image files, please wait...');
    for jj = 1:length(loopArray)
        ii = loopArray(jj);
        img(:,:,jj) = imread([pathToFolder '/' files(ii).name]);
        waitbar(jj/length(loopArray));
    end
end
close(w)
end