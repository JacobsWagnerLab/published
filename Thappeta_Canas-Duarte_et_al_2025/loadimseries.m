function img = loadimseries(pathToFolder,varargin)


if isdir(pathToFolder)
    files = dir([pathToFolder '/*.tif*']);
    if nargin==2
        fbounds=varargin{1};lb=fbounds(1);ub=fbounds(2);
        loopArray=lb:1:ub;
    elseif nargin==3
        fbounds=varargin{1};lb=fbounds(1);ub=fbounds(2);
        loopArray=lb:varargin{2}:ub;
    else
        lb = 1; ub = length(files);
        loopArray=lb:1:ub;
    end
    firstFile=imread([pathToFolder '/' files(1).name]);
    img=uint16(zeros(size(firstFile,1),size(firstFile,2),length(loopArray)));
    w = waitbar(0, 'Loading image files, please wait...');
    for jj=1:length(loopArray)
        ii=loopArray(jj);
        img(:,:,jj) = imread([pathToFolder '/' files(ii).name]);
        waitbar(jj/length(loopArray));
    end
end
close(w)
end