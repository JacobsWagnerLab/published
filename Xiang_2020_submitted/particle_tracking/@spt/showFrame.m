%% Show specific frames from the image stack.
%
% -About-
%   Show specific frames in the image stack associated with the object. By
%   default (no input) it shows all the frames.
%
% -Input-
%   obj: spt object
%
% -Varargin-
%   the indices of the frames to be showed
%
% -Output-
%   the specified frames are showed
%
% -Example-
%   %Show all frames in the stack
%   myParticle.showFrame();
%   %Show the first 10 frames in the stack
%   myParticle.showFrame(1:10);
%   %Show the 20th frame in the stack
%   myParticle.showFrame(20);
%
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

function showFrame(obj,varargin)

% Check if an index array is provided. If not, show the 1st frame only.
if isempty(varargin)
    index = 1 : obj.imNum;
else
    index = varargin{1};
end

if isempty(obj.stackPath)
    error('No image stack path found.');
end

% Show the indexed frame
for ii = 1 : length(index)
    figure(1); cla;
    im = double(imread(obj.stackPath, 'Index', index(ii)));
    % The display range is set to be between the lowest and highest
    % intensities
    imshow(im,[min(im(:)), max(im(:))]);
    title(['Frame #',num2str(index(ii))],'fontsize',14);
    colorbar;
    pause();
end