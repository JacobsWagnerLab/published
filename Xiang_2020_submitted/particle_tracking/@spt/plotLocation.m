%% Plot the localized particle positions with the corresponding frames.
%
% -About-
%   After the step of particle localization, the locations of all particles
%   identified are stored in the locs. In order to check the accuracy of
%   the localization, the original frame was plotted first and the particle
%   positions localized are plotted on top. This enables the user to make
%   visual inspection on the correctness as well as the performance of the
%   localization.
%
% -Input-
%   - obj: spt object
%
% -Output-
%   - Displays and loops through all original frames with localized
%   particle positions plotted on top. The user must press any key to
%   continue to the next frame. The plot title indicates the number of
%   total frames as well as the current frame.
%
% -Example-
%   % Check the particle positions localized by getLocation()
%   myParticle.plotLocation();
%
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

function plotLocation(obj,varargin)

if nargin < 2
    % Maximum frame number
    maxFrameNum = max(obj.locs(:,3));
    % Loop through each frame
    for ii = 1 : maxFrameNum
        % Collect all particle positions that appeared in this frame
        frame = obj.locs(:,3);
        index = find(frame == ii);
        % If no particle appeared in this frame, skip the frame
        if isempty(index), continue, end
        % Plot the particle position on top of the frame
        figure(1); clf;
        im = double(imread(obj.stackPath, 'Index', ii));
        imshow(im,[min(im(:)), max(im(:))]);
        hold on;
        plot(obj.locs(index,1), obj.locs(index,2),'gx','markersize',10);
        colormap('gray'); colorbar;
        title(['Frame ', num2str(ii)]);
        pause();
    end
    close all;
else
    target = varargin{1};
    frame = obj.locs(:,3);
    index = find(frame == target);
    % Plot the particle position on top of the frame
    figure(1); clf;
    im = double(imread(obj.stackPath, 'Index', target));
    imshow(im,[min(im(:)), max(im(:))]);
    hold on;
    plot(obj.locs(index,1), obj.locs(index,2),'gx','markersize',10);
    colormap('gray'); colorbar;
    title(['Frame ', num2str(target)]);
end
end