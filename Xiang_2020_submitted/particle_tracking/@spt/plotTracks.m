%% Plot individual tracks.
%
% -About-
%   The method plotTracks() plots the particle trajectories stored in the
%   object. The user can provide an index array to specify which tracks to
%   be plotted.
% 
% -Input-
%   - obj: spt object
% 
% -Varargin-
%   - index: An index array specifying which tracks to be plotted. If no
%            index provided, all tracks stored in the object are plotted
%   - color: color for plotting the tracks
%
% -Output-
%   A plot of particle trajectories
% 
% -Example-
%   % Plot all tracks
%   myParticle.plotTracks();
%   % Plot the 3rd track
%   myParticle.plotTracks(3);
%   % Plot the first 10 tracks
%   myParticle.plotTracks(1:10);
% 
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University


function plotTracks(obj,varargin)

% Check if an index array for which tracks to be plotted are given
% If no index is given, plot all tracks.
index = 1:length(obj.tracks);
color = '';

% Overwrite the default parameters based on the user input
for ii = 1:2:length(varargin)
    switch lower(varargin{ii})
        case 'color'
            color = varargin{ii+1};
        case 'index'
            index = varargin{ii+1};
        otherwise
            error('Invalid argument. Valid options: index|color.');
    end
end

if (isempty(color)) && (length(index) ~= 1)
    color = jet(length(index));
else
    color = obj.selectColor(color);
end

% Plot all tracks
for ii = 1: length(index)
    
    % Extract a single track
    thisTrack = obj.tracks{index(ii)};
    
    % Convert to physical units
    x = thisTrack(:,2).*obj.pixelLength;
    y = thisTrack(:,3).*obj.pixelLength;
    if size(color,1) == 1
        plot(x,y,'color',color);
    else
        plot(x,y,'color',color(ii,:));
    end
    hold on;
end

set(gca,'fontsize',14,'linewidth',1.5);
xlabel('x [µm]');
ylabel('y [µm]');

axis equal;
pbaspect([1 1 1]);
end
