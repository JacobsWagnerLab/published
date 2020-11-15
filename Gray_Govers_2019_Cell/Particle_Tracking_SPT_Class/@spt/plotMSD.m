%% Plot individual MSD curves.
%
% -About-
%   The method plotMSD() plots the individual MSD curves in a loglog scale.
%   The user can specify the indices of the MSDs to be plotted and the plot
%   styles using keywords.
% 
% -Input-
%   - obj: spt object
% 
% -Varargin-
%   - index: an index array used to specify which MSD curves to be plotted
%   - color: color used for the curves
% 
% -Output-
%   A plot of individual MSD curves in loglog scales.
%
% -Example-
%   % Plot all the MSD curves
%   myParticle.plotMSD();
%   % Plot the first 30 MSD curves
%   myParticle.plotMSD('index',1:30);
% 
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

function plotMSD(obj,varargin)

% Check if the MSD are calculated beforehand. If not, calculate now.
if isempty(obj.MSD)
    obj.getMSD();
end

% Default parameters
% Indicies of the MSDs to be plotted
index = 1:size(obj.MSD,1);
% No color specified by default
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

% If no color is specified and multiple MSD curves were to be plotted, use
% a default jet colormap
if (isempty(color)) && (length(index) ~= 1)
    color = jet(length(index));
else
    color = obj.selectColor(color);
end
    
% Go through the indexed MSD and plot each of them
for ii = 1:length(index)
    % Convert to physical units
    delay = obj.MSD{index(ii)}(:,1) .* obj.frameTime;
    msd = obj.MSD{index(ii)}(:,2) .* (obj.pixelLength)^2;
 
    % Plot the MSD curves
    if size(color,1) == 1
        plot(delay,msd,'color',[color,0.3],'linewidth',1.5);
    else 
        plot(delay,msd,'color',[color(ii,:),0.3],'linewidth',1.5);
    end
    hold on;
end

axis square;
grid on;
xlabel('Time [sec]');
ylabel('MSD [µm^2]');
set(gca,'xscale','log','yscale','log','linewidth',1.5,'fontsize',14);
end
