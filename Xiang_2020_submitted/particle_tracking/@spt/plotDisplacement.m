%% Plot the displacement/distance distribution.
%
% -About-
%   The method plotDisplacement() plots the distribution of displacements
%   stored in the object. If the displacement was not calculated
%   previously, it will be calculated by default at the smallest time delay
%   first before the plot can be generated. By default, plotDisplacement()
%   uses displacements from both dimensions (x and y). The user can specify
%   using keywords to chose only to use displacements in a certain
%   dimension for the calculation of the output distribution. Also, the
%   user can choose to plot the distribution using the absolute values of
%   the displacement or/and normalize by standard deviation.
%
% -Input-
%   - obj: spt object
%
% -Varargin-
%   - x:    displacement distribution using displacements only in the x
%           dimension
%   - y:    displacement distribution using displacements only in the y
%           dimension
%   - xy:   displacement distribution using displacements in both x and y
%           dimensions
%   - abs:  if the absolute of the displacement should be used. i.e.
%           distance. Boolean, either true or false (1 or 0)
%   - norm: if the displacement should be normalized by its std. dev..
%           Boolean, either true or false (1 or 0).
%
% -Output-
%   A histogram plot of displacement distribution in the specified
%   dimension. The histogram appeared similar to a scatter plot. Red
%   represents in x dimension, blue represents in y dimension and green
%   represents in collapsed data in dimensions.
%
% -Example-
%   % Plot the displacement distribution by default settings, i.e. use
%     displacements in both x and y dimensions
%   myParticle.plotDisplacement()
%   % Plot the displacement distribution using displacements only in the x
%     dimension
%   myParticle.plotDisplacement('x')
%   % Plot the displacement distribution using displacements only in the y
%     dimension
%   myParticle.plotDisplacement('y')
%   % Plot the displacement distribution using displacements in both x and
%     y dimensions
%   myParticle.plotDisplacement('xy')
%   % Plot the distance distribution
%   myParticle.plotDisplacement('xy','abs',1);
%   % Plot the displacement distribution normalized by std dev.
%   myParticle.plotDisplacement('xy','norm',1);
%
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

function plotDisplacement(obj, varargin)

% Check if the displacement has been calculated before. If not, calculate
% the displacements using the default, i.e. the smallest time delay
if isempty(obj.displacement)
    obj.getDisplacement();
    fprintf('No stored displacement found. Now calculating...\n');
    fprintf('Plotting displacement distribution at %d msec.\n', round(1e3*obj.frameTime));
    ds = obj.displacement;
else
    ds = obj.displacement;
end


% By default, no abs() or normalization are taken on the displacement
absBool = 0;
normBool = 0;

% If no user input, displacements in both x and y dimensions are combined
data = ds;

if ~isempty(varargin)
    for ii = 1:2:length(varargin)
        switch lower(varargin{ii})
            % Use displacement in x dimension only
            case 'x'
                data = ds(:,1);
                % Use displacement in y dimension only
            case 'y'
                data = ds(:,2);
                % Use displacement in both dimensions
            case 'xy'
                % Check if need to take absolute values
            case 'abs'
                if (varargin{ii+1} == 1)
                    absBool = 1;
                end
                % Check if need to normalize the data by std. dev.
            case 'norm'
                if (varargin{ii+1} == 1)
                    normBool = 1;
                end
            otherwise
                error('Invalid argument. Valid options: x|y|xy|abs|norm')
        end
    end
end

% Calculate distance
if absBool
    data = abs(data);
end

% Normalize by std. dev.
if normBool
    data = data ./ std(data);
end

% Calculate the number of bins
bin = floor(sqrt(length(data(:))));
% Calculate the center and height of each bin
[count, center] = hist(data(:),bin);
% Calculate the normalization factor
normFactor = sum(count.*(center(2)-center(1)));
% Normalize the bin
count = count / normFactor;
% Plot the normalized histogram
plot(center,count,'.'); hold on;

% Label the axis properly
% Default
xlabel('Displacement [µm]');

if absBool
    set(gca,'xscale','log');
    xlabel('Distance [µm]');
end

if (normBool && absBool)
    xlabel('Normalized distance');
elseif normBool
    xlabel('Normalized displacement')
end

ylabel('Prob. density')
set(gca,'linewidth',1.5,'fontsize',14,'yscale','log');
axis square;

end