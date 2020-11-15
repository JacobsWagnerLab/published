%% Plot the alpha distribution using the kernel density estimator.
%
% -About-
%   The method plotAlpha() uses the kernel density estimator to plot the
%   distribution of the alpha values stored in the object. A Gaussian
%   kernel is used, and the bandwidth of the kernel is optimized by the
%   MATLAB built-in functions ksdensity() automatically. The user can
%   change the line width as well as the line color in the plot by
%   specifying with keywords. The default is 1.5 for line width and black
%   for line color.
% 
% -Input-
%   - obj: spt object
% 
% -Varargin-
%   - linewidth: line width
%   - color: line color
% 
% -Output-
%   - A plot of probability density distribution of alpha values.
%
% -Example-
%   % Plot the probability density distribution of alpha values using
%     default line styles
%   myParticle.plotAlpha();
%   % Plot the probability density distribution of alpha values using red
%     color and 2.0 line width
%   myParticle.plotAlpha('color','red','linewidth',2);
%
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University


function plotAlpha(obj,varargin)

% Check if alpha values have been calculated. If not, calculate now.
if isempty(obj.alpha)
    fprintf('No alpha values found. Now calculating (first 5 delays)...\n');
    obj.fitMSD(2:6,'loglog');
end

lw = 2;
c = '';

% Parse user inputs and overwrite the defaults
if ~isempty(varargin)
    for ii = 1: 2: length(varargin)
        switch varargin{ii}
            case 'linewidth'
                lw = varargin{ii+1};
            case 'color'
                c = obj.selectColor(varargin{ii+1});
            otherwise
                error('Invalid options. Valid options: linewidth|color');
        end
    end
end

% Calculate the kernel density
[fi, xi] = ksdensity(obj.alpha);
if ~ isempty(c)
    plot(xi,fi,'linewidth',lw,'color',c);
else
    plot(xi,fi,'linewidth',lw);   
end

axis square;
set(gca,'fontsize',14,'linewidth',1.5);
xlabel('Alpha');
ylabel('Prob. density');

end