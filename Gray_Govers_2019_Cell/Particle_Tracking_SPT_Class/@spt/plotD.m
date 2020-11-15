%% Plot the diffusion coefficients using the kernel density estimator.
%
% -About-
%   The method plotD() uses the kernel density estimator to plot the
%   distribution of the diffusion coefficients stored in the object. A
%   Gaussian kernel is used, and the bandwidth of the kernel is optimized
%   by the MATLAB built-in functions ksdensity() automatically. The user
%   can change the line width as well as the line color in the plot by
%   specifying with keywords. The default is 1.5 for line width and black
%   for line color.
%
% -Input-
%   - obj: spt object
%
% -Varargin-
%   - linewidth: line width
%   - color:     line color
%   - log:       Boolean. True or false (1 or 0): whether to log-transform
%
% -Output-
%   - A plot of probability density distribution of diffusion coefficients.
%
% -Example-
%   % Plot the probability density distribution of diffusion coefficients
%     using default line styles
%   myParticle.plotD();
%   % Plot the probability density distribution of diffusion coefficients
%     using red color and 2.0 line width
%   myParticle.plotD('color','red','linewidth',2);
%
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

function plotD(obj,varargin)

% Check if diffusion coefficients have been calculated. If not, calculate now.
if isempty(obj.D)
    fprintf('No diffusion coefficients found. Now calculating (first 5 delays)...\n');
    obj.fitMSD(1:5,'linear');
end

lw = 2;
c = '';
logX = 0;

% Parse user inputs and overwrite the defaults
if ~isempty(varargin)
    for ii = 1:2:length(varargin)
        switch varargin{ii}
            case 'linewidth'
                lw = varargin{ii+1};
            case 'color'
                c = obj.selectColor(varargin{ii+1});
            case 'log'
                logX = varargin{ii+1};
            otherwise
                error('Invalid options. Valid options: linewidth|color|log');
        end
    end
end

% Extract all diffusion coefficients
x = obj.D;
if logX
    x = log10(obj.D(obj.D>0));
end

% Calculate the kernel density
[fi, xi] = ksdensity(x);

if ~ isempty(c)
    plot(xi,fi,'linewidth',lw,'color',c);
else
    plot(xi,fi,'linewidth',lw);   
end

axis square;
set(gca,'fontsize',14,'linewidth',1.5);
xlabel('Diffusion coefficients [µm^2/s]');
ylabel('Prob. density');

end