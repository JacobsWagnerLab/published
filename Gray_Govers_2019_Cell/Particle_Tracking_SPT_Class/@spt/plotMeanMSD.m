%% Plot the ensemble-averaged MSD.
%
% -About-
%   The method plotMeanMSD() plots the ensemble-averaged MSD in a loglog
%   scale. The slope of the loglog MSD plot indicates the mode of the
%   diffusion (e.g. normal, sub-diffusive, super-diffusive). The
%   y-intercept also scales with the population-averaged diffusion
%   coefficient.
% 
% -Input-
%   - obj: spt object
% 
% -Output-
%   A plot of the ensemble-averaged mean squared displacement in loglog
%   scales.
% 
% -Example-
%   % Plot the ensemble-averaged mean squared displacements
%   myParticle.plotMeanMSD();
% 
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University


function plotMeanMSD(obj)

% Check if the ensemble-averaged MSD was calculated beforehand. If not
% calculate now.
if isempty(obj.meanMSD)
    obj.getMeanMSD();
end

% Plot the ensemble-averaged MSD using physical units
plot(obj.meanMSD(:,1).*obj.frameTime,...
     obj.meanMSD(:,2).*(obj.pixelLength)^2,'s','markersize',10); 

axis square;
grid on;
xlabel('Time [sec]');
ylabel('<MSD> [µm^2]');
set(gca,'xscale','log','yscale','log','linewidth',1.5,'fontsize',14);
end