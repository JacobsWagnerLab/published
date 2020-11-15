%% Approximate the particle intensity as the volume below fitting Gaussian
%
% -About-
%   getIntensity() estimates the underlying particle intensity stored in
%   the object based on the fitting results (bivariate Gaussian). The
%   intensity was estimated as the volume below the Gaussian, that is, v =
%   2 * sigmaX * sigmaY * A.
%
% -Inputs-
%   - obj:      spt object
%
% -Output-
%   - Particle intensity in arbitrary units
%
% -Example-
%   % Get all the particle intensity
%   myParticle.getIntensity();
%
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

function getIntensity(obj)

% Multiple Gaussian fitting results stored in 1 frame
frames = length(obj.gaussians);
obj.intensity = zeros(size(obj.locs,1),1);
cumIx = 0;

% Go through each frame
for ii = 1:frames
    % Extract the sigmaX and sigmaY
    sigmas = cellfun(@(x)[x.sigmaX,x.sigmaY],obj.gaussians{ii},'UniformOutput',0);
    
    % No particle detected in that frame
    if isempty(sigmas)
        continue;
    end
    
    % Extract the amplitude/height of Gaussians
    As = cellfun(@(x)x.A,obj.gaussians{ii},'UniformOutput',0);
    As = unique(vertcat(As{:}),'stable');
    sigmas = vertcat(sigmas{:});
    sigmas = unique(sigmas,'rows','stable');
    
    obj.intensity(cumIx+1:cumIx+size(sigmas,1)) = 2*pi.*As.*sigmas(:,1).*sigmas(:,2);
    cumIx = cumIx+size(sigmas,1);
end

end