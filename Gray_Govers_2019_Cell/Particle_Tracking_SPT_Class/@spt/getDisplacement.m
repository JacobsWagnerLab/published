%% Calculate the displacement distributions from particle tracks.
%
% -About-
%   The diffusion coefficient is proportional to the time derivative of the
%   variance of the displacement distribution. A greater time delay
%   corresponds to a greater variance in the displacement. We can deduce if
%   any drift was present by checking if the displacement distribution is
%   centered around 0. The shape of the distribution is also important. For
%   a purely Brownian particle, its displacement distribution is always
%   Gaussian regardless of the delay. However, for a sub-diffusive
%   particle, its distribution turns from Gaussian to non-Gaussian at
%   longer delays. This function allows the user to plot the distribution
%   of displacement (default) or distance (absolute values of the
%   displacement) at a specific time delay (1 frame by default).
% 
% -Inputs-
%   - obj: spt object with tracks stored
% 
% -Varargin
%   - delay: a positive integer specifying the delay at which the
%            distribution will be calculated.
% 
% -Output-
%   - displacement: an array of displacement calculated from the tracks
%                   stored in the object. The unit is converted to physical
%                   units.[µm, µm]
% 
% -Example-
%    % Default: displacement distribution at 1-frame delay
%    myParticle.getDisplacement()
%    % Displacement distribution at 10-frame delay
%    myParticle.getDisplacement(10)
%
% -Author-
%    Yingjie Xiang, CJW Lab, Yale University

function getDisplacement(obj,varargin)

% Default: 1-frame delay displacement distribution.
if isempty(varargin)
    delay = 1;
else
    delay = varargin{1};
end

% To be sure, update the track length
obj.getTrackLen();

% Check if all tracks are longer than the delay specified
if (min(obj.trackLen) <= delay)
    error('At least one track is shorter than the specified delay.');
end

% Memeory for displacement in x and y dimensions (left and right columns)
% We have some complicated dimension calculation here.
% Pre-allocation of memory here helps to speed up the program.
displacement = zeros(sum(obj.trackLen-delay)-1,2);
ix = [1; cumsum(obj.trackLen-delay)+1];

% Loop through each track
for ii = 1: length(obj.tracks)
    % Calculate the displacement
    track = obj.tracks{ii};
    d = track(delay+1:end,2:3) - track(1:end-delay,2:3);
    % Append the displacement calculated from this track 
    displacement(ix(ii):ix(ii) + obj.trackLen(ii)-delay-1,:) = d;
end

% Convert the displacement into physical units and pass to the object
obj.displacement = displacement .* obj.pixelLength;

end