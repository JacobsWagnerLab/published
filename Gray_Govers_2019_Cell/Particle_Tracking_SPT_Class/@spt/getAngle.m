%% Calculate the angles between consecutive displacements.
%
% -About- 
%   The angles between consecutive displacements of a particle can indicate
%   the presence of the drift. For a Brownian particle, the displacement is
%   completely uncorrelated. The angles between displacement are therefore
%   uniformly distributed at 1/pi, with a mean of pi/2. If there was a
%   drift, the likelihood of having a smaller angle would be increased from
%   this base-line.
% 
% -Inputs- 
%   - obj: spt object with tracks stored
% 
% -Output- 
%   - angles:       angles between consecutive displacements
%   - meanAngles:   mean of angles per track
%
% -Example-
%   % Calculate the angles between consecutive displacements
%   myParticle.getAngle();
%
% -Author- 
%   Yingjie Xiang, CJW Lab, Yale University

function getAngle(obj)

% Preallocation for angles
obj.getTrackLen()
angles = zeros(sum(obj.trackLen-2)-1,1);
meanAngles = zeros(length(obj.tracks),1);
ix = [1; cumsum(obj.trackLen-2)+1];

% Loop through every track
for ii = 1:length(obj.tracks)
    thisTrack = obj.tracks{ii};
    % Extract the x and y coordinates
    x = thisTrack(:,2);
    y = thisTrack(:,3);
    % Calculate the displacements
    dx = x(2:end) - x(1:end-1);
    dy = y(2:end) - y(1:end-1);
    % Memory for storing angles
    angle = zeros(length(dx)-1,1);
    
    % Loop through the step
    for jj = 1:length(dx)-1
        % Preparing the two vectors representing the displacements
        a = [dx(jj),dy(jj)];
        b = [dx(jj+1),dy(jj+1)];
        % Normalize the vectors
        a = a ./ norm(a);
        b = b ./ norm(b);
        % Calculate the angles between the vectors
        angle(jj) = acos(dot(a,b));
    end
    
    angles(ix(ii):ix(ii) + obj.trackLen(ii)-3) = angle;
    meanAngles(ii) = mean(angle);
end

obj.angles = angles;
obj.meanAngles = meanAngles;

end