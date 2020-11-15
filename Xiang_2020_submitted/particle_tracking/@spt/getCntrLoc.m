%% Calculate the angles between consecutive displacements.
%
% -About- 
%   The method getCntrLocs() first centers all tracks stored in the object
%   to the origin (0,0), by substracting the entire track from the initial
%   position. Then, it collapses all particle locations into a single
%   matrix with 2 columns. These two columns are particle positions in x
%   any y dimensions. This summarizes all the particle locations and can be
%   used to construct a density plot to investigate the length scale of the
%   space explored by all particles given the corresponding track length.
% 
% -Inputs- 
%   - obj: spt object with tracks stored
% 
% -Output- 
%   - cntrLocs: collpased locations from the origin-centered tracks. The
%               unit is converted to physical unit. [µm]
%
% -Example-
%   % Collapse all locations from the origin-centered tracks
%   myParticle.getCntrLoc()
%
% -Author- 
%   Yingjie Xiang, CJW Lab, Yale University

function getCntrLoc(obj)

% Calculate the track length for pre-allocation of memory
obj.getTrackLen();
ix = cumsum([1; obj.trackLen]);
locs = zeros(sum(obj.trackLen),2);

% Collapse the positions from all origin-centered tracks
for ii = 1:length(obj.tracks)
    % Extact a single track
    thisTrack = obj.tracks{ii};
    % Substract the first position to center to the origin
    thisTrack = thisTrack - thisTrack(1,:);
    % Store the centered location based on the index
    locs(ix(ii):ix(ii)+obj.trackLen(ii)-1,:) =  thisTrack(:,2:3);
end

% Remove all the origin points
locs = locs(any(locs,2),:);

% Convert to physical units and pass the results back to the object
obj.cntrLocs = locs .* obj.pixelLength;

end