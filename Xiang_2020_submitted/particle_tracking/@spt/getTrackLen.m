%% Calculate the length (number of steps) of individual tracks.
%
% -About-
%   The method getTrackLen() loops through all tracks and finds out the
%   length of individual tracks.
%
% -Input-
%   - obj: spt object
% 
% -Output-
%   - trackLen: the lengths of all tracks. A Nx1 matrix, where N is the
%     number of tracks.
% 
% -Example-
%   % Calculate the lengths of all tracks
%   myParticle.getTrackLen();
% 
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University


function getTrackLen(obj)
% Memory for track length
obj.trackLen = zeros(length(obj.tracks),1);
% Loop through all tracks
for ii = 1:length(obj.tracks)
    % Extract a single track
    thisTrack = obj.tracks{ii};
    % Find its length and store to the trackLen
    obj.trackLen(ii) = size(thisTrack,1);
end
end