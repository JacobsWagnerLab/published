%% Link the particle locations to create tracks.
%
% -About-
%   The method getTrack() invokes the dependency function track() and
%   converts its output into a cell array. This data structure conversion
%   makes the analysis and calculations later easier. For full details
%   about the track() function, check its comments and the reference below.
% 
% -Input-
%   - obj: spt object
%   - maxDist: maximum distance traveled by a particle within 1 frame [px]
%   - minTrackLen: minimum length below which a track will be discarded
% 
% -Output-
%   - tracks: the tracks of individual particles. A Nx1 cell array, with
%             each cell representing a particle. Within each cell, there is
%             a Nx3 matrix, where N is the length of each track (number of
%             steps). The three columns (left to right) represents the
%             frame number, positions in x and y dimensions. The units are
%             not converted to physical units. [frame px px]
% 
% -References-
%   Crocker JC, Grier DG, Journal of Colloid Interface Science. 179 298-310
%   (1996).
% 
% -Example-
%   % Link particles that can travel no more than 5 pixels within 1 frame,
%     and eliminate any tracks shorter than 10 steps.
%   maxDist = 5;
%   minTrackLen = 10;
%   myParticle.getTrack(maxDist,minTrackLen);
% 
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

function getTrack(obj,maxDist,minTrackLen)

% Time of disappearance of a paritcle when it is forgotten
param.mem = 0; 
% Dimension of tracking
param.dim = 2;
% Whether to print warning
param.quiet = 1;
% Remove tracks that have fewer steps than defined value
param.good = minTrackLen; 

% Original output from the track function
% A Nx4 matrix containing all particle positions and their associated
% particle ID. N is the number of particle positions, the columns from left
% to right represents the particle positions in x and y dimensions, the
% frame at which it was located, the particle ID it has. That is, tks is in
% the format of [x y t id]
tks = track(obj.locs,maxDist,param); 

% Transform of the data structure
% Convert the Nx4 matrix into a Mx1 cell array, where M is the number of
% particles. Within each cell, there is a Kx3 matrix, where K is the number
% of positions. The three columns (left to right) represents the frame
% number, positions in x and y dimensions.

% [x y t id] --> [t x y] by particle id
ids = unique(tks(:,4));
tracks = cell(length(ids),1);
for ii = 1 : length(ids)
    tracks{ii} = tks(tks(:,4)== ids(ii),[3,1,2]); % t x y
end

% Pass the results back to the object
obj.tracks = tracks(~cellfun(@isempty,tracks));
end