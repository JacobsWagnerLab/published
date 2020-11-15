%% Filter the tracks based on specified crtieron, remove those outside the range.
%
% -About-
%   Filter out the trajectories based on the specified criterion and range.
%   Only the trajectories inside the range will be kept. Those outside the
%   range will be deleted and the MSD and meanMSD are re-calculated.
%
% -Inputs-
%   - obj:      spt object
%   - criterion char array. Either 'trackLen' or 'intercept' or 'roi'
%   - range     range within which the tracks will be kept
%   - reverse   boolean, if true, tracks outside the range are kept
%
% -Output-
%   - filtered tracks
%
% -Example-
%   % Filter out the tracks below 10 steps and longer than 50 steps
%   myParticle.filterTrack('trackLen',[10,50]);
%   % Filter out the tracks with y-intercepts larger than 1e-2;
%   myParticle.filterTrack('intercept',[-Inf,1e-2]);
%
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

function filterTrack(obj,criterion,range,reverse)

if nargin < 3
    error('Need to specify both the filtering criterion and range to keep');
end

if range(1) > range(2)
    error('Filtering range must be in a numeric ascending order');
end

if nargin < 4
    reverse = 0;
end

switch lower(criterion)
    case 'tracklen'
        obj.getTrackLen();
        if reverse
            obj.tracks(obj.trackLen > range(1) & obj.trackLen < range(2)) = [];
        else
            obj.tracks(obj.trackLen < range(1) | obj.trackLen > range(2)) = [];
        end
    case 'intercept'
        obj.getIntercept();
        if reverse
            obj.tracks(obj.intercepts > range(1) & obj.intercepts < range(2)) = [];
        else
            obj.tracks(obj.intercepts < range(1) | obj.intercepts > range(2)) = [];
        end
    case 'roi'
        mins = cellfun(@min, obj.tracks,'UniformOutput',0);
        mins = vertcat(mins{:});
        mins = mins(:,2:3).*obj.pixelLength;
        maxs = cellfun(@max, obj.tracks,'UniformOutput',0);
        maxs = vertcat(maxs{:});
        maxs = maxs(:,2:3).*obj.pixelLength;
        minx = range(1); maxx = range(2); miny = range(3); maxy = range(4);
        within_roi = mins(:,1) > minx & mins(:,2) > miny & ...
                     maxs(:,1) < maxx & maxs(:,2) < maxy;
        if reverse
            obj.tracks(within_roi) = [];
        else
            obj.tracks(~within_roi) = [];
        end
        
    otherwise
        fprintf('Invalid criterion: use either "trackLen" or "intercept" or "roi" followed by the range to keep \n');
end

obj.getTrackLen();
obj.getIntercept();
obj.getMSD();
obj.getMeanMSD();
obj.D = [];
obj.alpha = [];
obj.ensembleD = [];
obj.ensembleAlpha = [];

end