function [traj] = pos2trajectories(pos, maxdXY, gapFrames, minTrajectoryLength)

%{
-About-
Perform trivial linking of localizations into trajectories by searching for
nearest link in the subsequent frames, tolerant of temperal gaps in
localizations.

-Inputs-
pos: a structured cell array: 
  pos{cellN}.x = a matrix with number of columns == number of frames and
  each row is a localization event. It is OK to have missing
  localizations (which should be indicated with a 0. 

maxdXY is the maximum distance a point can travel between observations


gapFrames is the number of frames allowed to be absent to link 2
trajectories

minTrajectoryLength is the shortest allowed duration of a trajectory

-varargin-
na

-Outputs-

traj: A nested cell array containing trajectories. traj{1} contains
trajectories from cell #1 and traj{1}{1} is the 1st trajectory identified
in cell #1. Each trajectory is n-rows and columns:
1) frame #
2) image x coordinate
3) image y coordinate

-Example-
   
-Supplementary-

-Dependencies-

-References-
n/a

-Author-
Brad Parry, 2014 July 28
%}



traj = cell(1,length(pos));
for cellN = 1:length(pos)
    if isempty(pos{cellN}), continue, end
    if size(pos{cellN}.x,1) == 0, continue, end
    % allXpos and allYpos must have observations at the same locations
    % throughout. Check. if not,
    if sum(~((pos{cellN}.x(:) ~= 0) == (pos{cellN}.y(:) ~= 0))) ~= 0
    %     then fail
        error('X and Y suck')
    end

    for currentFrameINDEX = 1:size(pos{cellN}.x,2)
        %if this is the first iteration for pos01 or no trajectories have been
        %assigned yet, do some initialization
        if (currentFrameINDEX == 1) || isempty(traj{cellN})
            %Deal out positions to initialize trajectories
            ix = (pos{cellN}.x(:,currentFrameINDEX) ~= 0).*(1:size(pos{cellN}.x,1))';
            ix(ix == 0) = [];
            for t = 1:length(ix)
                traj{cellN}{end+1} = [currentFrameINDEX, pos{cellN}.x(ix(t),currentFrameINDEX), pos{cellN}.y(ix(t),currentFrameINDEX)];
            end
            continue
        end

        %compare to previously observed locations within gapFrames
        %get current unlinked observations
        ix = (pos{cellN}.x(:,currentFrameINDEX) ~= 0) .* (1:size(pos{cellN}.x,1))';
        if isempty(ix), continue, end
        ix(ix == 0) = [];
        currX = pos{cellN}.x(ix, currentFrameINDEX);
        currY = pos{cellN}.y(ix, currentFrameINDEX);
        %get previous linked obervations
        prevX = cellfun(@(x) x(end,2),traj{cellN});
        prevY = cellfun(@(x) x(end,3),traj{cellN});

        dx = repmat(prevX(:), [1, length(currX)]) - repmat(currX(:)', [length(prevX), 1]);
        dy = repmat(prevY(:), [1, length(currY)]) - repmat(currY(:)', [length(prevY), 1]);
        dxy = (dy.^2 + dx.^2).^(1/2);

        %prioritize by time....
        timeDiff = currentFrameINDEX - cellfun(@(x) x(end,1), traj{cellN});
        timeDiff = timeDiff(:);
        %try and link trajectories giving precendence to those with most recent
        %observations
        [~,checkOrder] = sort(timeDiff,'ascend');
        checkOrder(timeDiff(checkOrder) > gapFrames) = [];

        for chk = checkOrder(:)'
            dmin = (dxy(chk,:) < maxdXY) .* (dxy(chk,:));
            if sum(dmin) == 0, continue, end
            dmin(dmin == 0) = inf;
            [~, colN] = min(dmin);
            if isinf(dmin(colN)), continue, end
            if isempty(colN), continue, end
            traj{cellN}{chk}(end+1,:) = [currentFrameINDEX, pos{cellN}.x(colN,currentFrameINDEX), pos{cellN}.y(colN,currentFrameINDEX)];
            dxy(:,colN) = inf;
            if sum(isinf(dxy(:))) == length(dxy(:)), break, end
        end

        %which ones were not used???
        newTrajectories = find(sum(dxy,1) ~= inf);

        for nt = newTrajectories(:)'
            traj{cellN}{end+1} = [currentFrameINDEX, pos{cellN}.x(nt,currentFrameINDEX), pos{cellN}.y(nt,currentFrameINDEX)];
        end    
    end
    
    % split up trajectories with gaps of gap frames
    for t = length(traj{cellN}):-1:1
        d = diff(traj{cellN}{t}(:,1));
        splitIX = find(d > gapFrames);
        if isempty(splitIX), continue, end
        pairs2move = [splitIX(:)'+1;[splitIX(2:end)',size(traj{cellN}{t},1)]]';
        
        for r = size(pairs2move,1):-1:1
            a = pairs2move(r,1);
            b = pairs2move(r,2);
            traj{cellN}{end+1} = traj{cellN}{t}(a:b,:);
            traj{cellN}{t}(a:b,:) = [];
        end
    end
    
    % eliminate any trajectories that are too short
    tooShort = (cellfun(@length, traj{cellN}) < minTrajectoryLength).*(1:length(traj{cellN}));
    tooShort(tooShort==0) = [];
    traj{cellN}(tooShort) = [];
    
end