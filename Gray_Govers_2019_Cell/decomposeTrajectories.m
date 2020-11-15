function traj = decomposeTrajectories(SpotDat, traj)
%{
-About-
Transforms spot localizations in image coordinates to spot localizations in
cell coordinates

decomposes positions (x,y) found in
traj{cellNumber}{trajectoryNumber}(frameNumber,imageXPos,imageYpos)
relative the cellMesh found in SpotDat{Frame}{cellNumber}

-Inputs-
SpotDat: output from particleTracking.m
traj: output from pos2trajectories.m

-varargin-
na

-Outputs-
traj: A nested cell array containing trajectories. traj{1} contains
trajectories from cell #1 and traj{1}{1} is the 1st trajectory identified
in cell #1. Each trajectory is n-rows and columns:
1) frame #
2) image x coordinate
3) image y coordinate
4) fractional cell length location
5) half the cell length
6) fractional cell width
7) cell width

-Example-
   
-Supplementary-

-Dependencies-
particleTracking.m
SpotDat2Pos.m
pos2trajectories.m
decomposeTrajectories.m

-References-
n/a

-Author-
Brad Parry, 2014 July 28
%}


%polar region of the cell to exclude from centerline estimates
poleRegion = 6;

for C = 1:length(traj)
    for T = 1:length(traj{C})
        if isempty(traj{C}{T}), continue, end
        traj{C}{T}(end,7) = 0;
        for k = 1:size(traj{C}{T},1)
            Frame = traj{C}{T}(k,1);

            [percentCL, halflength, percentCW, cellwidth, ~] = cellProjection(SpotDat{Frame}{C}.mesh, traj{C}{T}(k,2:3), poleRegion);
            if isnan(percentCL)
                %print frames and cells that hvae an issue
                Frame
                C
            end
            traj{C}{T}(k,4:7) = [percentCL, percentCW, halflength, cellwidth];
            
        end
    end
end