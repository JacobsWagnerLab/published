%{
-About-
This script pipelines several functions to: 1) obtain particle
localizations from fluorescent images, 2) link the localizations into
trajectories 3) Project localizations from all trajectories into cellular
coordinates (i.e., transform coordinates to position along cell length,
width)

Parameters explicitly set below represent defauly parameters used in Gray
and Govers 2018

-Inputs-
The first fluorescent image location from a time-lapse, corresponding to
frame 1 of a cellList (generated from Oufti)
The cellList that corresponds to fluorIM

-varargin-
na

-Outputs-
SpotDat: nested structure of cell arrays {Frame}{Cell}.(tracking output)
with fields (the kth element of each field corresponds to the kth spot):
    Centroid: xy coordinates of spot peaks
    rmse: root mean-squared error of the fit
    sigma: sigma of the fitted gaussian function
    gtr: a measure of spot intensity over the background
    confint: confidence intervals of the
    fitdata: statistics of the fit
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

%before begining, load a cellList into the matlab workspace

%location of the first fluorescent image in a series and corresponding to
%cells in cellList.meshData{1}
fluorIM = '';

shift_cellList_meshes = 0; %def = 0
def_params.pixelsPerSpotRange = [9 50]; %def = [9 50]
def_params.spotlength = 9; %def = 9
def_params.frequency = .6; %def = 0.6
def_params.ithresh = .25; %def = 0.25
def_params.nSigmaCutoff = 3; %def = 3
def_params.fitTol = [-inf inf;-inf inf;-inf inf;-inf inf;-inf inf];

[SpotDat, def_params] = particleTracking_WTG(x, cellList, def_params);
pos = SpotDat2Pos(SpotDat); 

maxdXY = 30; 
gapFrames = 5;
minTrajectoryLength = 50; 
[traj] = pos2trajectories(pos, maxdXY, gapFrames, minTrajectoryLength);
traj = decomposeTrajectories(SpotDat, traj);

