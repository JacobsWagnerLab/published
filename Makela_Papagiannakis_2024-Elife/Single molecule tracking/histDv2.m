function [D, long_enough_tracks] = histDv2(tracks,params)
% calculates the diffusion coefficient
% for each track in "tracks" and plots a histogram of D.

% parameters
pixel = params.pixel;                   % length per pixel
dT = params.dT;                         % time per frame
sigmaNoise = params.sigmaNoise;         % localization noise
DhistMinSteps = params.DhistMinSteps;   % minimum number of steps for a track to be analyzed

% find the elements in all tracks which have a single frame step between the element and consecutive element
% i.e. ignore blinking steps in tracks.
single_frame_steps= find((tracks(2:end,3) - tracks(1:end-1,3))==1);
% single_frame_steps=1:length(tracks);
% find elements in tracks from the same molecules as consective element
same_mol_steps= find((tracks(2:end,4) - tracks(1:end-1,4))==0);
% intersect of both gives the single frame steps only from the same molecule
single_frame_same_mol = intersect( same_mol_steps, single_frame_steps);
% accept all tracks with blinking
% single_frame_same_mol = find((tracks(2:end,4) - tracks(1:end-1,4))==0);

% Calculate the squared dispacement form these elements in tracks
single_frame_sq_displacements = sum((tracks(single_frame_same_mol,1:2) - tracks(single_frame_same_mol+1,1:2)).^2,2);

% molecule numbers from only the desired track elements
single_frame_tracks = tracks(single_frame_same_mol,4);
% histogram of mol numbers gives the step number of each desired track
n = histc(single_frame_tracks,1:max(single_frame_tracks));
% find tracks with  enough steps
long_enough_tracks = find(n>=DhistMinSteps);

%pre allocate variables
MSD_all = zeros(length(long_enough_tracks),1);
track_lengths = zeros(length(long_enough_tracks),1);

% if tracks are truncated
if params.truncated_tracks == 1
    % loop over all long enough tracks
    for ii = 1:length(long_enough_tracks)
        % find the elements of the track
        xx = find(single_frame_tracks == long_enough_tracks(ii));
        
        % calculate the mean of the squared displacements and store values in
        % array. Truncated tracks so only take first DhistMinSteps steps
        MSD_all(ii,1) = mean(single_frame_sq_displacements(xx(1:DhistMinSteps)));
        
        % save the lenght of the track
        track_lengths(ii) = numel(xx(1:DhistMinSteps));
    end
    
else
    for ii = 1:length(long_enough_tracks)
        % find the elements of the track
        xx = find(single_frame_tracks == long_enough_tracks(ii));
        
        % calulate the mean of the squared displacements and store values in array
        MSD_all(ii,1) = mean(single_frame_sq_displacements(xx));
        
        %save the lenght of the track
        track_lengths(ii) = numel(xx);
    end
end

% convert into um
MSD = MSD_all * pixel^2; % convert from pixel to length units

% calculate D from MSD and correct for localization noise
D = MSD/(4*dT) -sigmaNoise^2*pixel^2/dT;
                           
end