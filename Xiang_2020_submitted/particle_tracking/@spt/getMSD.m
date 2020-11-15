%% Calculate the mean squared displacement for each track.
%
% -About-
%   The method getMSD() calculates the mean squared displacements for every
%   track stored in the object. The particle positions of a track are first
%   extracted. Based on its track length, all the possible delays are
%   defined. That is 0 to (track length -1). At each delay, the
%   displacement as well the mean squared displacements (MSD) are
%   calculated. All the MSDs are consolidated together and passed back to
%   the object.
% 
% -Input-
%   - obj: spt object
% 
% -Output-
%   - MSD: mean squared displacements for individual particles. A Nx1 cell
%          array, where N is the number of particles. Within each cell,
%          there is a Nx3 matrix, where N is the number of delays equal
%          exactly to the track length, since there is a 0 delay. The three
%          columns represents the time delay, the MSD at that delay and the
%          number of delays used to calculate this MSD. The units are not
%          converted to physical units. [frame, px^2, Unity]
% 
% -Example-
%   % Calculate the mean squared displacements based on tracks
%   myParticle.getMSD();
% 
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University


function getMSD(obj)

% Memory for all MSDs, a cell array is used since the track length can vary
MSD = cell(length(obj.tracks),1);

% Loop through individual tracks
for ii = 1 : length(obj.tracks)
    
    % Extract the 2nd and 3rd columns from a track
    % These are the positions of a particle in x and y dimensions
    track = obj.tracks{ii}(:,[2,3]);
    
    % Length of this track
    N = length(track);

    % All possible delays (0 to track length - 1)
    delays = 0:N-1;
    
    % Memory for MSD of this track
    msd = zeros(length(delays),1);
    
    % Memory for the number of delays at each time delay
    nDelay = zeros(length(delays),1);
    
    % Loop through individual delays
    for jj = 1 : length(delays)
        
        % Extract a specific delay
        delay = delays(jj);
        
        % Displacement Nx2 matrix, [x,y]
        dr = track(1:end-delay,:) - track(1+delay:end,:);
        
        % Squared displacement: dr^2 = x^2 + y^2
        % Mean squared displacement 
        msd(jj) = mean(sum(dr.^2,2));
        
        % Number of delays used to average MSD
        nDelay(jj) = size(sum(dr.^2,2),1); 
    end
    
    % Consolidate all MSDs calculated
    MSD{ii} = [delays' msd nDelay];
end

% Pass the result back to the object
obj.MSD = MSD;

end
