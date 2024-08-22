%% Track single molecules
% The algorithm tracks spots from localization data for each individual 
% cell area. Uses a maximum tracking window and allows 1 frame 
% disappearance inside a track. Multiple localizations within tracking
% window are removed. Saves data as .mat with each row corresponding to a
% single cell data:
%   localization data
%   tracks
%   number of molecules
%   cell area data from MicrobeTracker

% Following files are required:
% 1. MicrobeTracker file for cell meshes (*.mat)
% 2. Localization file from spot detection (*.out)
% 3. Brightfield image file (*.tif)

% 7/7/22 Jarno Makela

%% Choose files

[meshes_filename, meshes_pathname] = uigetfile('*.mat', 'Microbetracker meshes: ','MultiSelect', 'off');
if ~(isnumeric(meshes_filename)&&meshes_filename==0) % check the user has not pressed cancel
    % load files
    load([ meshes_pathname meshes_filename], 'cellList');
else
    return
end

% select localisation files
[appData.dataFilename, appData.dataPathname] =  ...
    uigetfile([meshes_pathname '.out'], 'localisation data to analyse:','MultiSelect', 'off');
if ~(isnumeric(appData.dataFilename)&&appData.dataFilename==0) % check the user has not pressed cancel
    data_filename = [appData.dataFilename appData.dataPathname];
else
    return
end

% select brightfield file
[appData.brightfield_filename, appData.brightfield_pathname] =  ...
    uigetfile([meshes_pathname '.tif'], 'Bright field files:','MultiSelect', 'off');
if ~(isnumeric(appData.dataFilename)&&appData.dataFilename==0) % check the user has not pressed cancel
    brightfield_filename = [appData.brightfield_pathname appData.brightfield_filename];
else
    return
end

%% Adjust pixel offset
% Fluorescence image and brightfield image might not be perfectly aligned.
% Especially if different z-positions used. This will only affect ROI where
% spots are tracked and not spot coordinates.

% This removes empty cells in cellList from MicrobeTracker
jj = 1;
for ll = 1:length(cellList{1})
    if ~isempty(cellList{1}{ll}) && cellList{1}{ll}.mesh(1,1) ~= 0
        cellList_allfull{1}{jj} = cellList{1}{ll};
        jj = jj + 1;
    end
end

% open gui for adjusting offset
appData.pixelshift = pixel_shift_gui(cellList_allfull{1}, appData.dataFilename, brightfield_filename);  

%% Track molecules

% parameters
appData.trackParams.maxDisp = 5;    % maximum distance (px) that a particle would move in a single time interval
appData.trackParams.mem = 1;        % number of time steps that a particle can be 'lost' and then recovered again
appData.trackParams.good = 0;       % eliminate all trajectories with fewer than X positions (all tracks for now)
appData.trackParams.dim = 2;        % set to equal to the dimensionality of the coordinate data

disp('Segmenting and tracking cells');
[movie_data] = ROI_tracks_microbetracker_v2(appData,cellList_allfull{1},appData.dataFilename,brightfield_filename);
disp('Done!');

% add pixelshift information to movie_data
if ~isempty(appData.pixelshift)
    movie_data.pixelshift = appData.pixelshift;
end

% save tracked data
save([ meshes_pathname meshes_filename(1:end-4) 'exp_locoli.mat' ],'movie_data');

%% Remove multiple localizations that are closer than a minimum distance

min_dist = 6;       % minimum distance in pixels between molecules 

% choose files
[filename,pathname] = uigetfile('*locoli.mat', 'Select','Select','MultiSelect', 'on');
if ~iscell(filename)
    name = filename;
    filename = cell(1);
    filename{1} = name;
end

for qq = 1:length(filename)   
    load([pathname filename{qq}])

    for kk = 1:length(movie_data.cellROI_data)
        tracks = movie_data.cellROI_data(kk).tracks;

        % find repeated elements in time frames
        edges = min(tracks(:,3)):(max(tracks(:,3))+1);
        [counts, values] = histcounts(tracks(:,3), edges);
        values = values(1:end-1); % remove last extra bin
        repElem = values(counts > 1);

        % loop from end to beginning
        for ii = length(repElem):-1:1
            % measure distances and keep indeces only with smaller than min_dist
            inds = find(repElem(ii) == tracks(:,3));
            if ~isempty(inds) && length(inds) > 1
                X = tracks(inds,1);
                Y = tracks(inds,2);
                locs = [X Y];
                d = pdist(locs,'euclidean');
                Z = squareform(d);
                Z(Z == 0) = NaN;
                % from end to beginning to avoid indexing problem
                for hh = size(Z,1):-1:1
                    if min(Z(hh,:)) < min_dist
                        % delete duplicated track index
                        tracks(inds(hh),:) = [];
                    end
                end
            end
        end
        % find missing mol ID by finding larger jumps in mol ID
        diff_mol_ID = diff(tracks(:,4));
        gap_ind = find(diff_mol_ID > 1);
        
        % loop over too large gaps
        for jj = 1:length(gap_ind)
            % subtract indeces after the gap by gap size - 1
            tracks((gap_ind(jj)+1):end,4) = tracks((gap_ind(jj)+1):end,4) ...
                - (diff_mol_ID(gap_ind(jj))-1);
        end
        
        % insert tracks back to movie_data
        movie_data.cellROI_data(kk).tracks = tracks;
    end
    % save movie_data
    fname_splitted = strsplit(filename{qq},'locoli.mat');
    name = char(strcat(fname_splitted{1},'filt_locoli.mat'));
    save([pathname name],'movie_data')
    disp('.')
end

