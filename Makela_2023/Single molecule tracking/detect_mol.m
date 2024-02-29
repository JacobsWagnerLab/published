%% Phasor based localization of fluorescent spots
% The algorithm converts the region of interest around a point spread function 
% to two phase vectors (phasors) by calculating the first Fourier coefficients 
% in both the x- and y-direction. The angles of these phasors are used to localize 
% the center of the fluorescent emitter, and the ratio of the magnitudes of the
% two phasors can be used for detection of assymetry in spots (e.g. astigmatism).

% File formats .tif (multi-page) or .nd2 (needs bioformats matlab package)

% Cite: Martens et al., J. Chem. Phys. 148, 123311 (2018)
% https://doi.org/10.1063/1.5005899

% 6/23/22 Jarno Makela

%% Determine detection threshold

% select a .nd2 or .tif file
[movieFilename, moviePathname] = ...
    uigetfile({'*.nd2';'*.tif'}, 'MatFile data:','MultiSelect', 'off');
appData.testFile.movieFilename = movieFilename;
appData.testFile.moviePathname = moviePathname;
filetype = movieFilename(end-3:end);

% parameters for detection
appData.localizationWindow = 3;         % window size for detection
appData.localizationThresh = 100;        % threshold for detection
appData.testDetectionFrame = 10000;      % test frame number

if ~(isnumeric(movieFilename)&&movieFilename==0) && strcmp(filetype,'.tif')
    % open image
    singleImage = double(imread([moviePathname movieFilename], appData.testDetectionFrame));

    % show scaled image
    figure
    hold all;
    imshow(singleImage,[1.02*min(1.02*min(singleImage)) 0.94*max(0.94*max(singleImage))]);
    if ~(isnumeric(movieFilename)&&movieFilename==0) %check the user has not pressed cancel  
        [locX, locY] = phasorLocalizationTiffSingle([moviePathname movieFilename],...
            appData.localizationThresh, appData.localizationWindow,...
            appData.testDetectionFrame);    
    end
    % plot spots
    if ~any(isnan([locX; locY]))
        hold all
        plot(locX,locY,'ro','MarkerSize',10);
        xlabel('x [pixels]');
        ylabel('y [pixels]');
        axis image;
    end

elseif ~(isnumeric(movieFilename)&&movieFilename==0) && strcmp(filetype,'.nd2')
    % open .nd2 file
	reader = bfGetReader([moviePathname movieFilename]);
	singleImage = double(bfGetPlane(reader, appData.testDetectionFrame));

    % show scaled image
    figure
    hold all;
    imshow(singleImage,[1.02*min(1.02*min(singleImage)) 0.94*max(0.94*max(singleImage))]);
    if ~(isnumeric(movieFilename)&&movieFilename==0) %check the user has not pressed cancel  
        [locX, locY] = phasorLocalizationND2Single([moviePathname movieFilename],...
            appData.localizationThresh, appData.localizationWindow,...
            appData.testDetectionFrame);    
    end
    % plot spots
    if ~any(isnan([locX; locY]))
        hold all
        plot(locX,locY,'ro','MarkerSize',10);
        xlabel('x [pixels]');
        ylabel('y [pixels]');
        axis image;
    end
end
 
%% Run detection for all frames

% parameters
appData.localizationWindow = 3;
appData.localizationThresh = 100;

% select .nd2 or .tif file(s)
[movieFilename, moviePathname] = ...
    uigetfile({'*.nd2';'*.tif'}, 'MatFile data:','MultiSelect', 'on');

if ~(isnumeric(movieFilename) && movieFilename == 0) % check for cancel
    % extract number of files
    info = whos('movieFilename');
    if strcmp(info.class,'char')
        appData.nFiles = 1;
    else
        appData.nFiles = numel(movieFilename);
    end

    for ii = 1:appData.nFiles % loop over files
        % path with either 1 or more files
        if appData.nFiles == 1
            loadname  = [moviePathname movieFilename];
        else
            loadname = [moviePathname movieFilename{1,ii}];
        end
        
        % get filetype
        filetype = loadname(end-3:end);

        % process either .tif or .nd2 files
        if strcmp(filetype,'.tif')
            phasorLocalizationTiff(loadname,...
                appData.localizationThresh, appData.localizationWindow);
 
        elseif strcmp(filetype,'.nd2')
            phasorLocalizationND2(loadname,...
                appData.localizationThresh, appData.localizationWindow);
        end
    end
end