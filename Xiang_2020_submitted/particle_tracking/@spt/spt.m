%% Class definition.
%
% -About-
%   The class SPT is used for the analysis of single particle tracking The
%   user should have the data gathered from the particle tracking
%   experiments, that is, images of particle locations with known
%   aquisition rates. This class comes with a libaray of useful methods
%   (functions) that can be directly applied to the data stored in an
%   instance. The main purpose of the SPT is to probe the material/physical
%   environment properties via studying the motion of tracers.
%
% -Inputs-
%   No required input. Without any input, a default instance will be
%   created.
%
% -Varargin-
%   Specified by corresponding keys. Keys listed below. Case insensitive.
%   When a key is used, all information (except optional) must be
%   specified. For example, if a stackPath is specified, particleSize,
%   pixelLength and frameTime must be specified altogether.
%   - date:         date of the experiment                (optional)
%   - about:        description about the experiment      (optional)
%   - stackPath:    full path to the image stack
%   - particleSize: size of particles used in the experiment [µm]
%   - pixelLength:  camera-specific length of a single pixel [µm/px]
%   - frameTime:    image aquisition rate [s]
%
% -Output-
%   - locs:         locations of particles identified [px px frame]
%   - tracks:       tracks deduced from particle locations [frame px px]
%   - MSD:          mean squared displacement [px^2]
%   - meanMSD:      ensemble averaged MSD [px^2]
%   - alpha:        slopes of individual log-log MSD curves
%   - D:            diffusion coefficents [µm^2/s]
%   - ensembleAlpha slope of the ensemble-averaged MSD in the log-log scale
%   - ensembleD     diffusion coefficient based on the ensemble MSD [µm^2/s]
%   - trackLen:     track lengths of individual tracks
%   - displacement: displacements of particles with a specific delay [µm]
%   - Rg:           radius of gyration of individual tracks [µm]
%   - intercepts:   y-intercepts of the log-log MSD curves [µm^2]
%   - angles:       angles between consecutive displacements [radian]
%   - meanAngles:   mean of angles between displacements per track [radian]
%   - cntrLocs:     particle locations from tracks centered at origin [µm]
%   - gaussians:    fitting parameters from the bivariate normal fit
%   - intensity:    particle intensities calculated from the Gaussian fits
%
% -Keywords-
%   particle tracking, diffusion coefficient, viscosity
%
% -Dependencies-
%   bpass.m; track.m; gaussfit2d.m;
%
% -Example-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   For a real experiment:
%   % Create a instance for real experiment
%   myParticle = spt('date',        '2017-AUG-01', ...
%                    'about',       'This is an example.', ...
%                    'stackPath',   '~/Desktop/myStack.tif',
%                    'particleSize', 0.1, ...
%                    'pixelLength',  0.0642, ...
%                    'frameTime',    0.01);
%   % Sample a frame to determine the intensity threshold and object size
%   myParticle.sampleFrame(10);
%   % Localization
%   myParticle.getLocation('minSpotSize',        10, ...
%                          'maxSpotSize',        30, ...
%                          'eccentricity',       0.5, ...
%                          'intensityThreshold', 10, ...
%                          'startFrame',         1, ...
%                          'endFrame',           1000);
%   % Check the localizations
%   myParticle.plotLocation();
%   % Link the particle locations
%   maxDisplacement = 5;
%   minTrackLen = 10;
%   myParticle.getTrack(maxDisplacement, minTrackLen);
%   % Plot the tracks to check linking
%   myParticle.plotTracks();
%   % Calculate the MSD
%   myParticle.getMSD();
%   % Plot the MSD cures
%   myParticle.plotMSD();
%   % Calculate the ensemble-averaged MSD
%   myParticle.getMeanMSD();
%   % Plot the ensemble-averaged MSD
%   myParticle.plotMeanMSD();
%   % Fit MSD curves with the 1st 10 points
%   myParticle.fitMSD(1:10);
%   % Plot the distribution of diffusion coefficients 
%   myParticle.plotD();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   You can also simulate the tracks:
%   % Create a default instance
%   myParticle = spt();
%   % Simulate some tracks
%   myParticle.simTrack();
%   % Calculate the MSD
%   myParticle.getMSD();
%
%
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

classdef spt < handle
    %% Private attributes
    % These can be set by the user when constructing the object but later
    % cannnot modify them. This helps to prevent unwanted changes. However,
    % the user can check these using the getInfo() specified below.
    
    properties (Access = private)
        
        % Date of the experiment (optional)
        date

        % Information regarding to the experiment (optional)
        about
        
        % Absolute path to the image stack
        stackPath
        
        % Size of particles used [µm]
        particleSize
        
        % Size of the pixel length of the camera used [µm]
        pixelLength
        
        % Time between two consecutive frames [s]
        frameTime
 
        % Calculated estimation of the object size based on the particle
        % size provided.
        estimatedSize
        
        % Number of frames in the image stack
        imNum
        
        % Intensity threshold used for bandpass [A.U.]
        iThreshold
        
        % Minimum spot size for localization [px]
        minSpotSize
        
        % Maximum spot size for localization [px]
        maxSpotSize
        
        % Eccentricity threshold for localization [A.U.]
        eccentricity
        
    end
    
    
    %% Public attributes
    % These data can be accessed and moified directly by the user. These
    % data are the basic results of the particle tracking.
    
    properties(Access = public)
        
        % Locations of particles identified [px px frame]
        % N by 3 arrays, where N depends on the number of particles
        % localized. The first two columns specify the x and y positions of
        % the localizations in that specific frame. The third column
        % specifies the frame number of the images, in which the particles
        % are localized.
        locs
        
        % Particle trajectories deduced from locations [frame px px] 
        % Cell arrays with a length equal the number of particles tracked.
        % Each cell, representing an individual track, is a N by 3 matrix
        % (N is the track length), with three columns (from left to right)
        % specifying the time, x positions and y positions of the particle.
        tracks
        
        % Mean squared displacement (MSD) [px^2]
        MSD
        
        % Ensemble-averaged MSD for the corresponding tracks [px^2]
        meanMSD
        
        % Slopes of individual log-log MSD curves
        alpha
        
        % Diffusion coefficents of individual linear MSD curves [µm^2/s]
        D
        
        % Slope of the ensemble-averaged MSD in the log-log scale
        ensembleAlpha
        
        % Diffusion coefficient fitted from the ensemble MSD [µm^2/s]
        ensembleD
        
        % Lengths of individual tracks
        trackLen
        
        % Particle displacements with a specific delay [µm]
        displacement
        
        % Radius of gyration of individual tracks [µm]
        Rg
        
        % Y-intercepts of the log-log MSD curves [µm^2]
        intercepts
        
        % The angles between consecutive displacements [radian]
        angles
        
        % Mean of angles between displacements per track [radian]
        meanAngles
        
        % Particle locations from tracks centered at origin [µm]
        cntrLocs
        
        % Fitting parameters from the bivariate normal fit
        gaussians
        
        % Particle intensities
        intensity
        
    end
    
    %% Function prototypes
    methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constructors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = spt(varargin)
            % Create an instance of the spt class.
            % If there is any user-defined input, catch this information to
            % specify how the object will be constructed. If no information
            % was given, the create an emtpy instance by default.
            if ~isempty(varargin)
                for ii = 1:2:length(varargin)
                    switch lower(varargin{ii})
                        % Date of the experiment (optional)
                        case 'date'
                            obj.date = varargin{ii+1};
                        % Basic information about the experiment (optional)
                        case 'about'
                            obj.about = varargin{ii+1};
                        % Full path to the .tif image stack
                        case 'stackpath'
                            obj.stackPath = varargin{ii+1};
                        % Known particle size [µm]
                        case 'particlesize'
                            obj.particleSize = varargin{ii+1};
                        % Camera-specific pixel size [µm/px]
                        case 'pixellength'
                            obj.pixelLength = varargin{ii+1};
                        % Time between consecutive frames, i.e. frame rates
                        case 'frametime'
                            obj.frameTime = varargin{ii+1};
                        % Any input key besides the ones specified above
                        % will be regarded as an error. Offer the user
                        % the keys available.
                        otherwise
                            error(['Invalid argument: %s. ',...
                                'Valid options: stackPath|particleSize|pixelLength|',...
                                'frameTime|date(optional)'],varargin{ii});
                    end
                end
                
                % Check input data types
                if ~(ischar(obj.stackPath) && isnumeric(obj.particleSize) && ...
                        isnumeric(obj.pixelLength) && isnumeric(obj.frameTime))
                    error('Invalid argument data type.')
                end
                
                % Provide an estimate of the particle size in terms of
                % pixel sizes.
                obj.estimatedSize = ceil(obj.particleSize / obj.pixelLength + 3);
                
                % Validate if the image stack can be opened.
                try
                    obj.imNum = length(imfinfo(obj.stackPath));
                catch
                    error(['Unable to open specified image stack @ ',...
                        obj.stackPath,'. Please verify the path.'])
                end
            end
        end
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%% Function prototypes %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Below are forward declarations of all the methods associated with
        % this class. I grouped them into 5 categories: accessories,
        % setters, getters, fitters and plotters. Their uses are briefly
        % described here. For more details, check the comments in
        % individual functions.
        
        %%%%%%% Accessories %%%%%%%%
        % Show speicfic frames from the image stack
        showFrame(obj,varargin);
        % Sample the 1st frame to help determine object intensity and size
        sampleFrame(obj,intensityThreshold);
        % Simulate Brownian tracks with a known diffusion coefficient
        simTrack(obj,varargin);
        % Remove short trajectories
        rmShort(obj,th);
        
        
        %%%%%%%%% Mutators %%%%%%%%%%
        % Set privately stored information about the instance
        setInfo(obj,varargin);
        
        
        %%%%%%%%% Accessors %%%%%%%%%%
        % Print the basic information about the instance
        varargout = getInfo(obj);
        % Calculate the length of all tracks
        getTrackLen(obj)
        % Determine the particle locations
        getLocation(obj, varargin);
        % Link particle locations into tracks
        getTrack(obj,maxDist,minTrackLen);
        % Calculate mean squared displacements
        getMSD(obj);
        % Calculate the ensemble-averaged MSD
        getMeanMSD(obj);
        % Calculate angles between consecutive displacements
        getAngle(obj);
        % Calculate the displacement distribution with a specific delay
        getDisplacement(obj,varargin);
        % Calculate the radius of gyration
        getRg(obj);
        % Calculate the y-intercepts in the log-log MSD plot 
        getIntercept(obj,delay);
        % Collapse all positions from origin-centered tracks
        getCntrLoc(obj);
        % Calculate the intensity of each particle
        getIntensity(obj);
        % Filter the tracks based on various criteria
        filterTrack(obj,criterion,range,reverse);
        % Choose a color for plotting
        color = selectColor(~,key);
        
        
        %%%%%%%%% Fitters %%%%%%%%%%
        % Fit the MSD curves based on methods (linear, loglog)
        fitMSD(obj,fitRange,method,ensemble)
        
        %%%%%%%%% Plotters %%%%%%%%%
        % Plot the particle locations on top of the frames
        plotLocation(obj,varargin)
        % Plot the individual MSD curves
        plotMSD(obj,varargin)
        % Plot the ensemble-average MSD curve
        plotMeanMSD(obj)
        % Plot the individual tracks
        plotTracks(obj,varargin)
        % Plot the displacement distribution (displacement or distance)
        plotDisplacement(obj,varargin)
        % Plot the kernel density distribution of alpha values
        plotAlpha(obj,varargin)
        % Plot the kernel density distribution of diffusion coefficients
        plotD(obj,varargin)
        
    end
end