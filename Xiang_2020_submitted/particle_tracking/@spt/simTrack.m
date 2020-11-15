%% Simulate Brownian trajectories with a known diffusion coefficient.
%
% -About-
%   The method simTrack() simulates trajectories of Brownian particles with
%   a known diffusion coefficient. The user can specify the number of
%   particles, the number of steps taken by each particle, the time stamps,
%   and the diffusion coefficient.
%
%   If any trajectory length distribution was stored beforehand in the
%   object, the method uses the this distribution for the number of
%   particles and the steps they take.
%
%   If no trajectory length distribution was found, the method simulates
%   100 particles that take the number of steps that is drawn from an
%   exponential distribution with a mean of 100. Both of the number of
%   particles and the mean of the number of steps taken can be set by the
%   user, as long as no trajectory length was stored. The method prevents
%   the user from changing the number of particles and the mean of the
%   steps they take, if a trajectory length distribution existed.
%
%   For the record, the method prints the number of particles, the average
%   number of steps they took, the time between the steps, and the known
%   diffusion coefficient to the command window.
%
% -Input-
%   - obj: spt object
%   
% -Varargin-
%   - np: number of particles
%   - ns: average number of steps taken
%   - dt: time between consecutive steps [s]
%   - D:  diffusion coefficient [µm^2/s]
%
% -Output-
%   - tracks: particle trajectories. A cell array with a length of np.
%             Each cell contains a matrix of ns x 3. The columns from left
%             to right presents the frame number, the postions in x and y
%             dimensions. The units are not converted to physical units.
%             [frame, px, px]
% -Example-
%   % Simulate 100 particles that take an average of 100 steps
%   myParticle.simTrack();
%   % Simulate 10 particles that take an average of 5000 steps
%   myParticle.simTrack('np',10,'ns',5000);
%   % Simulate 10 parrticles that take an average of 5000 steps with a
%   % diffusion coefficient of 2 µm^2/s
%   myParticle.simTrack('np',10,'ns',5000, 'D', 2);
%
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

function simTrack(obj,varargin)

% Default parameters for the simulation
% If the trajectory length was specified beforehand, use it. Otherwise,
% set the number of particles to be 100 and draw the number of steps from
% an expoential distribution with a mean of 100.
if ~isempty(obj.trackLen)
    % Number of particles
    np = length(obj.trackLen);
    % Number of steps each particle takes
    ns = obj.trackLen;
else
    % Number of particles
    np = 100;
    % Mean of the exponential distribution from which the number of steps
    % to be generated
    ns_mean = 100;
end

% Diffusion coeffcient of particles [µm^2/s]
D = 1;
% Time stamps [sec]
dt = 0.01;
% Pixel length [µm/px]
pixelLength = 0.0642;

% Overwrite the default with user input
% If the track length exists, warn the user about it and ignores their
% input
for ii = 1: 2: length(varargin)
    switch lower(varargin{ii})
        case 'np'
            if ~isempty(obj.trackLen)
                warning('Trajectory length distribution already existed, the input was ignored.');
            else
                np = varargin{ii+1};
            end
        case 'ns'
            if ~isempty(obj.trackLen)
                warning('Trajectory length distribution already existed, the input was ignored.');
            else
                ns_mean = varargin{ii+1};
            end
        case 'dt'
            dt = varargin{ii+1};
        case 'd'
            D = varargin{ii+1};
        otherwise
            error('Invalid argument. Valid options: np|ns|dt|D.');
    end
end

% If there is no track length input, random number of steps are taken
if isempty(obj.trackLen)
    % Finalize the number of takes each particle takes
    ns = round(exprnd(ns_mean,1,np));
end

% Displacement length per step [µm]
% Note: mathematically, this is the standard deviation of the displacement
% distribution.
k = sqrt(2 * D * dt);

% Memory allocation to store simulated tracks
sTracks = cell(np,1);

% For each particle, simulate Gaussian-distributed displacements in two
% dimensions and cumulatively sum the displacements up. The initial
% positions were specially treated to ensure a dispersion of the particles.

% The dispersion / spacing between initial positions of tracks. This
% dispersion is solely for visual purposes.
dispersion = 100*sqrt(4*D*dt*max(ns));

% Start simulation
for ii = 1 : np
    
    % Time stamps [frame]
    time = (0:ns(ii)-1)';
    
    % Displacements in two dimensions with a normal distribution and
    % calculated displacement length. [px]
    dx = k .* randn(ns(ii), 2) ./ pixelLength;
    
    % Replace with the newly generated intial positions.
    dx(1,:) = rand(1,2).*dispersion;
    
    % Integrate the displacements to generate tracks [px]
    x = cumsum(dx, 1);
    
    % Store the tracks in pre-defined memory.
    sTracks{ii}  = [time, x];
end

obj.pixelLength = pixelLength;
obj.frameTime = dt;
obj.tracks = sTracks;
obj.rmShort(3);

% If either MSD or meanMSD already existed, they store the results based on
% previous tracks. Update the MSD results below.
if ~(isempty(obj.MSD) || isempty(obj.meanMSD))
    obj.getMSD();
    obj.getMeanMSD();
end

fprintf('Simulated %d tracks. Each takes an average of %.2f steps.\n',length(obj.tracks),mean(ns));
fprintf('Frame rate: %.2f sec; Diffusion coefficient: %.2f um^2/s\n',dt,D);

end