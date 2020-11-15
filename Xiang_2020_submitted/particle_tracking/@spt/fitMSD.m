%% Fit MSD curves based on the specific methods (linear, loglog, ensemble).
%
% -About-
%    Fit the MSD curves (individual or ensemble-averaged) using various
%    methods. The results from the fitting are then used to calculate the
%    diffusion coefficients and alpha values depending on the fitting
%    method used.
%
% -Inputs-
%   - obj:      spt object
%   - fitRange: range of fit (e.g. 1:10, first ten points)
%   - method:   linear, loglog, specified as a string
%   - ensemble: boolean, if true, only fit the ensemble MSD, default: false
%
% -Output-
%   - D:      diffusion coefficients [µm^2/s]
%   - alpha:  alpha values
%
% -Example-
%   % Calculate the diffusion coefficient by fitting the linear MSDs
%   - myParticle.fitMSD(1:10,'linear');
%   % Calculate the alpha values by fitting the loglog MSDs
%   - myParticle.fitMSD(2:11,'loglog');
%   % Print out the alpha values by fitting the ensemble loglog MSD
%   - myParticle.fitMSD(2:11,1);
%
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

function fitMSD(obj,fitRange,method,ensemble)

% If no fitRange is specified, the range is determined by the shortest
% track. If no mode is specified, the default is linear.
if nargin < 3
    % Calculate the track length
    obj.getTrackLen();
    % Set the default fit range to be the shortest track length
    fitRange = 1:min(obj.trackLen);
    % Default fitting method
    method = 'linear';
    % Default to fit individual MSD
    ensemble = 0;
elseif nargin < 4
    % Default to fit individual MSD
    ensemble = 0;
end

% If no MSD or ensemble-averaged MSD were calculated and stored before,
% calculate them now.
if isempty(obj.MSD) || isempty(obj.meanMSD)
    % Calculate the MSD
    obj.getMSD();
    % Calculate the ensemble-averaged MSD
    obj.getMeanMSD();
end

% If ensemble is true, only fit the ensemble-averaged MSD
if ensemble
    fprintf('Fitting the ensemble-averaged MSD... \n');
    switch lower(method)
        case 'linear'
            fitRes = polyfit(obj.meanMSD(fitRange,1).*obj.frameTime,...
                obj.meanMSD(fitRange,2).*obj.pixelLength^2,1);
            obj.ensembleD = fitRes(1)/4;
        case 'loglog'
            if (fitRange(1) == 1)
                fitRange(1) = 2;
                warning('loglog fit can''t operate on delay 0, skip the 1st point now...')
            end
            fitRes = polyfit(log10(obj.meanMSD(fitRange,1).*obj.frameTime),...
                log10(obj.meanMSD(fitRange,2).*obj.pixelLength^2),1);
            obj.ensembleAlpha = fitRes(1);
        otherwise
            error('Invalid fitting method. Valid options: linear|loglog.');
    end

% Fit individual MSD curves
else
    fprintf('Fitting the individual MSD curves... \n');
    % Switching between different fitting methods
    switch lower(method)
        % Fit linear MSD to determine diffusion coefficients
        case 'linear'
            % Memory for diffusion coefficients
            D = zeros(length(obj.MSD),1);
            % Loop through MSD curves
            for ii = 1 : length(obj.MSD)
                % Extract a single MSD curve
                thisMSD = obj.MSD{ii};
                % Convert the time and MSD into physical units
                t = thisMSD(:,1).*obj.frameTime;
                msd = thisMSD(:,2).*obj.pixelLength^2;
                % Linear regression based on the fit range
                fitRes = polyfit(t(fitRange),msd(fitRange),1);
                % Store the quarter of the slope as the diffusion coefficient
                % Assumed 2D particle tracking here
                D(ii) = fitRes(1)/4; % [µm^2/s]
            end
            % Pass the result into the instance
            obj.D = D;
            
            % Fit loglog MSD to determine alpha values
        case 'loglog'
            % Memory for alpha values
            alpha = zeros(length(obj.MSD),1);
            % Loop through MSD curves
            for ii = 1 : length(obj.MSD)
                % Extract a single MSD curve
                thisMSD = obj.MSD{ii};
                % Convert the time and MSD into physical units
                t = thisMSD(:,1).*obj.frameTime;
                msd = thisMSD(:,2).*obj.pixelLength^2;
                % If the left bound of the fitRange is 1, warn the user that
                % it is mathematically incorrect and then update the left bound
                % to 2.
                if (fitRange(1) == 1)
                    fitRange(1) = 2;
                    warning('loglog fit can''t operate on delay 0, skip the 1st point now...')
                end
                % Linear regression based on the fit range
                fitRes = polyfit(log10(t(fitRange)),log10(msd(fitRange)),1);
                alpha(ii) = fitRes(1);
            end
            % Pass the result into the instance
            obj.alpha = alpha;
            
            % If any other key was used for specifying the fitting methods, raise the error.
        otherwise
            error('Invalid fitting method. Valid options: linear|loglog.');
    end
end

end