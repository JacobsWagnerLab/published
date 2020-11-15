%% Find the y-intercepts of all MSD curves.
%
% -About-
%   The diffusion coefficient may also be deduced from the y-intercept of
%   the loglog MSD curves, since MSD=4Dt (2D tracking), log(MSD) = log(4D)
%   + log(t). D = 10^(intercepts)/4. However, in real experiments, we never
%   measure the intercepts (requires infinitely short time delay). Instead,
%   we can either fit the loglog MSD curves using data at multiple delays
%   or simply assume the diffusion is Brownian and extrapolate just using
%   the point at the shortest delay. Here, we implement the latter. The
%   intercepts are deduced in assumption of the normal diffusion at small
%   delay. This assumption is often true for all practical purposes.
%
% -Inputs-
%   - obj: spt object
%   - delay: optional. Default to the smallest delay, obj.frameTime
%
% -Output-
%	- Intercepts: the y-intercepts of the MSD curves.
%                 (Note: the y-intercepts are in the log10 scale.)
%   If delay == 1, the y-intercepts returned are log(4D)
%       In order to calculate the diffusion coefficient: D = 10^(y-intercepts) / 4.
%       The y-intercepts calculated are the MSD at t = 1 sec, such that
%       logMSD = y-intercepts + log(1). The y-intercepts are equivalent to log(4D)).
%   If delay == obj.frameTime
%       y-intercepts at the smallest time delay will be returned in the log
%       scale.
%
% -Example-
%   % Calculate the y-intercepts of the loglog MSD at the smallest delay
%   myParticle.getIntercept();
%   % Deduce the diffusion coefficients from the intercepts
%   D = 10.^myParticle.getIntercept(1)./4;
%
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

function getIntercept(obj,delay)

if nargin < 2
    delay = obj.frameTime;
end

if delay ~= 1 && delay ~= obj.frameTime
    fprintf('Input argument delay must be either 1 or the smallest delay \n');
    fprintf('Now correcting the delay to be the smallest delay %.4f sec \n', obj.frameTime);
    delay = obj.frameTime;
end

% Memory for collecting the MSD at the smallest delay
x = zeros(length(obj.MSD),1);
% Loop through all MSD curves
for ii = 1 : length(obj.MSD)
    % Extract the MSD at the smallest delay
    x(ii) = obj.MSD{ii}(2,2).*obj.pixelLength^2;
end

if delay == 1
    % Assuming normal diffusion (loglog MSD slope 1), extrapolate the
    % y-intercepts and pass to the object
    obj.intercepts = log10(x) - log10(obj.frameTime) * 1;
else
    obj.intercepts = log10(x);
end

end