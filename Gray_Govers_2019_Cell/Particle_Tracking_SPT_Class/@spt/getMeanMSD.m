%% Calculate the ensemble-averaged MSD curve weighted by the track length.
%
% -About-
%   The method getMeanMSD() calculates the ensemble-averaged MSD based on
%   the track length. That is, the ensemble-average was weighted based on
%   the length of individual MSD curves. The ensemble-averaged MSD
%   represents a population average of all particle MSDs.
%   
%   The number of delays for a specific time delay is first extracted from
%   the MSD (the thrid column). All the delays are then consolidated
%   together. The total number of delays for a specific time delay is then
%   calculated. The weight is determined by dividing the number of delays a
%   particle has over the total number of delays (at that specific time
%   delay). The ensemble-averaged MSD is calculated as a dot product
%   between the MSD and its weight (again at that specific time delay).
% 
% -Input-
%   - obj: spt object
% 
% -Output-
%   - meanMSD: the ensemble-averaged mean squared displacement. A Nx2
%              matrix, where N is the number of time delays. The first
%              column represents the time delays, and the second column
%              represents the ensemble-average. The units are not converted
%              to physical units [frame, px^2].
% 
% -Example-
%   % Calculate the ensemble-averaged MSD
%   myParticle.getMeanMSD();
% 
% -Author-
%   Yingjie Xiang, CJW Lab, Yale Univeristy


function getMeanMSD(obj)

% Check if the MSD were calculated previously. If not, calculate now.
if isempty(obj.MSD)
   obj.getMSD();
end

% Number of trajectories/particles
np = size(obj.MSD,1); 
% Memory ength of each trajectory
len = size(obj.MSD,1); 
% Calculate the length of each MSD
for n = 1 : np
    len(n) = size(obj.MSD{n},1);
end
% Find the longest MSD
max_len = max(len);
% Each col is a trajectory, each row is a delay
allDelays = zeros(max_len,np); 

% Consolidate all delay counts
for i = 1 : np
    allDelays(1:len(i),i) = obj.MSD{i}(:,3);
end

% Find the weight for each delay point based on the number of delays
% Memory for the weight matrix
w = zeros(max_len,np);

% Total number of delays at a specific time delay
% For example, for a 1000-step long track, its number of delays with a
% 20-step delay is 980. While for a 25-step long track, its number of
% delays is only 5. At larger delay, the ensemble-average is dominated by
% only the long tracks. We are essentially calculating the dot product
% between the MSD and its number of delays at a certain delay.
sumDelays = sum(allDelays,2);

% Loop through individual MSD and find weights
for j = 1 : size(allDelays,1)
    w(j,:) = allDelays(j,:)./ sumDelays(j);
end

% Calculate mean MSD based on the weights above
% Memory for storing the ensemble-averaged MSD
meanMSD = zeros(max_len,2);

% Set the time delay column the same as the longest MSD
meanMSD(:,1) = obj.MSD{find(len==max_len,1,'first')}(:,1);

% Loop through every particle to append its contribution
% This is essentially a dot product.
for k = 1 : np
    % Extract a single MSD
    thisMSD = obj.MSD{k}(:,2);
    % Append its contribution based on its weight
    meanMSD(1:len(k),2) = meanMSD(1:len(k),2) + thisMSD.*w(1:len(k),k);
end

% Pass the final result back to the object
obj.meanMSD = meanMSD;
end