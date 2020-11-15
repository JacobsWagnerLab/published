%% Remove short trajectories
%
% -About-
%   The method rmShort() removes the short trajectories and re-calculates
%   the MSD and mean MSD using the updated longer trajectories. To avoid
%   any error, any fitted diffusion coefficients or alphas would be
%   deleted.
% 
% -Input-
%   - obj: spt object
%   - minLen: minimal length of the trajectories to keep
% 
% -Output-
%   Longer trajectories with length at least (including) the input minimum
%   threshold.
%
% -Example-
%   % Remove all trajectories shorter than 10 steps long
%   myParticle.rmShort(10);
% 
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University


function rmShort(obj,th)
    obj.getTrackLen();
    obj.tracks(obj.trackLen < th) = [];
    
    % Update and clean
    obj.getTrackLen();
    obj.getMSD();
    obj.getMeanMSD();
    obj.D = [];
    obj.alpha = [];
end