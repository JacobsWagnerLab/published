%% Calculate the radius of gyration of each track.
%
% -About-
%   The method getRg() calculates the radius of gyration for all tracks
%   stored in the object. The radius of gyration is calculated based on the
%   equation described in Ref. Parry (2016), on the supplemental S5. The
%   mean positions in both x and y dimensions are first calculated. The
%   radius gyration is the square root of the normalized squared distance
%   between the mean position and individual positions.
% 
% -Input-
%   - obj: spt object
% 
% -Output-
%   - Rg: radius of gyration for all tracks. Nx1 matrix, where N is the
%         number of tracks. The unit has been converted to physical units
%         [µm].
% 
% -References-
%   Parry et. al, Cell 156, 183-194 (2014).
% 
% -Example-
%   % Calculate the radius of gyration
%   myParticle.getRg()
% 
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University


function getRg(obj)

% Memory for storing the radius of gyration
obj.Rg = zeros(length(obj.tracks),1);

% Loop through individual tracks
for ii = 1:length(obj.tracks)
    
    % Extract a sigle track
    thisTrack = obj.tracks{ii};
    
    % Convert to physical units
    xs = thisTrack(:,2) .* obj.pixelLength;
    ys = thisTrack(:,3) .* obj.pixelLength;
    
    % Calculate the mean positions
    x_mean = mean(xs);
    y_mean = mean(ys);
    
    % Calculate the radius of gyration
    obj.Rg(ii) = sqrt(sum((xs - x_mean).^2 + (ys- y_mean).^2)./size(thisTrack,1));
end

end