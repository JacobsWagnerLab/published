%% Set private member data of the instance.
%
% -About-
%   The method setInfo() updates the private data member of the object with
%   the user input. The user specifies the member to be modified with a
%   key, followed by the updates. If the input key was not a member of the
%   object, the user get notified by an error.
%
% -Input-
%   - obj: spt object
%
% -Varargin-
%   - member key as a string followed by the value for that key
%
% -Example-
%   % Update the date
%   myParticle.setInfo('date','July-01-2017');
%   % Update the about and particle size
%   myParticle.setInfo('date','June-01-2017','particleSize',0.1);
%
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University

function setInfo(obj,varargin)

keys = {'date','about','stackPath','particleSize','pixelLength',...
        'frameTime','imNum','iThreshold','minSpotSize','maxSpotSize',...
        'eccentricity'};

% Check if the user input is a member of the private data members. If yes,
% updates the member with the user input. Else, raise an error.

for ii = 1:2:length(varargin)
    key = varargin{ii};
    val = varargin{ii+1};
    if any(ismember(lower(keys),key))
        switch lower(key)
            case 'date'
                obj.date = val;
            case 'about'
                obj.about = val;
            case 'stackpath'
                obj.stackPath = val;
            case 'particlesize'
                obj.particleSize = val;
            case 'pixellength'
                obj.pixelLength = val;
            case 'frametime'
                obj.frameTime = val;
            case 'imnum'
                obj.imNum = val;
            case 'ithreshold'
                obj.iThreshold = val;
            case 'minspotsize'
                obj.minSpotSize = val;
            case 'maxspotsize'
                obj.maxSpotSize = val;
            case 'eccentricity'
                obj.eccentricity = val;
        end
    else
        error('No field associates with the input: %s.',key);
    end
end
end