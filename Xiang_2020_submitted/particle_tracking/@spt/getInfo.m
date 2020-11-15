%% Print the basic information (private attributes) about the instance.
%
% -About-
%   When the constructor creates the instance, the input from user (e.g.
%   stackPath, particleSize et. al) is stored privately to prevent
%   accidental changes to the object. This privately stored fields can be
%   check using getInfo(). It prints all the privately stored fields. To
%   modify any of the field, the user make invoke setInfo() with the
%   specific filed name as a key.
% 
% -Inputs-
%   - obj: spt object
% 
% -Varargout-
%   Prints to the command window about all the private attributes of the
%   instance, if no output requested.
%   If at least 1 output requested, the result is the structure containing
%   all fields stored in the info.
% 
% -Example-
%   % Print out the information about the experiment/object
%   myParticle.getInfo();
% 
% -Author-
%   Yingjie Xiang, CJW Lab, Yale University


function varargout = getInfo(obj)
% Create a struct to assemble all the private attributes
info = struct('date',        obj.date,...
              'about',       obj.about,...
              'stackPath',   obj.stackPath,...
              'particleSize',[num2str(obj.particleSize*1e3), ' nm'],...
              'pixelLength', [num2str(obj.pixelLength*1e3),' nm'],...
              'frameTime',   [num2str(obj.frameTime*1e3),' msec'],...
              'imNum',       obj.imNum,... 
              'iThreshold',  obj.iThreshold,...
              'minSpotSize', [num2str(obj.minSpotSize),' px'],...
              'maxSpotSize', [num2str(obj.maxSpotSize),' px'],...
              'eccentricity', obj.eccentricity);
% Display in the command window
if (nargout == 0)
    disp(info);
else
    varargout{1} = info;
end