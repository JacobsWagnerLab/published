function [pos] = SpotDat2Pos(SpotDat)

%{
-About-
Convert SpotDat output from particleTracking.m into the pos array suitable
for localization linking (pos2trajectories.m)

-Inputs-
SpotDat: output from particleTracking.m

-varargin-
na

-Outputs-

pos: a structured cell array: 
  pos{cellN}.x = a matrix with number of columns == number of frames and
  each row is a localization event. It is OK to have missing
  localizations (which should be indicated with a 0. 


-Example-
   
-Supplementary-

-Dependencies-


-References-
n/a

-Author-
Brad Parry, 2014 July 28
%}

% convert SpotDat into a cell array of positions over time 
pos = cell(1,length(SpotDat{1}));

if iscell(SpotDat)
    %its a cell of cell arrays...

    for F = 1:length(SpotDat)
        for C = 1:length(SpotDat{F})
            if isempty(SpotDat{F}{C}) || ~isfield(SpotDat{F}{C},'Centroid'), continue, end
    %         add the shifts  back to the localizations
            if isfield(SpotDat{F}{C},'shift')
                x = SpotDat{F}{C}.Centroid(:,1) + SpotDat{F}{C}.shift(1);
                y = SpotDat{F}{C}.Centroid(:,2) + SpotDat{F}{C}.shift(2);
            else
                x = SpotDat{F}{C}.Centroid(:,1);
                y = SpotDat{F}{C}.Centroid(:,2);
            end
            pos{C}.x(1:length(x),F) = x(:)';
            pos{C}.y(1:length(y),F) = y(:)';
        end
    end

elseif istruct(SpotDat) && isfield(SpotDat,'meshData')

    
end