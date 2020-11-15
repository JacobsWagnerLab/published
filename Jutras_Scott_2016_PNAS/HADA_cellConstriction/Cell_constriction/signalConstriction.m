function [constriction] = signalConstriction(signals,field,sWin)
% calculate the degree of constriction for the signal indicated by field in
% all signals after smoothing with a windows the size of sWin
%
% signals is the output from read_cellLists
% field is the fieldname found in signals and should be a string
% sWin is a smoothing parameter and must be an odd integer
%
% constriction is a vector of constriction for all cells in signals where 0
% indicates complete constriction and 1 indicatios no constriction
%
% Author: Brad Parry

constriction = [];
for C = 1:length(signals)
    S = signals{C}.(field);
    Ssmoothed = conv(S(:)', ones(1, sWin)/sWin,'valid');
    constrictionDiameter = floor((length(Ssmoothed)/4)/2)*2+1;
    [percentConstricted, constrictionLocations, ~] = findValleys(Ssmoothed, constrictionDiameter);
    if isempty(constrictionLocations), constriction(end+1) =1;continue,end
    
    %find constriction nearest cell center
    [~, ix] = min((constrictionLocations - length(Ssmoothed)/2).^2);
    constrictionLocations = constrictionLocations(ix);
    percentConstricted = percentConstricted(ix);
    constriction(end+1) = percentConstricted;

end

function [percentConstricted, constrictionLocations, signal] = findValleys(signal, constrictionDiameter)

signal = signal(:)';

%trim off the pole regions with a derivative technique. This allows the
%cell poles to be of different lengths for each cell and each pole, which
%is much more realistic than specifying a fixed pole length for all cells
d = diff(signal);
firstInd = (d < 0) .* (1: (length(signal)-1) );
firstInd(firstInd==0) = [];
if isempty(firstInd)
    percentConstricted = 1;
    constrictionLocations = nan;
    return
end
firstInd = firstInd(1);

lastInd = (d > 0) .* (2:length(signal));
lastInd(lastInd == 0) = [];
if isempty(lastInd)
    percentConstricted = 1;
    constrictionLocations = nan;
    return
end
lastInd = lastInd(end);

%normalize cell widths by the full cell width which is estimated as the
%median of the 50th percentile widest points that are not part of the cell
%pole.
sortedWidths = sort(signal(firstInd:lastInd));
fullCellWidth = median( sortedWidths(ceil(length(sortedWidths)/2):end) );
signal = signal(firstInd:lastInd) ./ fullCellWidth;

%Routine to find local minima over the diameter defined by
%constrictionDiameter
Rplus1 = ceil(constrictionDiameter/2);
percentConstricted = [];
constrictionLocations = [];
for localCenter = Rplus1 : (length(signal) - Rplus1 + 1)
    inds = [(localCenter-(Rplus1-1)): (localCenter-1), (localCenter+1):(localCenter+(Rplus1-1))];

    if sum( signal(localCenter) < signal(inds) ) == (constrictionDiameter-1)
        percentConstricted(end+1) = signal(localCenter);
        constrictionLocations(end+1) = localCenter;
    end
end
constrictionLocations = constrictionLocations + firstInd - 1;
