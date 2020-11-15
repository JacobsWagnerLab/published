function signals = reformatCellList(cellList, signalNames)
%this is a helper function to reformat the cellList into a data format that
%is more amenable to signal analysis.reformatCellList finds all valid cells
%that also have complete signals and puts them into a linear cell area with
%their signals normalized by mesh step area.
%
% cellList.meshData{1}{1}    ----->   signals{1}.*
%
%signals will have fields that correspond to the signalNames provided.
%signalNames MUST but the same length as the number of signals in the
%cellList and MUST match their order. If not, this code will flip signals
%around....if signalNames are not provided, the behavior will be:
%       cellList.meshData{Frame}{CellInd}.signal0 ---> signals{N}.phs
%       cellList.meshData{Frame}{CellInd}.signal1 ---> signals{N}.gfp
%       cellList.meshData{Frame}{CellInd}.signal2 ---> signals{N}.dna
%       cellList.meshData{Frame}{CellInd}.signal3 ---> signals{N}.hada
%
%in addition to adjusting the signal data, reformatCellList will add
%cellLength and frame and cell index information to the signals. the field
%index is a linear index and is NOT the cell ID from cellList.cellId
%
%
%
%Brad Parry, Christine Jacobs-Wagner lab; 2016 April

if nargin == 1
    signalNames = {'phs','gfp','dna','hada'};
end


signals = [];
for F = 1:length(cellList.meshData)
    for C = 1:length(cellList.meshData{F})
        if isempty(cellList.meshData{F}{C}) || ~isfield(cellList.meshData{F}{C},'mesh') || length(cellList.meshData{F}{C}.mesh) < 6
            continue
        end
        
        if length(cellList.meshData{F}{C}.signal0) ~= (length(cellList.meshData{F}{C}.mesh(:,1))-1)
            continue
        end
        
        badSignal = false;
        for q = 1:length(signalNames)
            fn = ['signal',num2str(q-1)];
            if isempty(cellList.meshData{F}{C}.(fn)) || sum(cellList.meshData{F}{C}.(fn)) == 0
                badSignal = true;
            end
        end
        
        if badSignal
            continue
        end
        
        signals{end+1} = [];
        signals{end}.index = [F,C];
        
        x = mean(cellList.meshData{F}{C}.mesh(:,[1,3]),2);
        y = mean(cellList.meshData{F}{C}.mesh(:,[2,4]),2);
        signals{end}.cellLength = sum((diff(x).^2 + diff(y).^2).^(1/2));

        segmentAreas = zeros(size(cellList.meshData{F}{C}.mesh,1)-1,1);
        for q = 2:size(cellList.meshData{F}{C}.mesh,1)
            x = [cellList.meshData{F}{C}.mesh(q-1,1), cellList.meshData{F}{C}.mesh(q-1,3),...
                cellList.meshData{F}{C}.mesh(q,3), cellList.meshData{F}{C}.mesh(q,1),...
                cellList.meshData{F}{C}.mesh(q-1,1)];
            y = [cellList.meshData{F}{C}.mesh(q-1,2), cellList.meshData{F}{C}.mesh(q-1,4),...
                cellList.meshData{F}{C}.mesh(q,4), cellList.meshData{F}{C}.mesh(q,2),...
                cellList.meshData{F}{C}.mesh(q-1,2)];
            segmentAreas(q-1) = polyarea(x,y);
        end

        for q = 1:length(signalNames)
            fn = ['signal',num2str(q-1)];
            signals{end}.(signalNames{q}) = cellList.meshData{F}{C}.(fn) ./ segmentAreas;
        end

    end
end