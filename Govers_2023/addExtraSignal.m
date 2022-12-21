function [cellList] = addExtraSignal(cellList, cellList2)

%Function to add extra signal data to cellList obtained from Oufti
%CellList: contains info of first 2 fluorescent signals
%CellList2: contains info of subseauent 2 fluorescent signals
%Both cellLists contain same cells, only difference is fluorescent signal
%added in oufti



%% Clean out cell entries ahead of time to avoid 'if' testing in parallel loop
kk=0;
for ii=1:length(cellList.meshData)
    for jj=1:length(cellList.meshData{ii})
        if isempty(cellList.meshData{ii}{jj}) ||...
            ~isfield(cellList.meshData{ii}{jj},'mesh') ||...
            length(cellList.meshData{ii}{jj}.mesh)<=4
            cellList.meshData{ii}{jj}=[];
            kk=kk+1;
        end
    end
end
disp(['Number of cell entries cleaned out : ',num2str(kk)]);

%%Add extra signal data
w = waitbar(0,'Adding extra fluorescent signals...');
clear ii jj

for ii=1:length(cellList.meshData)
    clear meshData tmpCell
    
    if ~isempty(cellList.meshData{ii})
        for jj=1:length(cellList.meshData{ii})
            if ~isempty(cellList.meshData{ii}{jj})
                if isfield(cellList2.cellList.meshData{ii}{jj}, 'signal1')
                cellList.meshData{ii}{jj}.signal3= cellList2.cellList.meshData{ii}{jj}.signal1;
                end
                if isfield(cellList2.cellList.meshData{ii}{jj},'signal2')
                cellList.meshData{ii}{jj}.signal4= cellList2.cellList.meshData{ii}{jj}.signal2;
                end
            end
        end
    end 
    waitbar(ii/length(cellList.meshData));
%     disp(['Frame ',num2str(ii),' : ',num2str(toc)]);
end
close(w);
clear meshData tmpCell ii jj 