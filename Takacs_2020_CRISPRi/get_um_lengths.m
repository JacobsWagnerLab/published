for ii = 1:length(cellList.meshData)
    for jj = 1:length(cellList.meshData{ii})
        %Add missing data fields to the cellList
        cellList.meshData{ii}{jj} = getextradata(cellList.meshData{ii}{jj});
        %Iterate through all cells in the cellList that have a valid field
        %called 'area', since you will use area to normalize
    end
end

count = 0;

for ii = 1:length(cellList.meshData)
    for jj = 1:length(cellList.meshData{ii})
        if isfield(cellList.meshData{ii}{jj},'length')
            count = count + 1;
        end
    end
end

um_length = zeros(count,1);

zz = 0;

for ii = 1:length(cellList.meshData)
    for jj = 1:length(cellList.meshData{ii})
        if isfield(cellList.meshData{ii}{jj},'length')
            zz = zz + 1;
            um_length(zz,:) = cellList.meshData{ii}{jj}.length * 0.0643;
        end
    end
end
        