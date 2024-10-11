
function res = get_nuc_num(cellList)
%get extra data for cellList 
 if ~isfield(cellList.meshData{1}{1}, 'length')
        cellList = oufti_makeCellListDouble(cellList);
        for cells = 1:length(cellList.meshData)
            for frames = 1:length(cellList.meshData{cells})
                cellList.meshData{cells}{frames} = getextradata(cellList.meshData{cells}{frames});
         end
        end
    end

% creates res cell output with number of rows = frames, and 1 column
res = cell(length(cellList.meshData),1);
for frame = 1:length(cellList.meshData)
    res{frame,1} = nan(length(cellList.meshData{frame}),4); 
    for cc = 1:length(cellList.meshData{frame})
        res{frame,1}(cc,1) = frame; 
        res{frame,1}(cc,2) = cc; 
        % if a nucleoid is detected, find the number of objects 
        if ~isempty(cellList.meshData{frame}{cc}.objects) 
        res{frame,1}(cc,3) = length(cellList.meshData{frame}{cc}.objects.outlines); 
        else
            continue
        end 
         if ~isempty(cellList.meshData{frame}{cc}.length) 
        res{frame,1}(cc,4) = cellList.meshData{frame}{cc}.length*.065; 
        else
            continue
        end 
    end
end 
 res = vertcat(res{:})   


