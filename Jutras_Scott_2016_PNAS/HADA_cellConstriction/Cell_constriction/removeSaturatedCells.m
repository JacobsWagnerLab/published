function [cellList,removedInds] = removeSaturatedCells(cellList, ims)
% remove cells which have saturated pixels. Hardcoded to identify saturated
% pixels in 12 and 16 bit images only.
%
% cellList: an Oufti returned cellList
% ims: cell array of 3d matrices. Each matrix should be
% images(imRows,imCols,imFrames) so that each cell array of ims can be a
% channel of acquired images from a complete timelapse.
%
% the size(ims{N},3) must be equal to the length of cellList.meshData for
% all N
%
% Return a cellList with saturated cells removed and indices of removed
% cells.
% 
% Author: Brad Parry
removedInds = [];


saturatingValues = [2^16, 2^12]-1;
Frames2check = [];
for S = 1:length(ims)
    for F = 1:size(ims{S},3)
        
        if sum( max( max(ims{S}(:,:,F)) ) == saturatingValues ) > 0
            Frames2check(end+1) = F;
        end
        
    end
end
Frames2check = unique(Frames2check);
%check which frames have saturated pixels

bxlim = [1 size(ims{1},1) 1 size(ims{1},2)];
for F = Frames2check
    
    for C = length(cellList.meshData{F}):-1:1
        if isempty(cellList.meshData{F}{C}) || ~isfield(cellList.meshData{F}{C},'mesh') || length(cellList.meshData{F}{C}.mesh) < 6
            cellList.meshData{F}{C} = [];
            cellList.cellId{F}(C) = [];
            continue
        end
        
        cellMesh = double( cat(1, cellList.meshData{F}{C}.mesh(:,1:2), flipud(cellList.meshData{F}{C}.mesh(:,3:4))) );
        %find the smallest bounding box for the cell
        bx = [floor(min(cellMesh(:,2))), ceil(max(cellMesh(:,2))), floor(min(cellMesh(:,1))), ceil(max(cellMesh(:,1)))];
        bx([1 3]) = bx([1 3]); 
        bx([2 4]) = bx([2 4]);
        %if the box is outside of the image, bring it back
        ck = (bx - bxlim).*[1 -1 1 -1];
        bx(ck < 0) = bxlim(ck < 0);

        isSaturated = false;
        for q = 1:length(ims)
            maxVal = max( max( ims{q}(bx(1):bx(2),bx(3):bx(4),F) ) );
            isSaturated(sum(maxVal == saturatingValues) > 0) = true;
        end

        if isSaturated
           cellList.meshData{F}(C) = [];
           cellList.cellId{F}(C) = [];
           removedInds(end+1,:) = [F,C];
           continue
        end

        
    end
end