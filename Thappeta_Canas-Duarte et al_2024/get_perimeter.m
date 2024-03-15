%% Estimating the perimeter of cells defined in the Oufti mesh
% -Input-
%   cellList: Oufti cellList
% -Output-
%   cellList: Oufti cellList with the additional field
%       1. perimeter: the perimeter of the cell contour
% -Author-
%   Yingjie Xiang, 12/30/2020

function cellList = get_perimeter(cellList)
for frame = 1:length(cellList.meshData)
    for cc = 1:length(cellList.meshData{frame})
        % To calculate the perimeter, just add up the length of the sides
        % defining the polygon of the cell mesh. The "mesh" field has 4
        % columns, which correspond to the x and y coordinates of the
        % halves of the cell mesh
        mesh = cellList.meshData{frame}{cc}.mesh;
        left = sqrt((mesh(2:end,1)-mesh(1:end-1,1)).^2+(mesh(2:end,2)-mesh(1:end-1,2)).^2);
        right = sqrt((mesh(2:end,3)-mesh(1:end-1,3)).^2+(mesh(2:end,4)-mesh(1:end-1,4)).^2);
        cellList.meshData{frame}{cc}.perimeter = sum(left) + sum(right);
    end
end
end