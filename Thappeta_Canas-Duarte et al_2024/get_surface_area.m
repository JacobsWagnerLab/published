%% Estimating the surface area of cells defined in the Oufti mesh
% -Input-
%   cellList: Oufti cellList
% -Output-
%   cellList: Oufti cellList with two additional fields in each cell
%   structure:
%       1. step_surface_area: surface area of the cone defined by the
%       adjacent "ribs" of the cell mesh
%       2. surface_area: the sum of step_surface_area, that is, the total
%       surface area of this cell
% -Author-
%   Yingjie Xiang, 12/30/2020

function cellList = get_surface_area(cellList)
for frame = 1:length(cellList.meshData)
    for cc = 1:length(cellList.meshData{frame})
        mesh = cellList.meshData{frame}{cc}.mesh;
        % The length of each rib in the cell mesh
        rib_length = sqrt((mesh(:,3)-mesh(:,1)).^2+(mesh(:,4)-mesh(:,2)).^2);
        % The distance between two adjacent ribs
        ls = cellList.meshData{frame}{cc}.steplength;
        % Calculating the surface area of the cone defined by the adjacent
        % ribs, using the following equation:
        % SA = pi*(R+r)*sqrt(h^2+(R-r)^2), where R and r are the lengths of
        % the adjacent ribs, h is the distance between them, and SA is the
        % lateral surface area of the cone, i.e., without adding the areas
        % from the top and bottom circles
        r = rib_length(1:end-1);
        R = rib_length(2:end);
        cellList.meshData{frame}{cc}.step_surface_area = pi.*(R+r).*sqrt(ls.^2+(R-r).^2);
        cellList.meshData{frame}{cc}.surface_area = sum(cellList.meshData{frame}{cc}.step_surface_area);
    end
end
end