%{
-About-
Calculates the fluorescence of each cell, normalized by the area of that
cell
 
-Inputs-
cellList: an Oufti cellList 
 
-varargin-
N/A
 
-Outputs-
normalized_fluorescence: a matrix in which each entry corresponds to the
fluorescence (A.U.) of a single cell, normalized by the area of that cell
 
-Example-
normalized_fluorescence = CalculateFluorPerCell(cellList)
   
-Supplementary-
N/A
 
-Keywords-
cellList area fluorescence 
 
-Dependencies-
getextradata.m
 
-References-
 
-Author-
Molly Scott, 2018 July 08
%}


function normalized_fluorescence = CalculateFluorPerCell(cellList)

%Initialize a cell array to hold your normalized fluorescence data
normalized_fluorescence = [];

%Iterate through the cellList
for ii = 1:length(cellList.meshData)
    for jj = 1:length(cellList.meshData{ii})
        %Add missing data fields to the cellList
        cellList.meshData{ii}{jj} = getextradata(cellList.meshData{ii}{jj});
        %Iterate through all cells in the cellList that have a valid field
        %called 'area', since you will use area to normalize
        if isfield(cellList.meshData{ii}{jj},'area')
            %Calculate the fluorescence per unit area for each cell
            fluor_per_unit_area = sum(cellList.meshData{ii}{jj}.signal1) / cellList.meshData{ii}{jj}.area;
            %Add your calculated value to the growing list of normalized
            %fluorescence values
            normalized_fluorescence{end+1} = fluor_per_unit_area;
        end
    end
end

%Convert the cell array to a matrix for ease of further calculations 
normalized_fluorescence = cell2mat(normalized_fluorescence);
end