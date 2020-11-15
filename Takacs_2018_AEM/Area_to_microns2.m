function cellList = Area_to_microns2(cellList)
%{
-About-
converts the area field of a cellist from pixels^2 to microns^2 given the
conversion unit 0.0642

-Inputs-
cellList: an Oufti cellList with area added

-varargin-
N/A

-Outputs-
cellList: the input cellList with the area field converted to um^2

-Example-
converted_cellList = Area_to_micron2(cellList)
   
-Supplementary-
N/A

-Keywords-
cellList area 

-Dependencies-
N/A

-References-

-Author-
Brad Parry, 2018 June 25
%}

%the number of microns in one pixel
pixel_2_um = 0.0642;

no_conversion_performed = true;

for F = 1:length(cellList.meshData)
    for C = 1:length(cellList.meshData{F})
        if isempty(cellList.meshData{F}{C}) || ~isfield(cellList.meshData{F}{C},'area') || isempty(cellList.meshData{F}{C}.area)
            continue
        end
        
        cellList.meshData{F}{C}.area = cellList.meshData{F}{C}.area * (pixel_2_um^2);
        no_conversion_performed = false;
        
    end
end

if no_conversion_performed
    disp('no area conversion was performed. check that the cellList has the field area')
end

end