function [cell_properties] = Extract_Cell_Properties(cellList,objects_name, varargin)
%{
-About-
Function to extract properties of individual cells from cellLists generated
by Oufti.

-Inputs-
cellList:  cellList generated by Oufti, loaded into matlab. 
            If SeqA (for E.coli cells) or DnaN (for C. crescentus) 
            information needs to be added, run separate function 
            (Add_SeqA_Area or Add_DnaN_Area, respectively) first. The
            function will automatically look for this information.
objects_name:   name used to store information for detected objects (using
                Oufti's objectDetection module). Typically this is 'object'
                or 'nucleoiddata'. Enter 'none' in case no objects were
                detected
 
-varargin-
    'pixel_size':   pixel size, image resolution. If not specified, default
                    pixel size of 0.64 micron per pixel will be used

-Outputs-
cell_properties:    array in which each row corresponds to a cell and each
                    column to a different cellular property. The column
                    order is: 
                    1: frame
                    2: cell number
                    3: cell length (micron)
                    4: cell width (micron)
                    5: cell area (squared micron)
                    6: cell volume (cubic micron)
                    7: surface area (squared micron)
                    8: surface area over volume (micron-1)
                    9: nucleoid area (squared micron)
                    10: NC ratio
                    11: relative SeqA area
                    12: relative DnaN area

-Example-
NA

-Supplementary-
NA

-Keywords-
Cell properties, cell features, morphology, nucleoid, cellList, Oufti

-Dependencies-
edist (function to calculate euclidean distance between points in 2D)

-References-
NA

-Author-
Sander Govers, 22 October 2018
%}

%% Parse varargin input
tic
%if pixel size is not specified, use default pixel size of 0.064
%micron per pixel.
pixel_size=0.064;

%Examine whether pixel size was specified using varargin input
for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'pixel_size')
        %if pixel size is specified, use this to extract cellular
        %information
        pixel_size = varargin{k+1};
        %Check if specified pixel size is a single number, throw error if
        %this is not the case
        if length(pixel_size) ~= 1
            error('argument for pixel size was incorrectly formed')        
        end
    end
end

%% Clean out cell entries in cellList that do not contain correct information

%Initiate counter to keep track of the number of removed cells
counter=0;

%Loop through frames and cells in the cellList
for ii=1:length(cellList.meshData)
    for jj=1:length(cellList.meshData{ii})
        %Remove empty entries, entries that do not contain cell mesh
        %information or very short entries
        if isempty(cellList.meshData{ii}{jj}) ||...
                ~isfield(cellList.meshData{ii}{jj},'mesh') ||...
                length(cellList.meshData{ii}{jj}.mesh)<=4
            cellList.meshData{ii}{jj}=[];
            %Add 1 to counter if cell entry is removed
            counter=counter+1;
        end
    end
end
%Display the total number of cell entries that was removed from the
%cellList
disp(['Number of cell entries cleaned out : ',num2str(counter)]);


%% Extract properties from individual cells

w = waitbar(0,'Data coming soon...');

%Generate return structure
tab = struct('frame',[],'cell_id',[],'cell_length',[],'cell_width',[],'cell_area',[],...
    'cell_volume',[],'cell_surface_area',[],'cell_surface_area_over_volume',[],...
    'nucleoid_area',[],'NC_ratio',[],'relative_seqA_area',[],'relative_dnaN_area',[]);

%Loop through frames and cells in the cellList and extract all features
for ii=1:length(cellList.meshData)
    clear meshData tmpCell tmp_structure
    %Only do this for frames containing cells
    if ~isempty(cellList.meshData{ii})      
        %for shorter notation
        mesh_data = cellList.meshData{ii};
        %Generate temporary structure to store data of current frame
        tmp_structure = struct('frame',[],'cell_id',[],'cell_length',[],'cell_width',[],'cell_area',[],...
        'cell_volume',[],'cell_surface_area',[],'cell_surface_area_over_volume',[],...
        'nucleoid_area',[],'NC_ratio',[],'relative_seqA_area',[],'relative_dnaN_area',[]);
        %Loop through all cells in the current frame
        for jj=1:length(mesh_data)
            if ~isempty(mesh_data{jj})
                % Extract frame number and cell ids
                tmp_structure(jj).frame = ii;
                tmp_structure(jj).cell_id = cellList.cellId{ii}(jj);
            
                % Extract cell length from cell mesh information
                cell_mesh=double(mesh_data{jj}.mesh);

                %Calculate the distance between each cell segment
                step_length = edist(cell_mesh(2:end,1)+cell_mesh(2:end,3),cell_mesh(2:end,2)+cell_mesh(2:end,4),...
                    cell_mesh(1:end-1,1)+cell_mesh(1:end-1,3),cell_mesh(1:end-1,2)+cell_mesh(1:end-1,4))/2;
                %Use the sum of distances to calculate total cell length
                tmp_structure(jj).cell_length = sum(step_length)*pixel_size;

                % Cell width
                x1=cell_mesh(:,1);
                y1=cell_mesh(:,2);
                x2=cell_mesh(:,3);
                y2=cell_mesh(:,4);
                width = sort(sqrt((x1-x2).^2+(y1-y2).^2),'descend');
                %Calculate cell width based on the mean cell width of the
                %30% highest cell width measured across the cell contour
                tmp_structure(jj).cell_width = mean(width(1:floor(length(width)/3)))*pixel_size;

                % Cell area
                mesh_length = size(cell_mesh,1)-1;
                step_area=zeros(mesh_length,1);
                for counter=1:mesh_length
                    %Calculate the area of each cell segment
                    step_area(counter)=polyarea([cell_mesh(counter:counter+1,1);cell_mesh(counter+1:-1:counter,3)],...
                        [cell_mesh(counter:counter+1,2);cell_mesh(counter+1:-1:counter,4)]);
                end
                %Use the sum of areas to calculate total cell area
                tmp_structure(jj).cell_area = sum(step_area)*pixel_size*pixel_size;
                
                % Cell volume
                d = edist(cell_mesh(:,1),cell_mesh(:,2),cell_mesh(:,3),cell_mesh(:,4));
                %Calculate the volume of each cell segment
                step_volume = (d(1:end-1).*d(2:end) + (d(1:end-1)-d(2:end)).^2/3).*step_length*pi/4;
                %Use the sum of volumes to calculate total cell volume
                tmp_structure(jj).cell_volume = sum(step_volume)*pixel_size*pixel_size*pixel_size;
                
                % Surface area
                step_surface_area=zeros(length(cell_mesh)-1,1);
                for counter=1:length(cell_mesh)-1
                    %Calculate the surface area of each cell segment
                    step_surface_area(counter)=pdist2(cell_mesh(counter,1:2),cell_mesh(counter+1,1:2));
                end
                %Use the sum of surface area to calculate total cell
                %surface area
                tmp_structure(jj).cell_surface_area = sum(pi*step_surface_area.*(width(1:end-1)+width(2:end))./2)*pixel_size*pixel_size;
                
                % Surface area to volume ratio
                tmp_structure(jj).cell_surface_area_over_volume = tmp_structure(jj).cell_surface_area/tmp_structure(jj).cell_volume;
          
                
                % Extract nucleoid chracteristics
                % Set default values to NaN, in case no nucleoid were
                % detected or nucleoid detection failed for a certain cell
                tmp_structure(jj).nucleoid_area = NaN;
                tmp_structure(jj).NC_ratio = NaN;

                %Check if nucleoid information is present in cell
                %properties, this is done using the objects name specified
                %by the user in the input of the function
                if isfield(mesh_data{jj},objects_name)
                   tmp_nucleoid_area=[];
                   %Loop through individual nucleoids within a cell and
                   %store their area
                    if ~isempty(mesh_data{jj}.(objects_name).outlines)
                        for piece=1:length(mesh_data{jj}.(objects_name).outlines)
                            tmp_nucleoid_area=[tmp_nucleoid_area;double(mesh_data{jj}.(objects_name).area{piece})];
                        end
                        %Sum area of all segments to extract nucleoid area
                        tmp_structure(jj).nucleoid_area = nansum(tmp_nucleoid_area)*pixel_size*pixel_size;
                        %Divide nucleoid area by the cell area to obtain the NC
                        %ratio
                        tmp_structure(jj).NC_ratio = tmp_structure(jj).nucleoid_area/tmp_structure(jj).cell_area;
                    end
                end
                
                %Extract relative SeqA area is this field is present in the
                %cellList
                tmp_structure(jj).relative_seqA_area=NaN;
                if isfield(mesh_data{jj},'relative_seqA_area')
                    tmp_structure(jj).relative_seqA_area=mesh_data{jj}.relative_seqA_area;
                end
                %Extract relative DnaN area is this field is present in the
                %cellList
                tmp_structure(jj).relative_dnaN_area=NaN;
                if isfield(mesh_data{jj},'relative_dnaN_area')
                    tmp_structure(jj).relative_dnaN_area=mesh_data{jj}.relative_dnaN_area;
                end
            end
        end
        tab=[tab,tmp_structure];
    end
    
    waitbar(ii/length(cellList.meshData));

end
tab(1)=[];
toc
close(w);
clear ii jj kk steplength steparea stepvolume width d x1 x2 y1 y2 a b perDist
clear mBox imCell X Y contInt mesh lng w tmpStr piece perim Npixelvals
clear img imdb nucleoidarea nuCirc cell cellEntries meshData
out = tab;

%Convert structure to array
tmp_array(:,1)=cat(1,tab.frame);
tmp_array(:,2)=cat(1,tab.cell_id);
tmp_array(:,3)=cat(1,tab.cell_length);
tmp_array(:,4)=cat(1,tab.cell_width);
tmp_array(:,5)=cat(1,tab.cell_area);
tmp_array(:,6)=cat(1,tab.cell_volume);
tmp_array(:,7)=cat(1,tab.cell_surface_area);
tmp_array(:,8)=cat(1,tab.cell_surface_area_over_volume);
tmp_array(:,9)=cat(1,tab.nucleoid_area);
tmp_array(:,10)=cat(1,tab.NC_ratio);
tmp_array(:,11)=cat(1,tab.relative_seqA_area);
tmp_array(:,12)=cat(1,tab.relative_dnaN_area);

%Write output array
cell_properties=tmp_array;

end

function d=edist(x1,y1,x2,y2)
%{
-About-
Function that caculates euclidean distance between 2 points

-Inputs-
x1: x-coordinate of first point
y1: y-coordinate of first point
x2: x-coordinate of second point
y2: y-coordinate of second point

-Author-
Sander Govers
%}    
    d=sqrt((x2-x1).^2+(y2-y1).^2);
end