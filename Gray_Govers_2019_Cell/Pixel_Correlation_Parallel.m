function [allCorrs] = Pixel_Correlation_Parallel(imPaths, cellList_file, dx_from_center, pole_length)
%{
-About-
This function performs pixel correlation on two signals of interest in a
cellList and returns the correlations.


-Inputs-
imPaths = 1 x 2 cell array that contains the images that are to be compared
for correlations (along with their directory). Check
Pixel_Correlation_Multiple_Experiments_Scan.m for format.

cellList_file = location of the Oufti cellList that you want to analyze.

dx_from_center = 1;
This determines how many pixels away from the centerline you want to
calculate the correlation. Must add up to an integer. A value of one
means that the correlation area will be two pixels wide

pole_length = 8;
This is the number of pixels away from the pole that you want to
calculate the correlation for. This is to avoid the artificial positive
correlation at the pole that results from the decrease in volume.

-Outputs-
This function performs pixel correlation on two signals of interest in a
cellList and returns the correlations (see MATLAB corr).


-Keywords-
Pixel correlation
Ribosome segregation

-Dependencies-
getTIFsFromPaths.m
Cell_Pixel_Correlation.m
Extract_Cell_Pixels.m
Cell_Projection.m
Taylor_Smooth.m

-Author-
William Gray October 10, 2018
%}

%load in the cellList
cellListData = load(cellList_file);

% the imPaths directories must contain images from distinct channels
% (i.e., only GFP or only DAPI, etc). See documentation in getTIFsFromPaths
% for more information
files = getTIFsFromPaths(imPaths);

% The following loop measures image correlation for biological cells in the
% cellList. This loop assumes the data is a time course and files{1}{1}
% corresponds to cellList.meshData{1} and files{1}{2} corresponds to
% cellList.meshData{2}. Remember, the first index of files indicates the
% channel and is dictated by the order of imPaths (above). See
% documentation in getTIFsFromPaths and cellPxCorr for more information.
allCorrs = cell(1,length(cellListData.cellList.meshData));
cellIdList = cell(1,length(cellListData.cellList.meshData));

%Perform the correlation analysis in parallel to speed up computation
parfor Frame = 1:length(cellListData.cellList.meshData)
    allCorrs{Frame} = NaN(1,length(cellListData.cellList.cellId{Frame}));
    
    images = cell(1,length(imPaths));
    for Channels = 1:length(imPaths)
        %Read in all channels for the given frame of interest.
        images{Channels} = double(imread(files{Channels}{Frame}));
    end
    
    counter = 1;
    for Cell = cellListData.cellList.cellId{Frame}
        %Compute pixel correlations across all images for the indicated
        %biological cell. See cellPxCorr for more information.
        try
            [corrMatrix, ~] = Cell_Pixel_Correlation(images, cellListData.cellList, Frame, Cell, dx_from_center, pole_length);
            
            allCorrs{Frame}(counter) = corrMatrix(2,1);
            cellIdList{Frame}(counter) = Cell;
        catch
            %If there is an error with the correlation calculation display what
            %cell what improperly calculated
            disp(['Cell ID ',num2str(Cell),', meshData number ',num2str(counter),...
                ' in Frame ',num2str(Frame), ' correlation was not calculated'])
            allCorrs{Frame}(counter) = NaN;
        end
        counter = counter + 1;
    end
    
end
end