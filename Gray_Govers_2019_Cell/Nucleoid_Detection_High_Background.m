function Nucleoid_Detection_High_Background(cellList_directories,nucleoid_keyword,cellList_save_name)
%{
-About-
This code allows for thresholding images with high internal fluorescence 
in order to determine nucleoid size. This is especially useful for the
dnaC2 E. coli mutant that has high background of HU fluorescence as the
cells filament. It is done by using just the part of the image that is
within the Oufti cell mesh and then calculating multiple threshold levels
and using the highest one (default is two threshold levels).

-Inputs-
cellList_directories:  these are the directories and names of files that 
contain the cellList to be analyzed. This requires that the MATLAB dir
function is used to create the proper structure of directory file
(containing both the folder and file name). This enables the user to
analyze multiple datasets at once.

nucleoid_keyword:  this is the keyword used to find the raw images that
will be used to detect the nucleoids. An example is nucleoid_keyword = 'c2.tif';

cellList_save_name: the name of the file that you want to save that
contains the detected nucleoids.

-Outputs-
Saves new file in the directory of the cellList that contains the datset
with the objects detected.

-Example-
cellList_directories = dir('**\*cellLists_to_be_analyzed.mat');
nucleoid_keyword = 'c2.tif';
cellList_save_name = 'cellList_plus_nucleoid';
nucleoid_detection_new(cellList_directories,nucleoid_keyword,cellList_save_name)

-Keywords-
Nucleoid
Detection

-Dependencies-
Requres matlab's image processing toolbox

-Author-
William Gray October 10, 2018
%}

%Declare parameters for the segmentation
plot_figures = false;
frames2plot = 1:2;
close_nucleoids = true;
thresh_levels = 2;
image_bits = 16;
pixel_size = 0.064;
disk_radius = 5;

%Iterate through all datsets supplied
for datasets = 1:length(cellList_directories)
    cd([cellList_directories(datasets).folder])
    
    tic
    
    % load cellList
    cellList_data = load([cellList_directories(datasets).folder,'\',cellList_directories(datasets).name],'cellList');
    cellList = cellList_data.cellList;

    % load images
    images = dir(['**\*',nucleoid_keyword]);

    meshData = cellList.meshData;
    nucleoid_size = cell(1,length(meshData));
    nucleoid_length = cell(1,length(meshData));
    nucleoid_width = cell(1,length(meshData));
    nucleoid_diameter = cell(1,length(meshData));
    parfor ii = 1:length(meshData)

        %Load image to be thresholded
        im2thresh = imread([images(ii).folder,'\',images(ii).name]);

        %Create variables for calculations
        segmented_image_sum = ones(size(im2thresh));
        all_closed_nucleoid_images = zeros(size(im2thresh));

        %Go through each cell and create a mask and perform segmentation on
        %the raw image after selection of the region within that mask
        for jj = 1:length(meshData{ii})

            %Create a cellMesh array that can be input into the poly2mask
            %command
            if length(meshData{ii}{jj}.mesh) > 15
                cellMesh = double( cat(1, meshData{ii}{jj}.mesh(:,1:2),...
                    flipud(meshData{ii}{jj}.mesh(:,3:4))) );

                %Convert the outline of the cell into a mask
                current_cell_mask = poly2mask(cellMesh(:,1),cellMesh(:,2),2048,2048);

                %Select just the region of the maks in the raw image
                cell2thresh = immultiply(im2thresh,current_cell_mask);

                %Perform thresholding on the selected region that contains the cell
                %of interest
                thresh = multithresh(cell2thresh,thresh_levels);

                %Convert the threshold into a segmented image (for visualization
                %purposes)
                segmented_image = imquantize(cell2thresh,thresh);
                segmented_image_sum = segmented_image_sum + (segmented_image - 1);

                %This performs the actual single level thresholding on the nucleoid
                %based on the threshold level that is found using the multithresh
                %function.
                nucleoid_image = imbinarize(cell2thresh,max(double(thresh))./(2^image_bits-1));

                %This portion of the code will close nucleoids if they are nearby
                %(works well if there are heterogeneous nucleoids that are
                %thresholded into several nucleoids). You will want to try altering
                %the structuring element disk_radius to see what performs best.
                if close_nucleoids
                    se = strel('disk',disk_radius);
                    closed_nucleoid_image = imclose(nucleoid_image,se);
                else
                    closed_nucleoid_image = nucleoid_image;
                end

                %Add up all of the closed nucleoid images to determine how the
                %performance is for your dataset (for visualization purposes only)
                all_closed_nucleoid_images = all_closed_nucleoid_images + closed_nucleoid_image;

                area_holder = regionprops(logical(closed_nucleoid_image),'Area');
                length_holder = regionprops(logical(closed_nucleoid_image),'MajorAxisLength');
                width_holder = regionprops(logical(closed_nucleoid_image),'MinorAxisLength');
                diameter_holder = regionprops(logical(closed_nucleoid_image),'EquivDiameter');
                
                %Assign variables to be added to the cellList meshData
                nucleoid_size{ii}{jj} = [area_holder(:).Area]*(pixel_size^2);
                nucleoid_length{ii}{jj} = [length_holder(:).MajorAxisLength]*(pixel_size);
                nucleoid_width{ii}{jj} = [width_holder(:).MinorAxisLength]*(pixel_size);
                nucleoid_diameter{ii}{jj} = [diameter_holder(:).EquivDiameter]*(pixel_size);
                
                
                meshData{ii}{jj}.objectArea = nucleoid_size{ii}{jj};
                meshData{ii}{jj}.objectLength = nucleoid_length{ii}{jj};
                meshData{ii}{jj}.objectWidth = nucleoid_width{ii}{jj};
                meshData{ii}{jj}.objectDiameter = nucleoid_diameter{ii}{jj};
            else
                meshData{ii}{jj}.objectArea = 0;
            end
        end
        
        %In order to perform tests on how your code is working, you can
        %save and then view figures that display the thresholding and
        %closing of black-white images
        if plot_figures
            for kk = 1:length(frames2plot)
                if abs(frames2plot(kk)-ii) < 1e-5
                    figure;imshow(im2thresh,[min(im2thresh(:)),max(im2thresh(:))],'InitialMagnification','fit')
                    print(['raw_image',sprintf('%03d',ii)],'-dpng')
                    print(['raw_image',sprintf('%03d',ii)],'-depsc')
                    savefig(['raw_image',sprintf('%03d',ii)])
                    close
                    
                    RGB = label2rgb(segmented_image_sum);
                    figure;imshow(RGB,'InitialMagnification','fit')
                    print(['thresholded_image',sprintf('%03d',ii)],'-dpng')
                    print(['thresholded_image',sprintf('%03d',ii)],'-depsc')
                    savefig(['thresholded_image',sprintf('%03d',ii)])
                    close
                    
                    figure;imshow(all_closed_nucleoid_images,'InitialMagnification','fit')
                    print(['image_after_closing',sprintf('%03d',ii)],'-dpng')
                    print(['image_after_closing',sprintf('%03d',ii)],'-depsc')
                    savefig(['image_after_closing',sprintf('%03d',ii)])
                    close
                end
            end
        end
    end
    
    cellList.meshData = meshData;
    save(cellList_save_name,'cellList')
    elapsed_time = toc;
    disp(['cellList ',num2str(datasets),' completed in ',num2str(elapsed_time),' seconds']) 
end