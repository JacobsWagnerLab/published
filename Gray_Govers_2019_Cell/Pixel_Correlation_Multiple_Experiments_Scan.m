function Pixel_Correlation_Multiple_Experiments_Scan(directory_files, signal_names, dx_from_center, pole_length, save_name)
%{
-About-
This function uses raw images and an Oufti cellList to calculate the pixel
by pixel correlation between two fluorescent signals. The function enables 
the user to run pixel correlation on many different cellList files with 
minimal user input. The function takes advantage of MATLAB's parfor loop to
calculate correlations in parallel and perform the analysis more quickly. 
Any two signals can be analyzed but your Oufti cellList must have at least 
two signals for the analysis to work.


-Inputs-
directory_files = dir('**\Directory where your files are located')
this is the directory where all of your files are located. You can clean
this to contain any/all files you are interested.

signal_names = {'GFP','DAPI'};
These area the names of the folders that contain the two signals that you 
want to analyze. For the code to run properly, you need to have separate
image files for each image frame that are in two different folders with
corresponding numbers:
(i.e. GFP\image001c1.tif,GFP\image002c1.tif;DAPI\image001c2.tif,DAPI\image002c2.tif)

dx_from_center = 1 or numerical array;
This determines how many pixels away from the centerline you want to
calculate the correlation. Must add up to an integer. A value of one
means that the correlation area will be two pixels wide. In order to scan
multiple values, input an array (i.e. [1:0.5:4])

pole_length = 8 or numerical array;
This is the number of pixels away from the pole that you want to
calculate the correlation for. This is to avoid the artificial positive
correlation at the pole that results from the decrease in volume there.
Same as above, input an array if you want to scan parameters.

save_name = 'pixel_correlation_file';
This is the name of the file that will be saved in each folder that
contains your analysis files/images.

-Outputs-
Saved file with the pixel correlations calculated, corresponds to the
length and organization structure of the cellList. Also saves the
iteration_key which species the parameter combinations that were used.

-Example-
directory_files = dir('**\*cellLists_to_be_analyzed.mat');
signal_names = {'GFP','DAPI'};
dx_from_center = 1:0.5:4;
pole_length = 1:20;
save_name = 'pixel_correlation_file';
Pixel_Correlation_Multiple_Experiments_Scan(directory_files, signal_names,...
 dx_from_center, pole_length, save_name)

-Keywords-
Pixel correlation
Ribosome segregation

-Dependencies-
Pixel_Correlation_Parallel.m
getTIFsFromPaths.m
Cell_Pixel_Correlation.m
Extract_Cell_Pixels.m
Cell_Projection.m
Taylor_Smooth.m

-Author-
William Gray October 10, 2018
%}

%This will go through all of your different signals and find the folder
%that contains that signal (the folder must contain some portion of the
%same name as the signal_names variable).

for all_files = 1:length(directory_files)
    all_correlations = cell(1,length(dx_from_center)*length(pole_length));
    
    %Create the cell arrays to speed up processing
    signal_folder_name = cell(1,length(signal_names));
    
    %This is to hold the location of the images that you will be analyzing 
    %for the correlations
    imPaths = cell(1,length(signal_names));
    
    %Iterate through all of the signals and all of the images
    for ii = 1:length(signal_names)
        signal_folder_name{ii} = dir([directory_files(all_files).folder,'\*',signal_names{ii},'*']);
        signal_folder_holder = [];
        %Check to make sure that the dir command is finding only directories
        for jj = 1:length(signal_folder_name{ii})
            if signal_folder_name{ii}(jj).isdir
                signal_folder_holder = signal_folder_name{ii}(jj);
            end
        end
        %Reassign the signal folder variable
        signal_folder_name{ii} = signal_folder_holder;
        
        %This finds the actual paths of the tif files for each signal
        imPaths{ii} = [directory_files(all_files).folder,'\',signal_folder_name{ii}.name];
        
    end
    
    %Specify the cellList_file for the current analysis
    cellList_file = [directory_files(all_files).folder,'\',directory_files(all_files).name];
    
    %Perform the correlation analysis with the specified parameters
    iteration_key = zeros(length(dx_from_center)*length(pole_length),2);
    iteration_counter = 1;
    for center_iterations = 1:length(dx_from_center)
        tic
        for pole_iterations = 1:length(pole_length)
            all_correlations{iteration_counter} = Pixel_Correlation_Parallel(imPaths, cellList_file, dx_from_center(center_iterations), pole_length(pole_iterations));
            iteration_key(iteration_counter,:) = [dx_from_center(center_iterations),pole_length(pole_iterations)];
            iteration_counter = iteration_counter + 1;
        end
        toc
    end
    
    %Create a variable for the individual experiment and save it in the
    %same directory as its corresponding cellList
    all_correlations_ind_expt = all_correlations;
    save([directory_files(all_files).folder,'\',save_name],'all_correlations_ind_expt','iteration_key')
    disp(['Dataset ',num2str(all_files),' analyzed'])
end
end