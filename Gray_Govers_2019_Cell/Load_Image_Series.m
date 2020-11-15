function combined_image = Load_Image_Series(path_to_folder)
%{
-About-
Function that loads in a series of individual .tif images and
combines them into a stack

-Inputs-            
path_to_folder:   full path to folder containing separate images to be combined. 
 
-varargin-
NA

-Outputs-
combined_image:   stack of all .tif images in the specifed folder

-Example-
NA

-Supplementary-
NA

-Keywords-
Image stack, combine multiple images

-Dependencies-

-References-
NA

-Author-
Manuel Campos
edited by Sander Govers, 22 October 2018
%}



%Only perform is specified path is that to an actual folder
if isfolder(path_to_folder)
    %Search for all .tif files and store them in a directory
    files = dir([path_to_folder '/*.tif*']);
    %Identify first file
    firstFile=imread([path_to_folder '/' files(1).name]);
    %Define return stack of images, based on size of first image
    combined_image=uint16(zeros(size(firstFile,1),size(firstFile,2),length(files)));
    %Initiate progress bar
    w = waitbar(0, 'Loading image files, please wait...');
    %Loop through all .tif files and insert them into the stack of images
    %the function will return
    for jj=1:length(files)
        %Assign image to correct position in stack 
        combined_image(:,:,jj) = imread([path_to_folder '/' files(jj).name]);
        %Increase progress bar value
        waitbar(jj/length(files));
    end
end
close(w)
end