%%
% Author: Silvia, February 05 2024
% Importing segmentation results after SGO and constriction analysis

clc
clear all
close all
%For each FOV, load the labelled image (mask) from SuperSegger-Omnipose output
L=imread('F:\Shares\Data_03\SilviaCanasDuarte\DATA\InVivo_Glycogen\20220824_YT\Merged3\xy01\masks\220824_sb17006_m9glucaat_od283t00001xy001c1_cp_masks.png');
L0=L; %Creates a copy of the labelled image

%%
% Find the touching labels to identify cells that are constricting but were splitted (segmented) as two cells too early
touching_labels = [];
for i = 1:max(L(:))
    % Get the bounding box of the current label
    bbox = regionprops(L == i, 'BoundingBox');
    bbox.BoundingBox = bbox.BoundingBox + [-1 -1 2 2];
    
    % Check if the bounding box overlaps with any other label
    for j = i+1:max(L(:))
        bbox2 = regionprops(L == j, 'BoundingBox');
        if rectint(bbox.BoundingBox, bbox2.BoundingBox) > 0 
            touching_labels = [touching_labels; i j];
        end
    end
end
%%
% Iterate over the list of touching labels and ask the user if they should be merged
for i = 1:size(touching_labels, 1)
    % Get the labels to merge
    label1 = touching_labels(i, 1);
    label2 = touching_labels(i, 2);

    % Show the pair of labels

    figure(1);
    imshow(ismember(L, [label1 label2]));
    title(sprintf('Labels %d and %d', label1, label2));
    zoom on;
    pause;
    zoom off;

    
    % Ask the user if they want to merge the labels
    answer = questdlg(sprintf('Do you want to merge labels %d and %d?', label1, label2), 'Merge Labels', 'Yes', 'No', 'Yes');
    
    if strcmp(answer, 'Yes')
        % Merge the labels
        L(L == label2) = label1;
    end
end

% Display the merged image
RGB = label2rgb(L, 'jet', 'w', 'shuffle');
imshow(RGB);
imwrite(L, 'labeled_image.png', 'BitDepth', 8); %Overwrites the original labelled image with the updated labels (masks)

%%
