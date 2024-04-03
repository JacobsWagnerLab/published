%% 
% Script to count the number of spots in identified yersinia cells. 
% Made by Joshua McCausland in the CJW lab, 2024.
% Clear the workspace before starting to stay tidy.
clear;
clc;
close all;

% Variables to edit.
% CJW lab scopes have a 65nm pixel size.
px_size = 0.065; % in microns

%load directory
%Input the path to your thunderstorm files in dir(''). 
ts_dir = dir('N:\Collaborations\Vicki Auerbuch Stone\20230308 pcnb null non inducing\replicate 2\ts_results\*.csv');

% Load the meshes from oufti
% this assumes you are in the same directory as the oufti file.
% Put the oufti file name in load('')
load('rep2 mesh new +signal_curated.mat')

% Now count up all the spots per cell. 
regions = cellList.meshData;
spots = [];
for idx = 1:length(ts_dir)
    % Call the current region.
    region = regions{idx};
    % load the current region's thunderstorm results.
    ts_results = readtable([ts_dir(idx).folder '\' ts_dir(idx).name]);
    
    %Within the region, iterate through every identified cell and count
    %spots.
    for qq = 1:size(region, 2)    
        cell = region{qq}.mesh;
        %Create polygon from Oufti mesh.
        pgon = polyshape([cell(:,1);cell(:,3)]*px_size,[cell(:,2);cell(:,4)]*px_size);

        % Thunderstorm's coordinates are in nanometers. Divide that by 1000
        % to match coordinates for regions.
        % Matlab's inpolygon function counts spots in a particular polygon.
        points_in_cell = inpolygon(ts_results.x_nm_/1000,ts_results.y_nm_/1000,pgon.Vertices(:,1),pgon.Vertices(:,2));
        if sum(points_in_cell) > 0      
               spots = [spots; sum(points_in_cell)];
        end
    end
    
end

%%
% Make a quick histogram to view the spots.
histogram(spots, 10)
disp(['The average number of spots per cell in region 1 is ' num2str(mean(spots)) '.'])