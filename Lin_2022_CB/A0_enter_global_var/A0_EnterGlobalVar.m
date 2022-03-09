
% Generate global variables for the following data analysis pipeline

path = 'O:\Shares\Data_04\Wei-Hsiang Lin\';

%(1) Enter script_folder name
script_folder_name = 'Workflow_PhaseV5\';
script_dir = strcat(path, script_folder_name); 

%(2) Enter data_folder name
data_folder_name = 'WHLin_data_do_not_backup\20220124\';
data_folder =  strcat(path, data_folder_name);

%(3) Enter sub_folder name
input_subfolder = 'waf7_001';
main_folder = strcat(data_folder, input_subfolder);

%(4) Enter maximal digit of image index 
%    For example, for 100 to 999 images, enter '3'
%                 for 1000 to 9999 images, enter '4'
max_digit = 4;

%(5) Enter number of channels
ChNum = 3;

%(6) Enter flag name for each channels
InpChM = cellstr(['c1';'c2';'c3']);
ExpChM = cellstr(['c1';'c2';'c3']);
InpFdM = cellstr(['c1';'c2';'c3']);

%(7) Enter the time interval of each channels
ChannelFlag = [1 6 6];  % image interval for each channel 

