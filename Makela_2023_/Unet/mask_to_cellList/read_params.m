%% read_params
%
% -Purpose-
%   Parse in the .set Oufti parameter file to a cell array of char arrays.
%   Output from this function was then compiled in generate_params() to
%   create the parameter struct 'p' that is used by Oufti.
%
% -Input-
%   - filename (char array): path to the .set Oufti parameter file
%
% -Output-
%   - paramString (cell array): each cell is a single line read from the
%                Oufti .set parameter file
%
% -Author-
%   Yingjie Xiang, 2019-01-17
%
% -Patch Notes-
%   2019-01-17: created the function

function paramString = read_params(filename)

% Open the Oufti .set parameter file
fileID = fopen(filename);

% Create a buffer to store the lines in the file
buffer = cell(100,1);

% Read the first line
line = fgetl(fileID);

% Check on the first line
if ~ischar(line)
    error('Read parameter set file.');
end

% Store the first line and initiate the count
buffer{1} = line;
count = 2;

% Go through each line and store them into the buffer
while ischar(line)
    line = fgetl(fileID);
    buffer{count} = line;
    count = count + 1;
end

% Always remember to close the file
fclose(fileID);

% The end of file is marked by the number -1, find it here
end_of_file = find(cellfun(@isnumeric,buffer),1,'first');

% Output the paramString used by oufti
paramString = buffer(1:end_of_file-1);

end