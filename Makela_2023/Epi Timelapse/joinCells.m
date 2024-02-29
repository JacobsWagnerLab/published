function joinCells
% remove lines that have empty no_chrom and joins the cells into a single
% file

% get cells file
[filename,pathname] = uigetfile('*.mat', 'Select','Select','MultiSelect', 'on');
if ~iscell(filename)
    name = filename;
    filename = cell(1);
    filename{1} = name;
end

disp(num2str(pathname))

for hh = 1:length(filename)
    load([pathname filename{hh}])
    
    if isfield(cells,'no_chrom')
        % remove indeces without no_chrom
        indeces = [];
        for ii = 1:length(cells)
            if isempty(cells(ii).no_chrom) || cells(ii).no_chrom == 0
                indeces = [indeces ii];
            end
        end
        cells(indeces) = [];

        % concatenate structs
        if hh == 1
            new_cells = cells;
        else
            new_cells = [new_cells, cells];
        end
    else
        % concatenate structs
        if hh == 1
            new_cells = cells;
        else
            new_cells = [new_cells, cells];
        end
    end
end
cells = new_cells;

% remove not useful fields
if isfield(cells,'fluo1_corr')
    cells = rmfield(cells,'fluo1_corr');
end
    
% save combined cells with a new name '_allcells_cropped'
splitted_name = split(filename{hh},char(strcat('_XY', num2str(hh))));
new_name = char(strcat(splitted_name{1}, '_cells_all.mat'));
save([pathname new_name],'cells','-v7.3')

end

