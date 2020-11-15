%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%script: count_locs.m
%@author: Molly Scott
%@date: March 4, 2016
%==========================================================================
%************************Output**********************:
%percents:              Matrix containing the percent of cells with each
%                       number of zones (ranges 0-8)
%************************Input**********************:
%cells:                 The output from peakFinder_cellCheck.m.
%==========================================================================
%This simple script searches through the cell structure, cells, to count
%the number of zones, but excludes the "zones" that occur at poles of
%cells. Outputs the numbers as percents of the total number of cells.
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 


count = 0;
for ii = 1:length(cells)
    if ~isempty(cells{ii})
        if isfield(cells{ii},'locs')
            count = count + 1;
        end
    end
end

all_locs = zeros(count,1);
inner_locs = zeros(count,1);
pole_locs = zeros(count,1);
count = 0;

%determine whether the determined locs are at the poles 
for ii = 1:length(cells)
    if ~isempty(cells{ii})
        if isfield(cells{ii},'locs')
            count = count + 1;
            all_locs(count) = length(cells{ii}.locs);
            for jj = 1:length(cells{ii}.locs)
                rel_loc = cells{ii}.locs(jj)/cells{ii}.length;
                if (rel_loc > 0.1) && (rel_loc < 0.9)
                    inner_locs(count) = inner_locs(count) + 1;
                else
                    pole_locs(count) = pole_locs(count) + 1;
                end
            end
        end
    end
end

%count the number of zones contained within the cell, excluding the poles
for ii = 1:length(inner_locs)
    num_zero = sum(inner_locs(:) == 0);
    num_one = sum(inner_locs(:) == 1);
    num_two = sum(inner_locs(:) == 2);
    num_three = sum(inner_locs(:) == 3);
    num_four = sum(inner_locs(:) == 4);
    num_five = sum(inner_locs(:) == 5);
    num_six = sum(inner_locs(:) == 6);
    num_seven = sum(inner_locs(:) == 7);
    num_eight = sum(inner_locs(:) == 8);
end

total_cells = length(cells);

%get the percent of cells with each number of zones
perc_zero = num_zero/total_cells;
perc_one = num_one/total_cells;
perc_two = num_two/total_cells;
perc_three = num_three/total_cells;
perc_four = num_four/total_cells;
perc_five = num_five/total_cells;
perc_six = num_six/total_cells;
perc_seven = num_seven/total_cells;
perc_eight = num_eight/total_cells;

percents = [perc_zero,perc_one,perc_two,perc_three,perc_four,perc_five,perc_six,perc_seven,perc_eight];






