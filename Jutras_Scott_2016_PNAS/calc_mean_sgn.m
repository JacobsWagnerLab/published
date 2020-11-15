%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function mean_sgn_perCell = calc_mean_sgn(cellList,signal)
%@author: Molly Scott
%@date: March 4, 2016
%==========================================================================
%************************Output**********************:
%mean_sgn_perCell:      matrix containing the average fluorescence per
%                       cell, normalized by the length of each cell
%************************Input**********************:
%cellList:              output from Oufti
%signal:                'signal 1' or 'signal2'
%==========================================================================
%This function allows you to calculate the average fluorescence in each
%cell, normalized by the length of that cell.
%-------------------------------------------------------------------------- 
%--------------------------------------------------------------------------

function mean_sgn_perCell = calc_mean_sgn(cellList,signal)

count = 0;
for frames = 1:length(cellList.meshData)
    for cells = 1:length(cellList.meshData{frames})
        if ~isempty(cellList.meshData{frames}{cells}) && isfield(cellList.meshData{frames}{cells},(signal)) && isfield(cellList.meshData{frames}{cells},'steparea')
            count = count + 1;
        end
    end
end

count = 0;
normalized_sum_sgn = zeros(count,1);
for frames = 1:length(cellList.meshData)
    for cells = 1:length(cellList.meshData{frames})
        if ~isempty(cellList.meshData{frames}{cells}) && isfield(cellList.meshData{frames}{cells},(signal)) && isfield(cellList.meshData{frames}{cells},'steparea') 
            count = count + 1;
            sgn = cellList.meshData{frames}{cells}.(signal)./cellList.meshData{frames}{cells}.steparea;
            normalized_sum_sgn(count,1) = sum(sgn)/cellList.meshData{frames}{cells}.length;
        end
    end
end

mean_sgn_perCell = mean(normalized_sum_sgn);
sem_tot = std(normalized_sum_sgn)/sqrt(length(normalized_sum_sgn));

end


            