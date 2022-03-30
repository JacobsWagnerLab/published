
function [output] = bin_2D_data_V3(var1, var2, data, grid, minimal_sample_size)

bin_num = floor( (grid.max - grid.min) ./ grid.size ) + 2;

count_matrix = zeros(bin_num(1), bin_num(2));
data_matrix = zeros(bin_num(1), bin_num(2));

for j = 1:length(data)
    
    ind_var1 = floor( ( var1(j) - grid.min(1) ) / grid.size(1) );
    ind_var2 = floor( ( var2(j) - grid.min(2) ) / grid.size(2) );
    
    if (ind_var1 < 1)
        bin1 = 1;
    elseif (ind_var1 >= 1) && (ind_var1 < bin_num(1)-1)
        bin1 = ind_var1 + 1;
    elseif (ind_var1 >= bin_num(1)-1)
        bin1 = bin_num(1);
    end

    if (ind_var2 < 1)
        bin2 = 1;
    elseif (ind_var2 >= 1) && (ind_var2 < bin_num(2)-1 )
        bin2 = ind_var2 + 1;
    elseif (ind_var2 >= bin_num(2)-1)
        bin2 = bin_num(2);
    end    
        
    count_matrix(bin1, bin2) = count_matrix(bin1, bin2) + 1;    
    data_matrix(bin1, bin2) =  data_matrix(bin1, bin2) + data(j);        
    
end

data_matrix = data_matrix./ count_matrix;

%%%

filter_mask = (count_matrix > minimal_sample_size);

for r1 = 1:bin_num(1)
    
    for r2 = 1:bin_num(2)
        
        if (filter_mask(r1,r2) == 0)

            count_matrix(r1,r2) = 0;
            data_matrix(r1,r2) = NaN;

        end
        
    end 
    
end

%%% 

bin_tick1 = NaN(bin_num(1),2);
bin_tick1(1,:) = [-Inf grid.min(1)];

for r1 = 1:bin_num(1)-2    
    bin_tick1(1+r1,:) = grid.min(1) + grid.size(1)*[r1-1 r1]; 
end

bin_tick1(end,:) = [grid.max(1) Inf];

%

bin_tick2 = NaN(bin_num(2),2);
bin_tick2(1,:) = [-Inf grid.min(2)];

for r2 = 1:bin_num(2)-2    
    bin_tick2(1+r2,:) = grid.min(2) + grid.size(2)*[r2-1 r2]; 
end

bin_tick2(end,:) = [grid.max(2) Inf];

%%% 

output = {};

output.data = data_matrix;
output.count = count_matrix;
output.bin_tick1 = bin_tick1;
output.bin_tick2 = bin_tick2;


% ======================================================================



% =================================================================

% (1) Count 2D figure

figure('position', [1 1 600 325]);

subplot(121);

ytickslabels = output.bin_tick1(2:end,1);
yticks = (2:size(ytickslabels,1))';
xtickslabels = output.bin_tick2(2:end,1);
xticks = (2:size(xtickslabels,1))';

RYBmap = customcolormap(linspace(0,1,11), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd','#def4f9','#abd9e9','#73add2','#4873b5','#313691'}, 128);
RYBmap(1,:) = 0.6*[1 1 1];   % Making threshold level, below are labebled as gray

imagesc(output.count);
set(gca,'YDir','normal');

set(gca, 'XTick', xticks, 'XTickLabel', xtickslabels);
set(gca, 'YTick', yticks, 'YTickLabel', ytickslabels);
xtickangle(90);

colorbar('northoutside');
colormap(RYBmap);

% =================================================================
% (1) GR 2D figure

subplot(122);

ytickslabels = output.bin_tick1(2:end,1);
yticks = (2:size(ytickslabels,1))';
xtickslabels = output.bin_tick2(2:end,1);
xticks = (2:size(xtickslabels,1))';

imagesc(output.data, 0.001*[4.5 8]);
set(gca,'YDir','normal');

set(gca, 'XTick', xticks, 'XTickLabel', xtickslabels);
set(gca, 'YTick', yticks, 'YTickLabel', ytickslabels);
xtickangle(90);

colorbar('northoutside');
colormap(RYBmap);
