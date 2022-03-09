
function [ret1] = export_corrected_separator_data(kmat, peaks, separatorC, param, export_dir, plot_flag)

% Export peak figures (kymograph and separator plot) and .xls file.

% (A) Summary figures. Each figure has 1000 frames

fn_range = param.Nrange;

if ( plot_flag(1) == 1 ) 
    
fnA_per_figure = 1000;
tick_sizeA = 50;

numA = ceil ( (fn_range(2) - fn_range(1) + 1)/fnA_per_figure ); 
fnA = NaN(numA, 2);

for r = 1:numA
    
    fnA(r,1) = fn_range(1) + (r-1)* fnA_per_figure;
    fnA(r,2) = fnA(r,1) + fnA_per_figure + 5; %create some overlap between figures
    
    fnA(r,1) = max(1, fnA(r,1));
    fnA(r,2) = min(fn_range(2), fnA(r,2));
        
end

cd(export_dir);

for r = 1:numA
    
    name_temp = strcat('peak_data_corrected_A_', num2str(fnA(r,1)), '-', num2str(fnA(r,2)), '.jpg');
    
    fn_range_temp = fnA(r,:);
     
    plot_peak_kmat2(kmat, peaks, separatorC, fn_range_temp, tick_sizeA, name_temp);    
    
    close all;    
    
end

end


% (B) "Zoom in" figures. Each figure has 100 frames

if (plot_flag(2) == 1)

fnB_per_figure = 100;
tick_sizeB = 5;

numB = ceil ( (fn_range(2) - fn_range(1) + 1)/fnB_per_figure ); 
fnB = NaN(numB, 2);

for r = 1:numB
    
    fnB(r,1) = fn_range(1) + (r-1)* fnB_per_figure;
    fnB(r,2) = fnB(r,1) + fnB_per_figure + 5; %create some overlap between figures
    
    fnB(r,1) = max(1, fnB(r,1));
    fnB(r,2) = min(fn_range(2), fnB(r,2));
        
end

cd(export_dir);

for r = 1:numB
    
    name_temp = strcat('peak_data_corrected_B_', num2str(fnB(r,1)), '-', num2str(fnB(r,2)), '.jpg');
    
    fn_range_temp = fnB(r,:);
     
    plot_peak_kmat2(kmat, peaks, separatorC, fn_range_temp, tick_sizeB, name_temp);    
    
    close all;    
    
end


end

ret1 = 1;

end

%=========================================================================

function [ret1] = plot_peak_kmat2(kmat, peaks, separatorC, fn_range, tick_size, name)

hfig = figure('position', [1 1 1500 400]);

% (a) Plot kmat in the fn_range
imagesc( kmat( :,fn_range(1):fn_range(2) ) ); hold on;
colormap(gray);

% (b) Plot peaks

for fn = fn_range(1):(fn_range(2)-1)

    peak_num = size(peaks{fn}.dataC,1);
    
    if (peak_num > 0)    

        for j = 1:peak_num
            
            link_test = (peaks{fn}.dataC(j,3) > 0) + (peaks{fn}.dataC(j,4) > 0) ;            
            
            if (link_test == 2)
                
                %plot(fn, peaks{fn}.data(j,1), 'k-'); hold on;    
                
            elseif (link_test < 2) 
                
                ind_col = fn - fn_range(1) + 1;
                
                plot(ind_col, peaks{fn}.dataC(j,1), 'b.'); hold on;  
                                            
            end
            
        end
        
    end
    
end


for m = 1: length(separatorC)
    
    data = separatorC{m}.data;
    
    % Find if the track data overlapped with the fn_range
    % if yes, return the overlapped data region.
    
    data_temp = data;
    data_temp( (data_temp(:,1) < fn_range(1)) ,:) = [];
    data_temp( (data_temp(:,1) > fn_range(2)) ,:) = [];
    
    if ( size(data_temp, 1) > 0 ) && ( sum(separatorC{m}.colorB) < 3 )
        
        plot(data_temp(:,1)-fn_range(1)+1, data_temp(:,3), '-',...
            'color', separatorC{m}.colorA, ...
            'LineWidth', 2); hold on;        
            
    elseif ( size(data_temp, 1) > 0 ) && ( sum(separatorC{m}.colorB) == 3 )
    
        plot(data_temp(:,1)-fn_range(1)+1, data_temp(:,3), 'o-', ...
            'color', [1 0 0], ...
            'LineWidth', 1 ); hold on;
        
    end

end


xdata = fn_range(1):tick_size:fn_range(2);

xticks(xdata-xdata(1)+1);
xticklabels(xdata);

saveas(hfig, name);

ret1 = 1;


end

%=========================================================================



