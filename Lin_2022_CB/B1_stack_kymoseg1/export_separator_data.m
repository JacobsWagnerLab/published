
function [peak_linking_table] = export_separator_data(kmat, peaks, paramA, separator, export_dir, plot_flag)

% Export peak figures (kymograph and separator plot) and .xls file.
% The .xls file is used for manaul curation for wrong connection of separators

fn_range = paramA.Nrange;

% (1) Create the peak-linking table

mflag = 0;

for fn = fn_range(1):fn_range(2)

    if ( size(peaks{fn}.data, 1) >= mflag )
    
        mflag = size(peaks{fn}.data,1);
    
    end
    
end
    
% Peak-linking table: each row represent one frame. 
% Column 1:mflag representing linking data (the linked-peak index in the next frame)
% Column [mflag+1] is reserved for making note during manual change (optional)

peak_linking_table = [(-1)*ones(length(peaks), mflag) zeros(length(peaks),1)];

for fn = fn_range(1):fn_range(2)
    
    L = size(peaks{fn}.data,1);
    
    if (L > 0)
        peak_linking_table(fn, 1:L) = peaks{fn}.data(1:L,4);
    end
    
end

% (2) Write the peak-linking table into .xls data. 

% Typical way to use is to modify the "peaks_linking_corrected.xls" 
% and keep the original "peak_linking_original.xls" as reference. 

cd(export_dir);
writematrix(peak_linking_table, 'peak_linking_original.xls', 'WriteMode', 'overwritesheet');
writematrix(peak_linking_table, 'peak_linking_corrected.xls', 'WriteMode', 'overwritesheet');

% -----------------------------------------------------------------------

% (2) Plot the linking results

% (A) Summary figures. Each figure has 1000 frames

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
    
    name_temp = strcat('peak_data_A_', num2str(fnA(r,1)), '-', num2str(fnA(r,2)), '.jpg');
    
    fn_range_temp = fnA(r,:);
     
    plot_peak_kmat(kmat, peaks, separator, fn_range_temp, tick_sizeA, name_temp);    
    
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
    fnB(r,2) = fnB(r,1) + fnB_per_figure + 5;   %create some overlap between figures
    
    fnB(r,1) = max(1, fnB(r,1));
    fnB(r,2) = min(fn_range(2), fnB(r,2));
        
end

cd(export_dir);

for r = 1:numB
    
    name_temp = strcat('peak_data_B_', num2str(fnB(r,1)), '-', num2str(fnB(r,2)), '.jpg');
    
    fn_range_temp = fnB(r,:);
     
    plot_peak_kmat(kmat, peaks, separator, fn_range_temp, tick_sizeB, name_temp);    
    
    close all;    
    
end


end

ret1 = 1;


end


%=========================================================================

function [ret1] = plot_peak_kmat(kmat, peaks, separator, fn_range, tick_size, name)

hfig = figure('position', [1 1 1500 400]);

% (a) Plot kmat in the fn_range
imagesc( kmat( :,fn_range(1):fn_range(2) ) ); hold on;
colormap(gray);

% (b) Plot peaks

for fn = fn_range(1):(fn_range(2)-1)

    peak_num = size(peaks{fn}.data,1);
    
    if (peak_num > 0)    

        for j = 1:peak_num
            
            link_test = (peaks{fn}.data(j,3) > 0) + (peaks{fn}.data(j,4) > 0) ;
            link_test2 =  (peaks{fn}.data(j,5) > 0) + (peaks{fn}.data(j,6) > 0) ;
            
            if (link_test == 2)
                
                %plot(fn, peaks{fn}.data(j,1), 'k-'); hold on;    
                
            elseif (link_test < 2) && (link_test2 == 2 )
                
                ind_col = fn - fn_range(1) + 1;
                
                plot(ind_col, peaks{fn}.data(j,1), 'bo'); hold on;  
                
            elseif (link_test < 2) && (link_test2 < 2 )
                
                ind_col = fn - fn_range(1) + 1;
                
                plot(ind_col, peaks{fn}.data(j,1), 'b.'); hold on;  
                
            end
            
            
            
        end
        
    end
    
end


for m = 1: length(separator)
    
    data = separator{m}.data;
    
    % Find if the track data overlapped with the fn_range
    % if yes, return the overlapped data region.
    
    data_temp = data;
    data_temp( (data_temp(:,1) < fn_range(1)) ,:) = [];
    data_temp( (data_temp(:,1) > fn_range(2)) ,:) = [];
    
    if ( size(data_temp, 1) > 0 ) && ( sum(separator{m}.colorB) < 3 )
        
        plot(data_temp(:,1)-fn_range(1)+1, data_temp(:,3), '-',...
            'color', separator{m}.colorA, ...
            'LineWidth', 2); hold on;        
            
    elseif ( size(data_temp, 1) > 0 ) && ( sum(separator{m}.colorB) == 3 )
    
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



