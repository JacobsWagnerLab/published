
% Get peak data for each column of kymograph matrix (kmat)

function [peaks] = kmat_to_peaks(kmat, peaks, paramA)

% Define a column subrange of kmat for finding peaks
fn_range = paramA.Nsubrange;

% If "append_option" == 0, reset all peak data.
% Otherwise, partially overwrite the peak data.

if ( paramA.append_option == 0 )   % Reset "peaks" as an empty dataset
    peaks = {};
end

% Find peaks for each column, an write the data into "peaks{}".

for fn = fn_range(1):fn_range(2)
    
    data_temp = get_local_max(kmat(:,fn), paramA);  
    
    peaks{fn}.data = [data_temp 0*data_temp 0*data_temp];  
    % col 1: y-position
    % col 2: peak intensity
    % col 3,4,5,6 reserved for linking record
    
end

end

% ================================================================

function [peak_vec] = get_local_max(col,  param)

% Example for parameters:
% width = 10
% min_level = -6;
% offset = 0.5;

w = param.w;                    % Radius (in pixel) for peak finding 
min_level = param.min_level;    % Minimal level for peak finding    
offset = param.offset;          % Peak threshold parameter. Peak must larger than (mean + offset*SD )


% (1) Find local peak position

peak_vec1 = [];   % (1) ypos, (2) peak intensity

for j = 1:length(col)
    
    rgL = max(1, j-w);
    rgR = min(length(col), j+w);
    
    temp = col(rgL:rgR);    
    
    local_max = max(temp);
    local_max_pos = find(temp == local_max);
    local_max_count = length(local_max_pos);    
    local_threshold = mean(temp) + offset *std(temp);
    
    if ( local_max_count  == 1 )  % finding a unique local peak
    
        if ( col(j) == local_max ) 
        
            if ( ( col(j) > min_level ) && ( col(j) > local_threshold ) )
        
                peak_vec1 = [peak_vec1; j col(j)];
        
            end
            
        end
        
    elseif ( local_max_count == 2 )  % two points with same value. Check if they are adjacent    
        
        
        if  ( abs( local_max_pos(2)-local_max_pos(1) ) == 1 ) 
            
            if ( col(j) == local_max ) && ( j == local_max_pos(1)+rgL-1 )
            
                if ( ( col(j) > min_level ) && ( col(j) > local_threshold ) )
            
                    peak_vec1 = [peak_vec1; j col(j)];
                    
                end
            end
                
        end      
        
    end
    
end

% (2) Find plateau region 

plateau_value = -0.01;
plateau_vec = ( col > plateau_value ); 
plateau_seg = bwlabel(plateau_vec);
num = max(plateau_seg);

peak_vec2 = NaN(num, 2);

% Finding the initial and ending point of each plateau segment

if ( num > 0 )
    
    plateau_range = NaN(num, 3);  % initial, end, mid

    for m = 1:num
        
        temp = find(plateau_seg == m);
        
        temp_ind = (temp == 1);
        plateau_range(m,:) = [min(temp) max(temp) floor( (min(temp) + max(temp)) / 2 )];
                
    end
    
    for m = 1:num
        
        peak_vec2(m,:) = [plateau_range(m,3) col(plateau_range(m,3))] ;        
        
    end
    
end
    
% Remove the overlapped peaks between peak1 and peak2 

for m1 = 1:size(peak_vec1)
    
    for m2 = 1:size(peak_vec2)
        
        if ( abs(peak_vec2(m2) - peak_vec1(m1) ) <= 2 )
        
            peak_vec2(m2,:) = [];
            
        end
        
    end

end

peak_vec = [peak_vec1; peak_vec2];

peak_vec = sortrows(peak_vec, 1);

% Remove the peak at the exit border (not real peak but at boundary)
% 
% for r = 1:size(peak_vec)
%     
%     if ( peak_vec(r,1) == size(col,1) ) 
%     
%         peak_vec(r,:) = [];
%         
%     end
%     
% end


end

% ====================================================================== 

