
% Get mother cell lineage (MCL)

function [lineage_data, MCstat, mother_cell] = get_MCL(input)
%input = combined_dataset.cc_ensemble{4};

% Assembling lineage data from cc_ensemble

%%% (1) Brief statistics of lineage

[lineage_data] = get_lineage_data(input); 

% Each row is one cell cycle track (may be incomplete)
% [1] initial frame
% [2] ending frame
% [3] parent
% [4] offspring #1
% [5] offspring #2
% [6] checking if this cell cycle has proper division (start and end)


%%% (2) Constructing the mother-cell lienage

[mother_cell] = get_mother_cell_statistics(input, lineage_data);
[mother_cell] = remove_NaN_in_MC(mother_cell);

%%% (3) Perform basic statistics for cell cycle on MC

max_frame = size(mother_cell.vec, 1);

if ( max_frame > 10) 
    [MCstat] = get_MC_stat(mother_cell);
else
    MCstat = {};
end

%{
figure('position', [1 1 800 200]);

subplot(211);
plot(1:max_frame, 10*mother_cell.flag);  hold on;
plot(Rsc(:,1), Rsc(:,2), '-');
ylim([2 8]);

subplot(212);
semilogy(1:max_frame, 1000*(mother_cell.flag+0.1));  hold on;
semilogy(cellsize, 'k-');
ylim([100 1000]);
%}

end

% ================================================================== 

function [mother_cell2] = remove_NaN_in_MC(mother_cell)

first_NaN = find(isnan(mother_cell.vec(:,1)), 1, 'first');

mother_cell2 = mother_cell;

if ( ~isempty(first_NaN ) )  % truncate the lineage data
    
    mother_cell2.vec = mother_cell.vec(1:first_NaN-1,:); 
    mother_cell2.size = mother_cell.size(1:first_NaN-1,:); 
    mother_cell2.signal = mother_cell.signal(1:first_NaN-1,:); 
    mother_cell2.flag = mother_cell.flag(1:first_NaN-1,:); 

end

end

% ================================================================== 

function [mother_cell] = get_mother_cell_statistics(input, lineage_stat)

cc_num = size(lineage_stat,1);
max_frame = max(lineage_stat(:,2));

mother_cell = {};
mother_cell.vec = NaN(max_frame, 2);  % [cc, ypos]

for cc = 1:cc_num
    
    fn = input.data{cc}.data(:,1);    % time frame   
    ypos = input.data{cc}.data(:,8);  % y_position
    
    for r = 1:size(fn,1)        
        
        cc_current = mother_cell.vec(fn(r), 1);
        y_current = mother_cell.vec(fn(r), 2);
        
        if ( isnan(cc_current) )  % not assign yet
            
             mother_cell.vec(fn(r), 1) = cc;
             mother_cell.vec(fn(r), 2) = ypos(r);
        
        else  % not NaN
                
            if ( ypos(r) < y_current)  % Having smaller y
        
             	mother_cell.vec(fn(r), 1) = cc;
                mother_cell.vec(fn(r), 2) = ypos(r);                
                
            end
            
        end
        
    end
    
end

% For MC lineage, go through each frame to append (1) cell size and 
% (2) intensity info. 

mother_cell.size = NaN(max_frame, 1);

for fn = 1:max_frame
    
    cc = mother_cell.vec(fn,1);
    
    if (~isnan(cc))
    
        data_cc = input.data{cc}.data;
        
        ind = find(data_cc(:,1) == fn, 1, 'first');
        mother_cell.size(fn) = data_cc(ind,3);  % cell size
    
    end
    
end


mother_cell.signal = NaN(max_frame, 2);

for fn = 1:max_frame
    
    cc = mother_cell.vec(fn,1);
    
    if (~isnan(cc))
    
    signal = input.data{cc}.signal;
    
    if ( ~isempty(signal) )
        
        ind = find(signal(:,1) == fn, 1, 'first');
    
        if ( ~isempty(ind) )
        
           mother_cell.signal(fn,:) = signal(ind, 3:4);  % 405, 488 intensity
    
        end 
    
    end
    
    end
    
    
end

mother_cell.flag = NaN(max_frame, 1);

mother_cell.flag(1) = 1;

for fn = 1:(max_frame-1)
    
    if ( mother_cell.vec(fn+1,1) ~= mother_cell.vec(fn,1)  )
        
        mother_cell.flag(fn+1) = (-1)*mother_cell.flag(fn);        
        
    elseif ( mother_cell.vec(fn+1,1) == mother_cell.vec(fn,1)  )
        
        mother_cell.flag(fn+1) = mother_cell.flag(fn);    
        
    end
        
end

mother_cell.flag(mother_cell.flag == 1) = 0;
mother_cell.flag(mother_cell.flag == -1) = 1;

end

% ================================================================== 

function [lineage_data] = get_lineage_data(input)

cc_num0 = max(size(input.track_note));

lineage_data = NaN(cc_num0,7);
% [1] initial frame
% [2] ending frame
% [3] parent
% [4] offspring #1
% [5] offspring #2
% [6] checking if this cell cycle has proper division (start and end)

for cc = 1:cc_num0
    
   if ( ~isempty(input.data{cc}.data)  ) 
       lineage_data(cc,1) = input.data{cc}.data(1,1);
       lineage_data(cc,2) = input.data{cc}.data(end,1);   
   end
   
   lineage_data(cc,3) = input.track_note(cc,1); 
   
end


% Remove empty tracks 
lineage_data( isnan(lineage_data(:,1)), : ) = [];



% ==================  Finding offspring cells ======================= 

cc_num = size(lineage_data, 1);

for cc = 1:cc_num
    
    parent_id = lineage_data(cc,3);
    
    if (parent_id > 0)  % This cell is a parent. Find its offsprings
        
    if ( isnan(lineage_data(parent_id,4))  )
    
        lineage_data(parent_id, 4) = cc;

    elseif ( lineage_data(parent_id,4) > 0 )
        
        lineage_data(parent_id, 5) = cc;

    end
    
    end
    
end

% Checking cell cycle has proper start and ending, using track_check data

for cc = 1:cc_num
    
    if ( input.track_check(cc,:) == [1 1]  )
    
        lineage_stat(cc,6) = 1;

    else
        
        lineage_stat(cc,6) = 0;
        
    end

end

end

% ================================================================== 

function [MCstat] = get_MC_stat(mother_cell)

time_unit = 1;
ind = mother_cell.vec(:,1);
cellsize = mother_cell.size(:,1);

cc_list = unique(ind);
cc_num = length(cc_list);

% Have basic statistics for each cell cycle on MC

MCstat = {};
MCstat.cc_data = NaN(cc_num, 7);

% [1] initial time
% [2] ending time
% [3] cell cycle time
% [4] cell cycle GR
% [5] cell cycle GR CoD
% [6] initial size
% [7] added size
% [8] mean size

for c = 1:cc_num
    
    temp = find(ind == cc_list(c));  % a range of time points of this cell cycle    
    MCstat.cc_data(c,1:3) = [temp(1) temp(end) length(temp)];
    
end

for c = 1:cc_num
    
    time = ( MCstat.cc_data(c,1):MCstat.cc_data(c,2) )' ;
    area = cellsize( MCstat.cc_data(c,1) : MCstat.cc_data(c,2) );
    
    [GR, CoD, fit_size] = GR_slope(time, area, time_unit);
    MCstat.cc_data(c,4:5) = [GR CoD];
    
    MCstat.cc_data(c,6:8) = [area(1) area(end)-area(1) mean(area)];
    
end

%

MCstat.mean_ccT = mean(MCstat.cc_data(2:end-1, 3));
MCstat.std_ccT = std(MCstat.cc_data(2:end-1, 3));

MCstat.mean_GR = mean(MCstat.cc_data(2:end-1, 4));
MCstat.std_GR = std(MCstat.cc_data(2:end-1, 4));

end

% ===================================================================

function [GR, CoD, fit_size] = GR_slope(time_point, area, time_unit)

log_area = log(area);

% linear fit on log(area)

coef = polyfit(time_point, log_area, 1);
GR = coef(1)/time_unit;
fit_size = exp(coef(1)*time_point + coef(2));

% Calculating coefficient of determination (CoD, R-square)

log_area_mean = mean(log_area);
log_area_est = polyval(coef, time_point); 

SS_tot = sum( (log_area - log_area_mean).^2 );  % Total sum of square
SS_res = sum( (log_area_est - log_area).^2 );   % Residual sum of square

CoD = 1 - (SS_res/SS_tot);

end

% ===================================================================

