
% Append c2,c3 intensity data for each cell cycle

function [cc_ensemble] = append_cc_intensity(cc_ensemble, fluo_rec)

for c = 1 : length(cc_ensemble)
    
    data = cc_ensemble{c}.data;
    
    signal = [];
    
    for j = 1:size(data,1)
        
        frame = data(j,1);
        cellID = data(j,2);
    
        [signal_data] = get_signal(frame, cellID, fluo_rec);
        
        if ( ~isempty(signal_data) )
            
            write = [frame cellID signal_data];
            signal = [signal; write];
            
        end
        
    end
    
    cc_ensemble{c}.signal = signal; 
    
end

end

% ======================================================================= 

function [signal_data] = get_signal(fn, cellID, fluo_rec)

% Example: 
% data = cc_ensemble{1}.data;
% frame = 1;
% cellID = 1;

intensity_col = 4;  % Intensity at 4th column
cellID_col = 6;     % CellID    at 6th column

signal_data = [];

if (fn <= length(fluo_rec))
   
if ( ~isempty(fluo_rec{fn}) )
    
	signal_fn = fluo_rec{fn};
        
	c2_data_row = signal_fn.c2( (signal_fn.c2(:,cellID_col) == cellID), :  );         
    c2_signal = c2_data_row(intensity_col);
        
    c3_data_row = signal_fn.c3( (signal_fn.c3(:,cellID_col) == cellID), :  );         
    c3_signal = c3_data_row(intensity_col);
        
    signal_data = [c2_signal c3_signal];
                
end

end



end
