
% Convert the separator data into kymograph segmentation data structure

function [kseg, sep_filter] = separator_to_kseg_V2(peaks, sep_summary, paramA, paramB, data_dir)

y_range = paramA.Yrange;
t_range = paramA.Nrange;

% (1) Filter for the "good track group" based on the manaul-checked file

cd(data_dir);
sep_manual_check = readmatrix('sep_manual_check.xls');

sep_num = size(sep_manual_check, 1);
sep_filter = zeros(sep_num, 1);
 
for s = 1:sep_num
    
    if ( sep_manual_check(s,5) == 1 )
        
        sep_filter(s) = 1;
    
    end
    
end

%-----------------------------------------------------------------------
% (2a) Filter out separators that are not in (1)

dataF_array = {};

for fn = t_range(1):t_range(2)
    
    data0 = [peaks{fn}.dataC(:,1) peaks{fn}.groupC']; 
    dataF = [];

    for r = 1:size(data0,1)
        
        group_temp = data0(r,2);
            
        if ( sep_filter(group_temp) == 1 ) 
       
            dataF = [dataF; data0(r,:)];   
        
        end        
        
    end
    
    dataF_array{fn}.vec = dataF;
    
end

% (2b) Filter out the "initial portion" of separator
%      to correct time delay from cell constriction to division

dataG_array = {};

for fn = t_range(1):t_range(2)
    
    dataF = dataF_array{fn}.vec;
    dataG = [];
        
    for r = 1:size(dataF,1)
    
    % Checking if this peak belongs to initial portion of the separator
    % if yes, remove this peak from dataF
    
    sep_index = dataF(r,2);
    sep_ini_frame = sep_summary(sep_index, 1);
       
        if (sep_ini_frame == 1)  % start from 1st frame. Do not trimm the separator
        
            dataG = [dataG; dataF(r,:)];        

        elseif (sep_ini_frame > 1)
        
            if ( fn > sep_ini_frame + paramB.sep_delay )       
                dataG = [dataG; dataF(r,:)];
            end
        
        end
    
    end
            
    dataG_array{fn}.vec = dataG;
    
end


%-----------------------------------------------------------------------
% (3a) Generate kseg from separator data

kseg = NaN( y_range(2)-y_range(1)+1, t_range(2)-t_range(1)+1);

for fn = t_range(1):t_range(2)

    dataG = dataG_array{fn}.vec(:,1);
    
    for j = y_range(1):y_range(2)
    
        kseg(j,fn) = sum(j > dataG);
        
    end
        
end

% (3b) For the last group (the one that closest to chamber outlet, 
%      set it to value zero. This removes object in this boundary 

for fn = t_range(1):t_range(2)

    max_ind = max(kseg(:,fn));
    
    for r = 1:size(kseg,1)
        
        if (kseg(r,fn) == max_ind)
            
            kseg(r,fn) = 0;
            
        end
        
    end
        
end

%figure; imagesc(kseg); colormap(jet);

end





