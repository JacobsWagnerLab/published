
function [output] = PS_analysis_summary(MCL_data, tag1, tag2, tag3)

%tag1 = 'Rsc_data';
%tag2 = 'PS';
%tag3 = 'PS' or 'SH_mean';

output = {};
output.tag1 = tag1;
output.tag2 = tag2;
output.tag3 = tag3;

% (2) Fit all power spectrum data curves on the same frequency scale

PS_array = {};
frequency_scale = (0:0.05:5)';  % in (1/hr)

for m = 1:length(MCL_data)
    
    temp = MCL_data{m}.(tag1).(tag2);
    
    if ( ~isempty(temp) )
    
        freq = 60*(temp.freq);  % per hour
        spec = temp.(tag3);
        
        PS_fit.spec = interp1(freq, spec, frequency_scale);
        PS_fit.freq = frequency_scale;
        
        PS_array{m}.PS_fit = PS_fit;
        
    end        
        
end

% (3) Organize into table format for all power spectrum

PS_table = {};
PS_table.freq = frequency_scale';
PS_table.PSdata = [];
PS_table.info = [];


for m = 1:length(PS_array)

    temp = PS_array{m};
    
    if ( ~isempty(temp) )
        
        PS_data_temp = PS_array{m}.PS_fit.spec';            
        PS_table.PSdata = [PS_table.PSdata; PS_data_temp];
        
        PS_table.info = [PS_table.info; MCL_data{m}.dataset_indnum MCL_data{m}.chamber_index];

    end
    
end

PS_table.PSstat.mean = nanmean(PS_table.PSdata,1);
PS_table.PSstat.std = nanstd(PS_table.PSdata,1);
PS_table.PSstat.se = PS_table.PSstat.std / sqrt(size(PS_table.PSdata,1));


% (4) Organize into table format for different FOV

PS_FOV = {};
group_num = MCL_data{end}.dataset_indnum;

for g = 1:group_num
    
    PSdata_g = PS_table.PSdata;
    PSdata_g( (PS_table.info(:,1) ~= g), : ) = [];
    
    PS_FOV{g}.data = PSdata_g;
    
end

for g = 1:group_num
    
    PS_FOV{g}.freq = PS_table.freq;
    PS_FOV{g}.stat.mean = mean(PS_FOV{g}.data,1);
    PS_FOV{g}.stat.std = std(PS_FOV{g}.data,1);
    PS_FOV{g}.stat.se = std(PS_FOV{g}.data,1) / sqrt(size(PS_FOV{g}.data,1));
    
end

output.PS_array = PS_array;
output.PS_table = PS_table;
output.PS_FOV = PS_FOV;

end

