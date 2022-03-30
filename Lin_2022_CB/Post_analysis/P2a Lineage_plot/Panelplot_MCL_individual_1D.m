
% Plot a panel of Rsc figures with cell division and media flag

path = {};
path.load = '/Volumes/Data_04/Wei-Hsiang Lin/uF_dataset/stationary/MCL/';
path.dataset_name_prefix = 'MCL_combined_';
path.dataset_name = 'M9only';

x_range = [0 1500];
y_range1 = [0.5 3];
y_range2 = [1 3];
media_switch = [];     %Example: media_switch = [1 960 1920 2880];

pixel_to_um2 = 1/272;  % convert pixel to um^2 

% ======================================================================

index_list = [8];
    
cd(path.load);    
input = load( strcat(path.dataset_name_prefix, path.dataset_name) ).lineage_data;
    
for fm = 1:length(index_list)
        
    m = index_list(fm);    
    data = input{m};
          
    MCLdata = {};
        
    MCLdata.size = data.MC.size* pixel_to_um2;  % convert pixel to um^2 
    MCLdata.Rsc = data.MC.Rsc.data;  
    MCLdata.flag = data.MC.flag;
        
	% (1) Constructing media flag vector
        
	media_state = NaN(length(MCLdata.size),1);
        
	if ( ~isempty(media_switch) ) % Enter media flag
            
        flag_state = -1;
        flag_count = 1;        
        
        for j = 1:length(MCLdata.size)
            
            if ( j == media_switch(flag_count) )
                flag_state = (-1)*flag_state;
                flag_count = flag_count + 1;                
            end
                
            media_state(j) = flag_state;         

        end
            
    end
        
        
	% (2) Plot figures 
	hfig = figure('position', [1 1 1000 300]);
       
	% (2-1) Cell size in log scale
    
	subplot(2,1,1);      
    
    time = (1:size(MCLdata.size,1))' ;
    semilogy(time, MCLdata.size, '.b', 'MarkerSize', 4); hold on;
    semilogy(y_range1(2)*MCLdata.flag+0.1, 'r-'); hold on;
    
    if ( ~isempty(media_switch) )          
        semilogy(y_range1(2) *(media_state+1)+1, 'g-'); 
    end
    hold off;                
    
        
    xlim(x_range);        
    ylim(y_range1);       
    xlabel('time (min)');   
    ylabel('cell size (\mu m2)');
    
    title_info = strcat(data.dataset_index, ',CB', num2str(data.chamber_index));   
    title(title_info);
    
    % (2) R-score
    
	subplot(2,1,2);
	plot(MCLdata.Rsc(:,1), MCLdata.Rsc(:,2), '-k.');  hold on;   
	plot( (y_range2(2)+1) *MCLdata.flag, 'r-'); hold on;
    
	if ( ~isempty(media_switch) )        
        plot(y_range2(2) *media_state+1, 'g-');    
    end        
	hold off;
    
    xlim(x_range);   
    ylim(y_range2);   
    xlabel('time (min)');       
    ylabel('R-score ');

end

% ======================================================================

