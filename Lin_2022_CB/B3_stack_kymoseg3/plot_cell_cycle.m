
% Plot cell cycle as kymograph. 

function [ret1] = plot_cell_cycle(cc_ensemble, cc_check, save_path)

cc_quality = ( (cc_check(:,1) + cc_check(:,2)) == 2 );

FH1 = figure('position', [1 1 1200 400]);

for c = 1:length(cc_ensemble)

    if ( cc_quality(c) < 1 ) 
        
        temp = cc_ensemble{c}.data;        
        frame = temp(:,1);    
        y_pos = temp(:,8); 
        
        gray_color = 0.9*[1 1 1];
        plot(frame, y_pos, '.', 'color', gray_color); hold on;
    
    elseif (cc_quality(c) == 1)  % a good track, plot the area curve

        temp = cc_ensemble{c}.data;
        frame = temp(:,1);
        y_pos = temp(:,8);
        
        dot_color = rand(3,1);
        plot(frame, y_pos, 'o', 'color', dot_color); hold on;
    
    else        

    end
    
end

xlabel('time (min)');
ylabel('y-position (pixel)');

cd(save_path);
saveas(FH1, 'cc_ensemble_figureN.tiff');

close all;

ret1 = 1;

%}
