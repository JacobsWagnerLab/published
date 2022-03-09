
function [peaks, separator] = peak_to_separator(peaks, save_dir, param)

% ======================================================================

% (A) Find separators as tracks.

tag = 'data';
[separator, terminal_marks, peaks] = partition_separators(peaks, tag, param);

% ======================================================================

% (B) Filter the track based on its ending 
% (need either hittting the last frame of being closed to exit position)   

% (i) Test if the track hit last frame

plot_flag = param.plot_flag;

sep_num = length(separator);

for m = 1:sep_num
        
    if ( separator{m}.data(end,1)  == param.max_frame )

        separator{m}.ending = 1;
        
    else
        
        separator{m}.ending = 0;
    
    end        

end

% (ii) Test if the track ending is closed to the exit position

for m = 1:sep_num
        
    if ( separator{m}.data(end,3) > param.exit_y )

        separator{m}.exit_y = 1;
        
    else 
        
        separator{m}.exit_y = 0;
    
    end        

end

% (iii) Test if the track length is longer than minimal length

for m = 1:sep_num
        
    if ( size( separator{m}.data, 1 )  > param.minimal_length )

        separator{m}.length = 1;
        
    else
        
        separator{m}.length = 0;
        
    end        

end

% ================================================================== %

% (C) Plot separator

if (plot_flag == 1)

hfig = figure('position', [1 1 1500 800]);

colorvecA = 0.9*rand(sep_num, 3);
colorvecB = colorvecA;

% For display purpose, set the "bad track" to transparent. 

for m = 1:sep_num

    if (separator{m}.length == 0)
        
        colorvecA(m,:) = [1 1 1];
        colorvecB(m,:) = [1 1 1];
    end
    
    if (separator{m}.ending == 0) && (separator{m}.exit_y == 0)        
        
        colorvecB(m,:) = [1 1 1];
        
    end
    
end

for m = 1:sep_num
    
    separator{m}.colorA = colorvecA(m,:);
    separator{m}.colorB = colorvecB(m,:);
    
end


% --------------------------------------------------------------------
% (1) Plot all tracks

subplot(211);

for m = 1:length(separator)
    
    data = separator{m}.data;    
    plot(data(:,1), data(:,3), 'color',  colorvecA(m,:));  hold on;
        
end

% (2) Plot good tracks 

subplot(212);

for m = 1:length(separator)
    
    data = separator{m}.data;    
    plot(data(:,1), data(:,3), 'color',  colorvecB(m,:));  hold on;
        
end


cd(save_dir);
saveas(hfig, 'separator.jpg');
close all;


end



end


% =======================================================================

function [separator, terminal_marks, peaks] = partition_separators(peaks, tag, param)

max_frame = param.max_frame;

% (1) Finding terminal marks (peaks that has no pre-link) in the dataset

terminal_marks = [];

for fn = 1:max_frame
    
    for j = 1:size(peaks{fn}.(tag), 1)        
        
        if ( peaks{fn}.(tag)(j,3) == 0 )  % no prelink 
            
            terminal_marks  = [terminal_marks; fn j];
            
        end
        
    end
    
end

% (2) For each terminal marks, try to find the post-link repeatly. 
% Each terminal marks is extended into a "separator" across frames

sep_num = size(terminal_marks, 1);
separator = {};

for g = 1:sep_num
    
    % tracking all separator marks from a initial node
    
    ini_frame = terminal_marks(g,1);
    ini_ind = terminal_marks(g,2);
    ypos = peaks{ini_frame}.(tag);
    ypos = ypos(ini_ind,1);
    
 %   [g ini_frame ini_ind]
 
    ini_node = [ini_frame ini_ind ypos];  % frame, ind, ypos

    current_node = ini_node;
    track_data = [];
    track_data = [track_data; current_node];

    
    for fn = ini_node(1) : (max_frame - 1)
    
        next_frame = fn+1;
        next_ind = peaks{fn}.(tag)(current_node(2), 4);    
    
        if (next_ind > 0)  % track continues
        
            next_ypos = peaks{next_frame}.(tag)(next_ind,1);
        
            current_node = [next_frame next_ind next_ypos];
            track_data = [track_data; current_node];
        
        elseif (next_ind == 0) % track stopped
        
            break;
        
        end
    
    end
    
    separator{g}.data = track_data;    
    
end

% (3) Back to separator and label the track group 

for g = 1:sep_num
    
    for m = 1:size(separator{g}.data, 1)
        
        fn = separator{g}.data(m,1);
        ind = separator{g}.data(m,2);
        
        peaks{fn}.group(ind) = g;        
        
    end
    
end

end

% =======================================================================