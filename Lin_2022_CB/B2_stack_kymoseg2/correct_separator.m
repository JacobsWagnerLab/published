
function [peaks, pair_ind_arrayC, separatorC, peak_linking_tableC0] ...
    = correct_separator(peaks, data_dir, param)

% (0) Load the manually corrected peak data, link into seperator

cd(data_dir);
peak_linking_tableC0 = readmatrix('peak_linking_corrected.xls');
peak_linking_tableC = peak_linking_tableC0(:,1:end-1);  % The last column is modification record. 

for fn = 1:param.max_frame
    
    data_temp = peaks{fn}.data(:,1:2);    
    peaks{fn}.dataC = [data_temp 0*data_temp];
    
end

% (1) Based on the corrected peaks data, write into "peaks" in dataC

for fn = 1: (param.max_frame-1)
    
    temp = peak_linking_tableC(fn,:);
    
    for j = 1:size(temp,2)
        
        if ( temp(j) > 0  )
    
            pre_ind = j;
            post_ind = temp(j);
            
            peaks{fn}.dataC(pre_ind,4) = post_ind;
            peaks{fn+1}.dataC(post_ind, 3) = pre_ind;
        
        end
        
    end
    
end

tag = 'dataC';
[separatorC, terminal_marks, peaks] = partition_separators(peaks, tag, param.max_frame);


% (2) Based on the corrected peaks data, write into "pair_ind_arrayC"

pair_ind_arrayC = {};

for fn = 1: (param.max_frame-1)
    
    temp = peaks{fn}.dataC;
    list = [];
    
    for j = 1:size(temp,1)
        
        if ( temp(j,4) > 0 )  % connected to some post-peak
        
            list = [list; j temp(j,4) peaks{fn}.groupC(j)];
            
        end
        
    end
    
    pair_ind_arrayC{fn}.data = list;
    
end


% ======================================================================

% (3) Filter the separator based on its ending 
% (need either hittting the last frame of being closed to exit position)   

% (i) Test if the separator hit last frame
separator_num = length(separatorC);

for m = 1:separator_num
        
    if ( separatorC{m}.data(end,1) == param.max_frame )

        separatorC{m}.ending = 1;
        
    else
        
        separatorC{m}.ending = 0;
    
    end        

end

% (ii) Test if the separator ending is closed to the exit position

for m = 1:separator_num
        
    if ( separatorC{m}.data(end,3) > param.exit_y )

        separatorC{m}.exit_y = 1;
        
    else 
        
        separatorC{m}.exit_y = 0;
    
    end        

end

% (iii) Test if the separator length is longer than minimal length

for m = 1:separator_num
        
    if ( size( separatorC{m}.data, 1 )  > param.minimal_length )

        separatorC{m}.length = 1;
        
    else
        
        separatorC{m}.length = 0;
        
    end        

end

% ================================================================== %
% (4) Export the separator plot

plot_flag = param.plot_flag;

if (plot_flag == 1)

hfig = figure('position', [1 1 1500 800]);

colorvecA = 0.9*rand(separator_num, 3);
colorvecB = colorvecA;

% For display purpose, set the "bad separator" to transparent. 

for m = 1:separator_num

    if (separatorC{m}.ending == 0) && (separatorC{m}.exit_y == 0)        
        
        colorvecB(m,:) = [1 1 1];

    end
    
    if (separatorC{m}.length == 0)
        
        colorvecA(m,:) = [1 1 1];
        colorvecB(m,:) = [1 1 1];
        
    end
    
end

for m = 1:separator_num
    
    separatorC{m}.colorA = colorvecA(m,:);
    separatorC{m}.colorB = colorvecB(m,:);
    
end

% --------------------------------------------------------------------
% (1) Plot all separators

subplot(211);

for m = 1:length(separatorC)
    
    data = separatorC{m}.data;    
    plot(data(:,1), data(:,3), 'color',  colorvecA(m,:));  hold on;
        
end

% (2) Plot good separators 

subplot(212);

for m = 1:length(separatorC)
    
    data = separatorC{m}.data;    
    plot(data(:,1), data(:,3), 'color',  colorvecB(m,:));  hold on;
        
end


cd(data_dir);
saveas(hfig, 'separator_corrected.jpg');


end


close all;


end


% =======================================================================

function [separatorC, terminal_marks, peaks] = partition_separators(peaks, tag, max_frame)

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

separator_num = size(terminal_marks, 1);
separatorC = {};

for g = 1:separator_num
    
    % separatoring all separator marks from a initial node
    
    ini_frame = terminal_marks(g,1);
    ini_ind = terminal_marks(g,2);
    ypos = peaks{ini_frame}.(tag);
    ypos = ypos(ini_ind,1);
    
 %   [g ini_frame ini_ind]
 
    ini_node = [ini_frame ini_ind ypos];  % frame, ind, ypos

    current_node = ini_node;
    separator_data = [];
    separator_data = [separator_data; current_node];

    for fn = ini_frame : (max_frame - 1)
    
        next_frame = fn+1;
        next_ind = peaks{fn}.(tag)(current_node(2), 4);    
    
        if (next_ind > 0)  % separator continues
        
            next_ypos = peaks{next_frame}.(tag)(next_ind,1);
        
            current_node = [next_frame next_ind next_ypos];
            separator_data = [separator_data; current_node];
        
        elseif (next_ind == 0) % separator stopped
        
            break;
        
        end
    
    end
    
    separatorC{g}.data = separator_data;    
    
end

% (3) Back to separator and label the separator group 

for g = 1:separator_num
    
    for m = 1:size(separatorC{g}.data, 1)
        
        fn = separatorC{g}.data(m,1);
        ind = separatorC{g}.data(m,2);
        
        peaks{fn}.groupC(ind) = g;        
        
    end
    
end


end

% =======================================================================
