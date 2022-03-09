
% Add linking and cell division data to cell cycles

function [obj_link_list, cc_note, cc_check, cc_ensemble] = ... 
          construct_cell_cycle(separatorC,  sep_filter, segmatN)

% (1) Filter the "good separator tracks"
[separatorG, separator_section] = get_separator(separatorC,  sep_filter, segmatN);

% (2) Link cell object to the next frame. 
[bounding_data] = get_bounding_data(segmatN, separator_section);
[link_rec] = get_link_rec(bounding_data);

% Compile all data into a large table "obj_link_list"
[obj_link_list] = get_obj_link_list(link_rec, bounding_data);

% Convert the large table into struct ('obj_link_struct')
obj_link_struct = {};

row_index = [1:size(obj_link_list,1)]';
obj_link_listA = [obj_link_list row_index];
    
for fn = 1: size(segmatN, 3)   %size(obj_link_list)

    temp = obj_link_listA( (obj_link_listA(:,1) == fn ), :  );    
    obj_link_struct{fn}.data = temp(:,1:8);    
    obj_link_struct{fn}.index = temp(:,9);
    
end

% (3) Construct cell cycle and check the cell cycle quality
[cc_note, obj_link_list] = connect_cell_cycle(obj_link_list, obj_link_struct);
[cc_check] = check_cell_cycle(cc_note);
[cc_ensemble] = ensemble_cell_cycle(obj_link_list, cc_note, cc_check);

end


% ======================================================================

function [separatorG, separator_section] = get_separator(separatorC,  sep_filter, segmatN)

separatorG = {};
flag = 1;

for s = 1:length(separatorC)

    if  ( sep_filter(s) == 1 )  % Good separator
    
        separatorG{flag} = separatorC{s};
        separatorG{flag}.ind = s; 
        flag = flag + 1;
        
    end
    
end
 
separator_section = {};

for fn = 1:size(segmatN,3)

    separator_section{fn} = [];
    
    for s = 1:length(separatorG)
    
        temp = separatorG{s}.data;
    
        if ( fn >= temp(1,1) ) && (fn <= temp(end,1) )
            
            ind = find (temp(:,1) ==  fn);
            separator_section{fn} = [separator_section{fn}; s temp(ind,3)];
            
        end    
    end
end

end

% ======================================================================

function [bounding_data] = get_bounding_data(segmatN, separator_section)

% (1) For each frame, find  "flanking separators" for each cell

bounding_data = {};

for fn = 1:size(segmatN,3)

segdata = segmatN(:,:,fn);

num = max(max(segdata));
bounding_box = zeros(num, 7);

% (i) Finding bounding box

for c = 1:num
    
    temp = (segdata == c);
    yvec = sum(temp,2);
    vec_range = (yvec > 0);

    bounding_box(c,1) = min(find(vec_range == 1));
    bounding_box(c,2) = max(find(vec_range == 1));   
    bounding_box(c,7) = sum(sum(temp)); % cell area
    
end

% (ii) Find the largest lower bound from separator section

for c = 1:num
    
    valL = bounding_box(c,1);    
    sep_loc = separator_section{fn};
    
    sep_loc( (sep_loc(:,2) > valL+1 ), : ) = [];
        
    if ( size(sep_loc,1) > 0 )
        sep_loc = sortrows(sep_loc, 2, 'descend');
        bounding_box(c,3:4) = sep_loc(1,:);
    end
    
end


% (iii) Find the smallest upper bound from separator section

for c = 1:num
    
    valH = bounding_box(c,2);    
    sep_loc = separator_section{fn};
    
    sep_loc( (sep_loc(:,2) < valH-1 ), : ) = [];
        
    if ( size(sep_loc,1) > 0 )
        sep_loc = sortrows(sep_loc, 2, 'ascend');
        bounding_box(c,5:6) = sep_loc(1,:);
    end
    
end

bounding_data{fn}.box = bounding_box;

end


end

% ======================================================================

function [link_rec] = get_link_rec(bounding_data)

link_rec = {};

for fn = 1: (length(bounding_data) - 1)
    
numA = size(bounding_data{fn}.box, 1);
numB = size(bounding_data{fn+1}.box, 1);

    
link_data = [];

for c1 = 1:numA
    
    flanking_separatorA = [bounding_data{fn}.box(c1,3) bounding_data{fn}.box(c1,5)];
    
    for c2 = 1:numB
        
        flanking_separatorB = [bounding_data{fn+1}.box(c2,3) bounding_data{fn+1}.box(c2,5)];
                        
        checkA = sum(flanking_separatorA > 0);
        checkB = sum(flanking_separatorB > 0);
        check_match = sum(flanking_separatorA == flanking_separatorB);        
                
        
        if  ( ( checkA == 2 ) && (checkB == 2) && (check_match == 2) )
            
            % Case 1: Perfect 1-1 linking
            write = [c1 c2 1 flanking_separatorA flanking_separatorB];
            link_data = [link_data; write];            
            
        elseif ( ( checkA == 2 ) && (checkB == 2) && (check_match == 1) )
            
            % Case 2: Cell division event
            write = [c1 c2 2 flanking_separatorA flanking_separatorB];
            link_data = [link_data; write]; 
            
        elseif ( (checkA <= 1) || (checkB <= 1) ) && (check_match == 1)
            
            % Case 3: Incomplete linking (may happen in boundary of the chamber)           
            write = [c1 c2 3 flanking_separatorA flanking_separatorB];
            link_data = [link_data; write]; 
                        
        end  
        
    end
       
end

link_rec{fn}.data = link_data;

end

end

% ======================================================================

function [obj_link_list] = get_obj_link_list(link_rec, bounding_data)

obj_link_list = [];

% [1] frame index
% [2] object index
% [3] object area
% [4] link type
% [5,6] linked post object(s) (if any)
% [7] track label
% [8] mid_y_position

for fn = 1: (length(bounding_data)-1)
           
    obj_data = bounding_data{fn}.box;
    link_data = link_rec{fn}.data; 
    
    write = NaN(size(obj_data,1), 7);
    
    if ( size(link_data,1) > 0  )

    for c = 1:size(obj_data,1)
        
        write(c,1) = fn;        
        write(c,2) = c;
        write(c,3) = obj_data(c,7);        
        write(c,8) = mean(obj_data(c,1:2)); 
 
        link_data_temp = link_data(link_data(:,1) == c, :) ;
        
        if ( size(link_data_temp,1) > 0 )
    
            write(c,4) = link_data_temp(1,3);  % link type    
        
            for r = 1:size(link_data_temp,1)
            
                write(c,4+r) = link_data_temp(r,2);  
            
            end
            
        end
        
    end
        
    obj_link_list = [obj_link_list; write];
    
    end
    
end

end

% ======================================================================

function [row_index] = get_postobj(frame, obj_ind, obj_list)

% Example
%frame_p = 13;
%obj_ind = 3;
%obj_list = obj_full_list;

row_index = 0;

for r = 1:size(obj_list,1)
    
    if ( obj_list(r,1:2) == [frame obj_ind] )
        
        row_index = r;
        break;
    end
        
end

end

% ======================================================================

function [row_index] = get_preobj(frame, obj_ind, obj_list)

% Example:
%frame = 10;
%obj_ind = 3;
%obj_list = obj_full_list;

row_index = 0;

% First, need to find if pre-object exist or not

for r = 1:size(obj_list,1)
    
    Check1 = ( obj_list(r,1) == frame-1 ); 
    Check2a = (obj_list(r,5) == obj_ind);
    Check2b = (obj_list(r,6) == obj_ind);
    Check = Check1 * (Check2a + Check2b);
    
    if ( Check == 1 )     % find the object        
        row_index = r;        
    end
    
end

end

% ======================================================================

function [track_note, obj_link_list] = connect_cell_cycle(obj_link_list, obj_link_struct)

% Note that here "track" means cell cycle tracks, not separators. 

track_note = [];

% Each row in 'track_note' represent a track

% Columns:
% [1] previous track number (p) associated with it
% p == 0   : initial track (tracks start at frame #1)
% p > 0    : previous track label 
% p ==(-1) : did not found associated previous track

% [2] type of ending (q) for this track
% q == 2   : end by cell division
% q == 3   : end by 'did not find any post track'
% q == 4   : end by other error


% (i) Find post-seg and assign same label.
%     otherwise, create new label and repeat similar procedure

track_label = 0;

max_frame = max(obj_link_list(:,1));

for rk = 1 : size(obj_link_list,1)
              
        
        % Test if this object already assiged to a track. 
        
        if ( isnan( obj_link_list(rk,7) ) )             
              
                %%% (a) Create a new track_label
                track_label = track_label + 1;   
                
                % Initiation frame of this track
                frame_ini = obj_link_list(rk, 1);   
                
                %%% (b) Find pre-object for this track. 
                %       Write pre-track type into track_note column[1]
                
                if (frame_ini > 1)
                    
                    current_frame = obj_link_list(rk,1);
                    candidate_data = obj_link_struct{current_frame-1}.data;
                    candidate_index = obj_link_struct{current_frame-1}.index;
                    
                    % Looking for pre-object in the previous frame
                    r_previous_sub = get_preobj( obj_link_list(rk, 1), obj_link_list(rk, 2), candidate_data);

                    if (r_previous_sub > 0)                        
                        r_previous = candidate_index(r_previous_sub);
                    else
                        r_previous = 0;
                    end
                    
                    
                    if (r_previous > 0)  % found the previous track

                        % Add the previous track index into track_note
                        track_note(track_label, 1) = obj_link_list(r_previous, 7);
                        
                    else
                        
                        % Does not found previous track index
                        track_note(track_label, 1) = (-1) ;

                    end
                    
                elseif (frame_ini == 1)
                    
                    % Track start as 1st frame of the dataset
                    track_note(track_label, 1) = 0;
                    
                end
                   
                
                %%% (c) Find post-object
                %       Write post-track type into track_note column[2]
                
                r_next = rk;
                
                for frame_p = frame_ini : (max_frame-1)
                    
                    % (1) label the current object
                    obj_link_list(r_next, 7) = track_label;
                    
                    % (2) find the next simple-link object
                    
                    if ( obj_link_list(r_next, 4) == 1 )  % exist a simple link
                        
                        post_obj_ind = obj_link_list(r_next, 5);
                        
                        current_frame = frame_p;
                        candidate_data = obj_link_struct{current_frame+1}.data;
                        candidate_index = obj_link_struct{current_frame+1}.index;
                                                
                        r_next_sub = get_postobj(frame_p+1, post_obj_ind, candidate_data(:,1:2) );
                        
                        if (r_next_sub > 0 )
                            r_next = candidate_index(r_next_sub);
                        else
                            r_next = 0;
                        end
                        
                        if (r_next == 0)  % cannot find next object
                            
                            track_note(track_label, 2) = 3;
                            break;                            
                        end
                        
                    elseif ( obj_link_list(r_next, 4) == 2 )  % end by cell division

                        track_note(track_label, 2) = 2;                        
                        break;

                    else
                        track_note(track_label, 2) = 4;
                        break;
                    
                    end                    
                    
                end                
                
        end
        
end

end

% ======================================================================

function [track_check] = check_cell_cycle(track_note)

% Note that here "track" means cell cycle tracks, not separators. 

track_num = size(track_note,1);

for r = 1:track_num
    
    % Looking for its offspring 
    
    temp_flag = 3;  
    % column index of writing, start from column 3,
    % can extended to column 4
    
    for sp = 1: track_num
        
        if ( track_note(sp,1) == r )
            
            track_note(r, temp_flag) = sp;
            temp_flag = temp_flag + 1;
        end
        
    end
    
end


track_check = zeros(track_num, 2);

% Column [1]: Check if the beginning type is a proper division
% 1: proper division
% 0: has pre-track, but not proper division (brocken track)
% (-1): no pre-track or track start in frist frame

for c = 1:track_num 
    
    begin_type = NaN;
    
    % Check the pre-track. If not, label (-1) for begin type
    
    if ( track_note(c,1) > 0)  % Found a pre-track
        
        pre_track = track_note(c,1);
        
        if ( track_note(pre_track,2) == 2)  % pre-track ends by cell division 
        
            begin_type = 1;  % proper cell division
            
        else
            
            begin_type = 0;  % not started with cell division
            
        end
            
    else  % does not find a pre-track
        
        begin_type = (-1);
        
    end
    
    track_check(c,1) = begin_type;
    
end

% Column [2]: Check if the ending type is a proper division
% 1: proper division
% 0: has pre-track, no proper division

for c = 1:track_num

    end_type = NaN;
    
    if (track_note(c,2) == 2)  % proper division
        
        end_type = 1;
        
    else 
        
        end_type = 0;
        
    end
    
    track_check(c,2) = end_type;
    
end


end

% ======================================================================

function [cc_ensemble] = ensemble_cell_cycle(obj_link_list, cc_note, cc_check)

% Partition the data matrix (obj_link_list) into cell cycle data structure
% and also add the cell cycle check data.

cc_num = size(cc_note, 1);
cc_ensemble = {};

for m = 1:cc_num  
    
    % Get sub_matrix (temp) of mth track
    temp = obj_link_list;
    temp( temp(:,7) ~= m, :) = [];
    
    cc_ensemble{m}.data = temp; 
    cc_ensemble{m}.cc_check = cc_check(m,:);
    
end

end

