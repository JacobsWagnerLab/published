
% Connect peak data of adjacent columns into separators

% pair_ind_array{fn}.ind    : index pair for peaks 
% pair_ind_array{fn}.indW   : index pair for weak peaks

% peaks{fn}.data: 
% col 1: y-position
% col 2: peak intensity
% col 3: index of linked peaks in fn
% col 4: index of linked peaks in fn+1
% col 5: index of linked weak peaks in fn
% col 6: index of linked weak peaks in fn+1

function [pair_ind_array, peaks] = connect_peaks(peaks, param)

fn_range = param.Nrange;   % Range for all frames

pair_ind_array = {};

% (1) Connect peaks from consecutive frames

for fn = fn_range(1):fn_range(2)-1
    
    if ( ~isempty(peaks{fn}.data) ) && ( ~isempty(peaks{fn+1}.data) )
    
    pair_ind = align_peaks_between_frames(peaks, fn, param);
    pair_ind_array{fn}.ind = pair_ind;
    
    for j = 1:size(pair_ind, 1)
        
        pre_ind = pair_ind(j,1);
        post_ind = pair_ind(j,2);
        
        peaks{fn}.data(pre_ind, 4) = post_ind;   
        peaks{fn+1}.data(post_ind, 3) = pre_ind;
        
    end
    
    end
    
end

% (2) For remaining unaligned peaks, use a more relaxed criteria
%     Try to connect the peaks that are mutually closet to each other
%     in metric space

for fn = fn_range(1):fn_range(2)-1
    
    if ( ~isempty(peaks{fn}.data) ) && ( ~isempty(peaks{fn+1}.data) )
        
    pair_indW = align_peaks_between_framesW(peaks, pair_ind_array, fn, param);
    pair_ind_array{fn}.indW = pair_indW;    
    
    for j = 1:size(pair_indW, 1)
        
        pre_indW = pair_indW(j,1);
        post_indW = pair_indW(j,2);        
        
        peaks{fn}.data(pre_indW, 4) = post_indW;   
        peaks{fn+1}.data(post_indW, 3) = pre_indW;
        
        peaks{fn}.data(pre_indW, 6) = post_indW;   
        peaks{fn+1}.data(post_indW, 5) = pre_indW;
                        
    end
    
    end
    
end

end


% ======================================================================

function [pair_indW] = align_peaks_between_framesW(peaks, pair_ind_array, fn, param)

% Align for "weak peak pairs" with more relaxed parameters


% ---------------------------------------------------------------

% (1) Make copies of peak data of fn and fn+1 into peakC and peakD;
% remove the linked peaks.

peakC0 = peaks{fn}.data(:,1:4);
peakC0(:,5) = [1:size(peakC0,1)]';

peakD0 = peaks{fn+1}.data(:,1:4);
peakD0(:,5) = [1:size(peakD0,1)]';

peakC = peakC0;
peakD = peakD0;

peakC( (peakC0(:,4) > 0), : ) = [] ;
peakD( (peakD0(:,3) > 0), : ) = [] ;

% (2) For the remaining unlinked peaks, use weaker criteria to link them

pair_indW = [];  % pair index for "weak pairs"

n1 = size(peakC,1);
n2 = size(peakD,1);

if ( n1 > 0 ) && ( n2 > 0 )     % Has unlinked peak
    
    % (2-1) Calculate the distance matrix
    dist_mat = NaN(n1, n2);
    
    for j1 = 1:n1
        
        for j2 = 1:n2
            
            d1 = param.dist_weight(1) * abs(peakC(j1,1) - peakD(j2,1)); 
            d2 = param.dist_weight(2) * abs(peakC(j1,2) - peakD(j2,2));  
        
            dist_mat(j1,j2) = d1 + d2; 
                        
        end
        
    end   
    
    % (2-2) For each pair, test if the distance are mutually minimized. 

    candidateC = NaN(n1,2);  % distance, absolute index
    candidateD = NaN(n2,2); 

    for j1 = 1:n1
    
        temp = dist_mat(j1,:);
        [val, relative_ind] = min(temp);

        candidateC(j1,:) = [val relative_ind];
    
    end

    for j2 = 1:n2
    
        temp = dist_mat(:,j2);
        [val, relative_ind] = min(temp);

        candidateD(j2,:) = [val relative_ind];
    
    end
    
    candidate_pair = [];
    
    for j1 = 1:n1
    
        candidate_ind_post = candidateC(j1,2);
        candidate_ind_pre = candidateD(candidate_ind_post, 2);
    
        if ( j1 == candidate_ind_pre ) % mapped back to itself
        
            if ( dist_mat(j1, candidate_ind_post) < param.dist_cutW )
            
                absolute_ind1 = peakC(j1,5);
                absolute_ind2 = peakD(candidate_ind_post,5);             
               
                write = [absolute_ind1 absolute_ind2];
                candidate_pair = [candidate_pair; write] ; 
                
            end
        
        end
        
    end
    
    % (2-3) Test that the newly added pair cannot cross with pre-exist pair
    
    for k = 1:size(candidate_pair)
        
        yboundC = peakC0(candidate_pair(k,1), 1);
        yboundD = peakD0(candidate_pair(k,2), 1);
        
        ybound_min = min(yboundC, yboundD);
        ybound_max = max(yboundC, yboundD);
        
        cross_flag = 0;
        
        indP = pair_ind_array{fn}.ind;
        
        for pr = 1:size(indP, 1)
            
            yA = peaks{fn}.data (indP(pr,1), 1);
            yB = peaks{fn+1}.data (indP(pr,2), 1);
            
            %[yA yB]
            
            ymax = max(yA, yB);
            ymin = min(yA, yB);
            
            % Crossing criteria:
            
            if ( (ymax < ybound_min) || (ymin > ybound_max) )
                
                % No crossing. 
                cross_flag = cross_flag + 1;
                
            end
            
        end
        
        if ( cross_flag == size(indP,1) )
                
             pair_indW = [pair_indW; candidate_pair(k,:)];
            
        end
        
        
    end
    
    
        
end


end

% ================================================================

function [pair_ind] = align_peaks_between_frames(peaks, fn, param)

% Align peak pairs between frame fn and fn+1.

% Example for parameters
%fn = 274;
%peakA= peaks{fn}.data;    
%peakB = peaks{fn+1}.data;
%col1 = kmat(:,fn);
%col2 = kmat(:,fn+1);

% ---------------------------------------------------------------
% (1) Start from small y position (chamber end) 
%     connect pairs of peaks between two frames

peakA = peaks{fn}.data;
peakB = peaks{fn+1}.data;

peak_num = [size(peakA,1) size(peakB,1)];
max_pair = max(peak_num);

tempA = peakA;  % data vector, with y-positions of peaks of frame fn
tempB = peakB;  % data vector, with y-positions of peaks of frame fn+1

rec = {};

for j = 1:max_pair
    
    % If there are unaligned peak remains, find more next peak pairs
    
    if (size(tempA,1) > 0) && (size(tempB,1) > 0)
        
        [pair_data, tempA, tempB] = get_one_peak_pair(tempA, tempB, param);    
        rec{j} = pair_data;

    end
    
end


% (2) Convert relative index into absolute index of peak

pair_relative_ind = [];

for j = 1:length(rec)    
    pair_relative_ind(j,:) = rec{j}.ind; % Relative index returned from each round of get_one_peak_pair    
end

pair_ind = [];  % Absolute index of peak pairs
write = [0 0];

for j = 1:size(pair_relative_ind,1)
    
    temp = pair_relative_ind(j,:);
    
    if ( ~isnan(temp(1)) ) && (j == 1)
        
        pair_ind = [pair_ind; temp];
    
    elseif ( ~isnan(temp(1)) ) && (j > 1)
        
        if ( size(pair_ind,1) > 0 )
        
            write = write + pair_ind(end,:) + temp;
            pair_ind = [pair_ind; write];
            write = [0 0];
        
        else
            
            write = write + temp;
            pair_ind = [pair_ind; write];
            write = [0 0];
        
        end
        
    elseif ( isnan(temp(1)) ) && (j == 1)
        
        write = [1 1];
        
    elseif ( isnan(temp(1)) ) && (j > 1)
        
        write = write + [1 1];

    end


end


end

% ================================================================

function [pair_data, tempAe, tempBe] = get_one_peak_pair(tempA, tempB, param)

% For all peak pairs between two frames, find the closet peak pair into "pair_data"
% and return the remaining unaligned peaks data as tempAe and tempB2 (with relative y-position)

% (1) For all peaks pairs, calculate the distance matrix

dist = NaN(size(tempA,1), size(tempB,1));  % distance matrix

for j1 = 1:size(tempA,1)
    
    for j2 = 1:size(tempB,1)
    
    	d1 = param.dist_weight(1) * abs(tempA(j1,1) - tempB(j2,1)); 
        d2 = param.dist_weight(2) * abs(tempA(j1,2) - tempB(j2,2));  
        
        dist(j1,j2) = d1 + d2; 
        
    end
    
end

% (2) Looking for peak pairs from the smallest y-index

min_val = 10000; 
min_ind = [NaN NaN];

for j1 = 1:size(tempA,1) 

    if (dist(j1,1) < min_val)
    
        min_val = dist(j1,1);
        min_ind = [j1 1];
        
    end    
end

for j2 = 1:size(tempB,1) 

    if (dist(1,j2) < min_val)
    
        min_val = dist(1,j2);
        min_ind = [1 j2];
        
    end
end

% (3) Remove the index smaller than min_ind

tempAe = tempA; 
tempBe = tempB;

if ( min_val < param.dist_cut )  % Find an accetable pairs
        
    tempAe(1:min_ind(1),:) = [];   
    tempBe(1:min_ind(2),:) = [];    

    % (i) calculate the relative y-postion for remaining peaks
    y_offsetA = tempA(min_ind(1),1);
    y_offsetB = tempB(min_ind(2),1);    
    
    tempAe(:,1) = tempAe(:,1) - y_offsetA;  
    tempBe(:,1) = tempBe(:,1) - y_offsetB; 
    
    % (ii) return the peak pair data
    pair_data.ind = min_ind;
    pair_data.dist = min_val;

else    % No peak pairs found
    
    tempAe(1,:) = [];  
    tempBe(1,:) = [];    
   
    % (i) Reset relative y_pos with respect to the first peak of fn
    y_offsetA = tempA(1,1);
    y_offsetB = tempA(1,1);

    tempAe(:,1) = tempAe(:,1) - y_offsetA;
    tempBe(:,1) = tempBe(:,1) - y_offsetB;
        
    % (ii) return the peak pair data
    pair_data.ind = [NaN NaN];
    pair_data.dist = NaN;
    
end


% (4) Record the pair data

pair_data.tempA = tempA;  
pair_data.tempB = tempB;
pair_data.dist = dist;

end

% ================================================================

