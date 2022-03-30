
% Manaully load 'tempogram_summary'

ind_list = [1]';
ind_num = size(ind_list, 1);

% ====================================================================
% Caculate averaged curve for each matR, matZ
% ====================================================================

% (0) Creating RYB colormap 
               
RYBmap = customcolormap(linspace(0,1,11), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd','#def4f9','#abd9e9','#73add2','#4873b5','#313691'}, 128);

% ====================================================================
% 
% (1) Plot figures (Z-transformed)
% 
% h1 = figure('position', [1 1 1000 400]);
% 
% for b = 1:ind_num
%     
%     matZ = tempogram_summary{ind_list(b)}.data.sorted_matZ;    
%     subplot(1,ind_num,b);
%     imagesc(matZ, [-2.5 2.5]);
%     axis off;
%         
%     colorbar('location', 'northoutside');
% 
% end
% 
% colormap(RYBmap);


% (2) Plot figures (original)

h2 = figure('position', [1 1 350 400]);

b = 1;
    
matR1 = tempogram_summary{ind_list(b)}.data.sorted_matR_mean;    
matR2 = tempogram_summary{ind_list(b)}.data.sorted_matR_SD;    
matR0 = tempogram_summary{ind_list(b)}.data.sorted_matR; 
matZ0 = tempogram_summary{ind_list(b)}.data.sorted_matZ; 

subplot(141);
imagesc(matR1, [3 7]);
axis off;            
%colorbar('location', 'westoutside');
colormap(RYBmap);

subplot(142);
imagesc(matR2, [3 7]);
axis off;  

subplot(143);
imagesc(matR0, [3 7]);
axis off;  

subplot(144);
imagesc(matZ0, [-2.5 2.5]);
axis off;  

% (3) Plot min and max distribution
% 
% figure;
% 
% for b = 1:ind_num
%     
%     matZ = tempogram_summary{ind_list(b)}.data.sorted_matR; 
%     [output] = find_min_max_dist(matZ);
%     
%     subplot(1,ind_num,b);    
%     bar(output.max_dist);
%     ylim([0 0.3]);
%     
% %    subplot(2,ind_num,ind_num+b);
% %    bar(output.min_dist);
% %    ylim([0 0.3]);
%     
%     
% end

% 
% % (4) Plot average curve
% 
% figure;
% 
% for b = 1:ind_num
%     
%     tempR = tempogram_summary{ind_list(b)}.data.sorted_matR;
%     avg_curveR = mean(tempR,1);
%     
%     tempZ = tempogram_summary{ind_list(b)}.data.sorted_matZ;
%     avg_curveZ = mean(tempZ,1);
%     
%     subplot(2,ind_num,b); 
%     plot(avg_curveR, 'ko-');
%     ylim([4.25 5.5]);
%     
%     subplot(2,ind_num,b+ind_num); 
%     plot(avg_curveZ, 'ko-');    
%     ylim([-0.5 0.5]);
%     
% end
% 



% ====================================================================

function [output] = find_min_max_dist(matZ)

bin_size = 10;
bin_num = 100/bin_size;
cc_num = size(matZ, 1);

ATP_bin_average = NaN(bin_num, cc_num);

for B = 1:bin_num
    
    bin_ind_temp = [(B-1)*bin_size+1  B*bin_size]; 
    matZ_temp = matZ(:, bin_ind_temp(1):bin_ind_temp(2));
    ATP_bin_average(B,:) = mean(matZ_temp,2);     

end

ATP_min_max_bin = NaN(cc_num, 2);

for cc = 1:size(ATP_bin_average,2)
    
    [val1, ind1] = max(ATP_bin_average(:,cc));
    ATP_min_max_bin(cc,1) = ind1;
    
    [val2, ind2] = min(ATP_bin_average(:,cc));
    ATP_min_max_bin(cc,2) = ind2;
    
end

ATP_min_max_hist = zeros(bin_num, 2);

for cc = 1:cc_num 
    
    ind1 = ATP_min_max_bin(cc,1);
    ATP_min_max_hist(ind1, 1) = ATP_min_max_hist(ind1, 1) + 1;
    
    ind2 = ATP_min_max_bin(cc,2);
    ATP_min_max_hist(ind2, 2) = ATP_min_max_hist(ind2, 1) + 1;    
    
end

ATP_max_dist = ATP_min_max_hist(:,1)/sum(ATP_min_max_hist(:,1));
ATP_min_dist = ATP_min_max_hist(:,2)/sum(ATP_min_max_hist(:,2));

output = {};
output.max_dist = ATP_max_dist;
output.min_dist = ATP_min_dist;

%figure;
%plot(ATP_max_dist, 'ro-'); hold on;
%plot(ATP_min_dist, 'bo-'); hold off;

end

% ====================================================================
% 
% figure('position', [1 1 1000 100]);
% 
% Zcutoff = 1;
% 
% for b = 1:ind_num
%     
%     matZ = tempogram_summary{ind_list(b)}.data.sorted_matZ;  
%     
%     signature1 = 100*(matZ > Zcutoff);
%     signature2 = 100*(matZ < (-1)*Zcutoff);
%     
%     matZ_mean1 = mean(signature1, 1);
%     matZ_std1 = std(matZ, 0, 1);
%     
%     matZ_mean2 = mean(signature2, 1);
%     matZ_std2 = std(matZ, 0, 1);
%     
%     matZ_meanD = matZ_mean1-matZ_mean2;
%     
%     subplot(2,ind_num,b);
% 
%     plot(matZ_mean1, 'r-'); hold on;
%     plot(matZ_mean2, 'b-'); hold on;    
%     ylim([0 40]);
% 
%     subplot(2,ind_num,8+b);
%     plot(matZ_meanD, 'k-'); 
%     ylim([-30 30]);
% 
% end
