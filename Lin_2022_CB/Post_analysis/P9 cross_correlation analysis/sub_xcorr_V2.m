
function [output] = sub_xcorr_V2(X,Y,tau_vec)

% Perform cross correlation between two time-series data X(t) and Y(t)
% X(t) and Y(t) does not necessarily have same time points

%X = Rsc_data;
%Y = GR_rec;
%tau_vec = [-60:2:60]';

% Using X(t) as reference, find if Y(t+tau) has data. 
% If Y(t+tau) data exist, write into temporally matrix for calculating 
% cross correlation for this specific tau.

all_corr_data = {};
xy_corr_data = NaN(size(tau_vec,1), 1);

for m = 1:size(tau_vec, 1)
    
    tau = tau_vec(m);
    
    for rx = 1:size(X,1)
        
        X(rx,3:4) = [NaN NaN];
        
        % Searching for corss-correlation data        
        ty = X(rx,1) + tau;
        
        for ry = 1: size(Y,1)
        
            if ( Y(ry,1) ==  ty )
                
                X(rx, 3:4) = Y(ry, 1:2);   
            end
            
        end
        
    end
    
    % Remove rows with NaN
    
    X_data = [];
    
    for rx = 1:size(X, 1)
    
        temp = isnan(X(rx,:));

        if (sum(temp) == 0) 
            
            X_data = [X_data; X(rx,:)];
            
        end
        
    end
    
    all_corr_data{m}.tau = tau;
    all_corr_data{m}.X_data = X_data;
        
    % Normalize the data column (z-transform) and calculate the covariance    
    % Find the baseline by fitting the cubic function

    poly_order = 3;
    param_x = polyfit(X_data(:,1), X_data(:,2), poly_order);    
    baseline_x = polyval(param_x, X_data(:,1));
    X_adj = X_data(:,2) - baseline_x;
    
    param_y = polyfit(X_data(:,3), X_data(:,4), poly_order);    
    baseline_y = polyval(param_y, X_data(:,3));
    Y_adj = X_data(:,4) - baseline_y;    
    
    X_norm = ( X_adj - mean(X_adj) ) / std(X_adj);
    Y_norm = ( Y_adj - mean(Y_adj) ) / std(Y_adj);
    
    xy_corr_data(m) = sum(X_norm.*Y_norm)/size(X_norm,1);
   
    
end

output = {};
output.XCF = [tau_vec xy_corr_data];
output.all_data = all_corr_data;


%figure;  plot(tau_vec, xy_corr_data, '-o');
%plot(X_norm, Y_norm, '.');

% figure;
% subplot(211);  plot(X(:,2), 'o');
% subplot(212);  plot(Y(:,2), 'o');
