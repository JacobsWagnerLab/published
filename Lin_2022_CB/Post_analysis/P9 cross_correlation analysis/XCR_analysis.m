
% This is the code for cross-correlation analysis
% Need data from lineage analysis (under /MCL/ folder)

path = {};
path.load = '/Volumes/Data_04/Wei-Hsiang Lin/uF_dataset/exponential/MCL/';

dataset_nametag = {'M9GlcCA_lac', 'M9GlcCA_ldhA', 'M9GlcCA_adhE', 'M9GlcCA_pta', 'M9GlcCA_ackA',...
                   'M9Glc_lac', 'M9Xyl_lac', 'M9Gly_lac'};

DB_num = length(dataset_nametag);


% =====================================================================
% Parameter for XCR analysis:

w = 15;  % Using 15 minutes as sliding window of instantaneous growth rate
tau_vec = [-200:1:200]';  % window for calculating XCR

% =====================================================================

XCR_data = {};

for DB = 1:DB_num

    cd(path.load);
    dataset_name = strcat('MCL_combined_', dataset_nametag{DB}, '.mat');   
    temp = load(dataset_name).lineage_data;
    num = length(temp);  
    
    for m = 1:num
        
        cell_area = temp{m}.MC.size;
        cc_flag = temp{m}.MC.flag;
        
        % (1) Get instantenous growth rate
        GR_t = instant_GR_time_series(cell_area, cc_flag, w);       

        % (2) Get Rsc and sensor trajectory (excluded the last cell cycle)
        Rsc_t = temp{m}.MC.Rsc.dataS;
        sensor_t = temp{m}.MC.sensor.dataS;
        
        % (3) Calculate XCR        
        X = Rsc_t;
        Y = GR_t;
        Rsc_GR_xcorr = sub_xcorr_V2(X,Y,tau_vec);
        
        X = sensor_t;
        Y = GR_t;
        sensor_GR_xcorr = sub_xcorr_V2(X,Y,tau_vec);
        
        % =========== Record data =============

        XCR_data{DB}.data{m}.GR_t = GR_t;
        XCR_data{DB}.data{m}.Rsc_t = Rsc_t;
        XCR_data{DB}.data{m}.Rsc_GR_xcorr = Rsc_GR_xcorr;
        XCR_data{DB}.data{m}.sensor_GR_xcorr = sensor_GR_xcorr;        
        XCR_data{DB}.dataset_name = dataset_nametag{DB};
        
        XCR_data{DB}.param.w = w;
        XCR_data{DB}.param.tau_vec = tau_vec;
        
    end
    
end


% =====================================================================

function [GRt] = instant_GR_time_series(cell_area, cc_flag, w)

% Find suitable timepoints for sliding window
% (1) If the window overlapped with two cell cycle, return NaN
% (2) If the window is at the terminus of time-series, return NaN.

L = length(cell_area);
test_array = NaN(L,2);

for t = 1:L
    
    t1 = t - floor(w/2);
    t2 = t + floor(w/2);
        
    test1 = ( (t1 >= 1) && (t2 <= L) );
    test2 = 1;
    
    if (test1 == 1)
        
        cell_area_sub = cell_area(t1:t2);
        cc_flag_sub = cc_flag(t1:t2);
        
        flag0 = cc_flag_sub(1);

        for j = 1:length(cc_flag_sub)
            
            if ( cc_flag_sub(j) ~= flag0 )
                test2 = test2 * 0;
            end
            
        end
        
    end
    
    test_array(t,:) = [test1 test2];
    
end

% Calculate instant.GR

instantGR = NaN(L,1);

for t = 1:L

    if ( (test_array(t,1) == 1) && (test_array(t,2) == 1) )
            
        t1 = t - floor(w/2);
        t2 = t + floor(w/2);
        cell_area_sub = cell_area(t1:t2);
        
        instantGR(t) = get_instantGR(t1:t2, cell_area_sub); 
        
    end
    
end

% Col[1] time
% Col[2] instantGR

time = (1:L)';
GRt = [time instantGR];

%figure; 
%plot(instantGR, '.-'); hold on;
%plot(cc_flag); hold off; 
%ylim(0.001*[-10 20]);


end

% =====================================================================

function [instantGR] = get_instantGR(taxis, cell_area_sub)

coef = polyfit(taxis, log(cell_area_sub), 1);
instantGR = coef(1);

end

% =====================================================================



