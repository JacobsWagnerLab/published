
% Manaully load ATP_traj_data.mat.
% This code analyzes the correlation decay curve, the data is summarized 
% in "corr_data_ensemble". For the decay plot, use another file "Plot_corr_decay.m".

path.save = '/Users/wei-hsiang/Desktop/uF_local/Tree_data/sibling_decorrelation/';

dataset_name = {'M9GlcCA_lac', 'M9GlcCA_ldhA', 'M9GlcCA_adhE', 'M9GlcCA_pta', 'M9GlcCA_ackA', ...
                'M9Glc_lac', 'M9Xyl_lac', 'M9Gly_lac'};

% ======================= Script starts ================================

mtype_array = [4 4 4 4 4 4 2 2]; % Specify the maximal mtype to use for each conditions 

corr_data_ensemble = {};

for DB = 1:length(dataset_name)
    
    for mtype = 1:mtype_array(DB)
        
        paramA = {};
        paramA.DB = DB;
        paramA.mtype = mtype;
        paramA.Tmax = 60;             % maximal time frame
        paramA.frame_interval = 6; 
        paramA.cutoff = 20;
        
        % Calculate decorrelation curve
        [output] = sub_decorrelation_ensemble(ATP_traj_data, paramA, path.save, dataset_name);       
        corr_data_ensemble{DB}{mtype} = output;

    end
        
end

%  Filter the decorelation curve by fraction of cell cycle

cc_fraction = 0.25;

for DB = 1:8
    
    for mtype = 1:mtype_array(DB)
    
        data = corr_data_ensemble{DB}{mtype}.corr_decay;
        count = corr_data_ensemble{DB}{mtype}.sibpair_count;
    
        flag = 1 + count*0;

        for j = 1:length(count)
        
            if ( count(j)/count(1) < cc_fraction )  % Too few cell cycle remains            
                flag(j) = NaN;            
            end
        
        end
    
        corr_data_ensemble{DB}{mtype}.corr_decayF = data.*flag;    
    
    end
    
end

% ========================================================================

function [output] = sub_decorrelation_ensemble(ATP_traj_data, paramA, save_dir, dataset_name)

% Decorrelation using trajectory ensemble for each time points

DB = paramA.DB;
mtype = paramA.mtype;
Tmax = paramA.Tmax;   % maximal time frame
frame_interval = paramA.frame_interval; 
cutoff = paramA.cutoff;

% ---------------------------------------------------------------------

sib_data_total = {};
  
data0 = ATP_traj_data{DB}.data.ATP_traj;
kmax = length(data0{mtype});

for k = 1:kmax

    data_traj = data0{mtype}{k};    
    sib_data_total{k} = get_sib_traj(data_traj, frame_interval);
    
end

corr_data = {};

for t = 1:Tmax
    
    corr_data{t} = [];
    
    for k = 1:kmax
        
        if ( sib_data_total{k}.L >= t )
            
            write = [sib_data_total{k}.c1(t) sib_data_total{k}.c2(t)];
            corr_data{t} = [corr_data{t}; write];            
            
        end
        
    end
    
end


% ====================================================================
% Calculating correlation decay 

corr_decay = NaN(Tmax, 1);
sibpair_count = NaN(Tmax, 1);

for t = 1:Tmax
    
    sibpair_count(t) = size(corr_data{t}, 1);
    
    if ( size(corr_data{t}, 1) > cutoff )        
        corr_decay(t,1) = corr(corr_data{t}(:,1), corr_data{t}(:,2), 'type', 'Pearson');
        corr_decay(t,2) = corr(corr_data{t}(:,1), corr_data{t}(:,2), 'type', 'Kendall');
        corr_decay(t,3) = corr(corr_data{t}(:,1), corr_data{t}(:,2), 'type', 'Spearman');
    end
    
end

% Save correlation decay data

output = [];
output.corr_data = corr_data;
output.corr_decay = corr_decay;
output.sibpair_count = sibpair_count;
output.paramA = paramA;

end


% ====================================================================

function [sib_data] = get_sib_traj(data_traj, frame_interval)

    sib_data = {};
    sib_data.c1 = data_traj.cc1;
    sib_data.c2 = data_traj.cc2;
    sib_data.L = min(length(sib_data.c1), length(sib_data.c2) );
    sib_data.time = frame_interval *(1:sib_data.L);

end

% ====================================================================


