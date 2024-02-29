
dataset_num = 1;
temp = load_optimized_param_V329(dataset_num);
param0 = temp{dataset_num};

% Initial cell size
Amean_data = [2.48 1.51 1.40];
Aini = Amean_data(dataset_num)/(2*log(2));

% =======================================================================

% (0) Create different parameter sets

paramT = {};
rel_vec = [10 1 0.1]';

num = length(rel_vec);

tag_array = { 'r1','r2','d','K1','K2' };
tag_ind = [0 0 1 0 0];  % Choosing the parameter to be varied

for m = 1:num
    
    paramT{m} = param0;
    
    for r = 1:length(tag_ind)
            
        if (tag_ind(r) == 1)
            
            tag = tag_array{r};       
            paramT{m}.(tag) = rel_vec(m)* param0.(tag);

        end
        
    end
    
end

% =======================================================================


% (1A) Perform simulation for various parameters
simuT = {};

for m = 1:num
    simuT{m} = sub_DNA_model_V4(paramT{m});
end
    
% (1B) Interpolate for cell area grid (obtain simuTI)
simuTI = {};

Aq = [0:0.25:15]';
simuTI = simu_interp(simuT, Aq, num); 


% (2A) Calculate the active fractions of RNAP and ribosome

AF = {};

for m = 1:num

    temp = simuTI{m};
    param = paramT{m};
    
    DNA_concWT = param.Zs;      % genome copy / um^3
    DNA_num_1N = 1;             % genome copy / cell

    AF{m}.RNAP.F = ones(length(Aq),1) * ( DNA_concWT / ( DNA_concWT + param.K1 ) );
    AF{m}.RNAP.G = DNA_num_1N ./ ( DNA_num_1N + (param.K1 * param.c * temp.YG) );

    AF{m}.ribo.F = temp.XF ./ ( param.K2 * param.c * temp.YF + temp.XF ) ; 
    AF{m}.ribo.G = temp.XG ./ ( param.K2 * param.c * temp.YG + temp.XG ) ;

end


% (2B) Calculate decay curve for 1N cells

decay_curve = {};
indI = 11;  % initial bin used as "WT size"   

for m = 1:num

    temp = simuTI{m};    
    DNA_conc = DNA_num_1N./temp.VG;
    mRNA_conc = temp.XG./temp.VG;
    
    decay_curve{m}.Anu = temp.AnuG ./ temp.AnuG(indI);
    decay_curve{m}.GR = temp.GRG ./ temp.GRG(indI);

    decay_curve{m}.Zconc =  DNA_conc./DNA_conc(indI);
    decay_curve{m}.Xconc = mRNA_conc ./ mRNA_conc(indI);

    decay_curve{m}.RNAP_AF = AF{m}.RNAP.G ./ AF{m}.RNAP.G(indI);
    decay_curve{m}.ribo_AF = AF{m}.ribo.G ./ AF{m}.ribo.G(indI);

end

% =====================================================================

XLM = 10;

color_vec = [0 0.75 1; 0 0.5 1; 0 0 1];

figure('position', [1 1 600 150]);

subplot(141);

for m = 1:num
    
    % Plot DNA decay
    plot(Aq(indI:end), decay_curve{m}.Zconc(indI:end), 'k-', 'LineWidth',1.5); hold on;
    
    % Plot GR decay
    plot(Aq(indI:end), decay_curve{m}.GR(indI:end), '-','LineWidth',1.5, 'color', color_vec(m,:) ); hold on;
    
end

xlim([2 XLM]);
ylim([0 1]);
 
subplot(142);

for m = 1:num
    
    % Plot DNA decay
    plot(Aq(indI:end), decay_curve{m}.Zconc(indI:end), 'k-', 'LineWidth',1.5); hold on;
    
    % Plot [mRNA] decay
    plot(Aq(indI:end), decay_curve{m}.Xconc(indI:end), '-','LineWidth',1.5, 'color', color_vec(m,:) ); hold on;
    
end

xlim([2 XLM]);
ylim([0 1]);
 

subplot(143);

for m = 1:num
    
    % Plot DNA decay
    plot(Aq(indI:end), decay_curve{m}.Zconc(indI:end), 'k-', 'LineWidth',1.5); hold on;
    
    % Plot RNAP AF decay
    plot(Aq(indI:end), decay_curve{m}.RNAP_AF(indI:end), '-','LineWidth',1.5, 'color', color_vec(m,:) ); hold on;
    
end

xlim([2 XLM]);
ylim([0 1]);


subplot(144);

for m = 1:num
    
    % Plot DNA decay
    plot(Aq(indI:end), decay_curve{m}.Zconc(indI:end), 'k-.', 'LineWidth',1.5); hold on;
    
    % Plot ribo AF decay
    plot(Aq(indI:end), decay_curve{m}.ribo_AF(indI:end), '-','LineWidth',1.5, 'color', color_vec(m,:) ); hold on;
    
end

xlim([2 XLM]);
ylim([0 1]);


% ==========================================================

function [dataI] = simu_interp(simuT, Aq, num); 

    dataI = {};

    for m = 1:num
    
        dataI{m} = {};  % interpolated data

        dataI{m}.Aq = Aq;
                
        dataI{m}.AnuF = interp1(simuT{m}.F.A2, simuT{m}.F.Anu, Aq);
        dataI{m}.AnuG = interp1(simuT{m}.G.A2, simuT{m}.G.Anu, Aq);
        dataI{m}.GRF = dataI{m}.AnuF./dataI{m}.Aq;
        dataI{m}.GRG = dataI{m}.AnuG./dataI{m}.Aq;
        
        dataI{m}.XF = interp1(simuT{m}.F.A, simuT{m}.F.x, Aq);
        dataI{m}.XG = interp1(simuT{m}.G.A, simuT{m}.G.x, Aq);
        
        dataI{m}.YF = interp1(simuT{m}.F.A, simuT{m}.F.y, Aq);
        dataI{m}.YG = interp1(simuT{m}.G.A, simuT{m}.G.y, Aq);
        
        dataI{m}.VF = interp1(simuT{m}.F.A, simuT{m}.F.V, Aq);
        dataI{m}.VG = interp1(simuT{m}.G.A, simuT{m}.G.V, Aq);
        
        %

        dataI{m}.XcF = dataI{m}.XF./dataI{m}.VF;
        dataI{m}.XcG = dataI{m}.XG./dataI{m}.VG;
        
        dataI{m}.AnuR = dataI{m}.AnuG./dataI{m}.AnuF;
        dataI{m}.XR = dataI{m}.XG./dataI{m}.XF;
        dataI{m}.XcR = dataI{m}.XcG./dataI{m}.XcF;

    end

end

% ==========================================================



