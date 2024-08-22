
function [data_export] = sub_calculate_error(exp_data, AF_exp_data, simudata, flag, err_weight, err_range)

% Compare experiment and simulation
% Manaully load "GR_exp_data.mat" and "simu_data.mat"

% (1) Calculate the error between measured GR vs model GR
ind1 = (flag*2)-1;
ind2 = flag*2;

nu = {};
nu{1} = get_interp_data(exp_data{ind1}.bindata, simudata.traj.F);
nu{2} = get_interp_data(exp_data{ind2}.bindata, simudata.traj.G);

% Penalize the NaN bins in simulation by setting it to (-10)
nu{1}.S2 = panelize_nan(nu{1}.S2, nu{1}.E, -10);
nu{2}.S2 = panelize_nan(nu{2}.S2, nu{2}.E, -10);

GR_err_vec = {};
GR_err_vec{1} = abs(nu{1}.S2 - nu{1}.E)./nu{1}.E;
GR_err_vec{2} = abs(nu{2}.S2 - nu{2}.E)./nu{2}.E;

GR_err = {};
GR_err{1} = nanmean( GR_err_vec{1}(err_range(1,1):err_range(1,2)) );
GR_err{2} = nanmean( GR_err_vec{2}(err_range(1,1):err_range(1,2)) );

outputGR = {};
outputGR.nu = nu;
outputGR.err = GR_err;

% Load simulation data and interpolate at binned area
[AF_simu] = get_AF_bindata(simudata);

% Load experimental data from SMT
AF_data = {};

if (flag == 1)
    
    AF_data{1} = AF_exp_data{flag}.ftsZ;

elseif (flag == 2)
    
    AF_data{1} = AF_exp_data{flag}.WT;

elseif (flag == 3)
    
    AF_data{1} = AF_exp_data{flag}.WT;

end

AF_data{2} = AF_exp_data{flag}.oriC;

% (2) Calculate the error between measured AF_RNAP vs model AF_RNAP

AF_RNAP = {};

AF_RNAP.exp{1} = AF_data{1}.RNAP.stat(:,2);
AF_RNAP.exp{2} = AF_data{2}.RNAP.stat(:,2);
AF_RNAP.simu{1} = AF_simu.F.bin_AF_RNAP;
AF_RNAP.simu{2} = AF_simu.G.bin_AF_RNAP;

% Penalize the NaN bins in simulation by setting it to (-1)
AF_RNAP.simu{1} = panelize_nan(AF_RNAP.simu{1}, AF_RNAP.exp{1}, -1);
AF_RNAP.simu{2} = panelize_nan(AF_RNAP.simu{2}, AF_RNAP.exp{2}, -1);

AF_RNAP_err_vec = {};
AF_RNAP_err_vec{1} = abs(AF_RNAP.exp{1} - AF_RNAP.simu{1})./AF_RNAP.exp{1};
AF_RNAP_err_vec{2} = abs(AF_RNAP.exp{2} - AF_RNAP.simu{2})./AF_RNAP.exp{2};

AF_RNAP_err = {};
AF_RNAP_err{1} = nanmean( AF_RNAP_err_vec{1}(err_range(2,1):err_range(2,2)) );
AF_RNAP_err{2} = nanmean( AF_RNAP_err_vec{2}(err_range(2,1):err_range(2,2)) );

% (3) Calculate the error between measured AF_RNAP vs model AF_RNAP

AF_ribo = {};

AF_ribo.exp{1} = AF_data{1}.ribo.stat(:,2);
AF_ribo.exp{2} = AF_data{2}.ribo.stat(:,2);

AF_ribo.simu{1} = AF_simu.F.bin_AF_ribo;
AF_ribo.simu{2} = AF_simu.G.bin_AF_ribo;

% Penalize the NaN bins in simulation by setting it to (-1)
AF_ribo.simu{1} = panelize_nan(AF_ribo.simu{1}, AF_ribo.exp{1}, -1);
AF_ribo.simu{2} = panelize_nan(AF_ribo.simu{2}, AF_ribo.exp{2}, -1);

AF_ribo_err_vec = {};
AF_ribo_err_vec{1} = abs(AF_ribo.exp{1} - AF_ribo.simu{1})./AF_ribo.exp{1};
AF_ribo_err_vec{2} = abs(AF_ribo.exp{2} - AF_ribo.simu{2})./AF_ribo.exp{2};

AF_ribo_err = {};
AF_ribo_err{1} = nanmean( AF_ribo_err_vec{1}(err_range(3,1):err_range(3,2)) );
AF_ribo_err{2} = nanmean( AF_ribo_err_vec{2}(err_range(3,1):err_range(3,2)) );


% (4) Summarize error from different experimental types

err = {};
err.table = NaN(3,2); 
% Rows: GR, AF_RNAP, AF_ribo
% Cols: FtsZ, oriC

for c = 1:2
    
    err.table(1,c) = GR_err{c};
    err.table(2,c) = AF_RNAP_err{c};
    err.table(3,c) = AF_ribo_err{c};
    
end

err.weight = err_weight;
err.sum = nansum( nansum(err.table .* err.weight) );  % weighted average

% Export analysis

data_export = {};

data_export.err = err;
data_export.AF_RNAP = AF_RNAP;
data_export.AF_ribo = AF_ribo;
data_export.nu = nu;

% ---------------------------------------------------------------------

%{
figure;

subplot(131);
plot(nu{1}.E, 'bo-'); hold on;
plot(nu{1}.S2, 'b.-'); hold on;
plot(nu{2}.E, 'ro-'); hold on;
plot(nu{2}.S2, 'r.-'); hold on;

subplot(132);
plot(AF_RNAP.exp{1}, 'bo-'); hold on;
plot(AF_RNAP.simu{1}, 'b.-'); hold on;
plot(AF_RNAP.exp{2}, 'ro-'); hold on;
plot(AF_RNAP.simu{2}, 'r.-'); hold on;

xlim([0 15]);
ylim([0 1]);

subplot(133);
plot(AF_ribo.exp{1}, 'bo-'); hold on;
plot(AF_ribo.simu{1}, 'b.-'); hold on;
plot(AF_ribo.exp{2}, 'ro-'); hold on;
plot(AF_ribo.simu{2}, 'r.-'); hold on;

xlim([0 15]);
ylim([0 1]);
%}
%{
figure;
plot(nu{1}.axis, nu{1}.E, 'bo'); hold on;
%plot(nu{1}.axis, nu{1}.S, 'b.-'); hold on;
plot(nu{1}.axis, nu{1}.S2, 'b.-'); hold on;
plot(nu{2}.axis, nu{2}.E, 'ro'); hold on;
%plot(nu{2}.axis, nu{2}.S, 'r.-'); hold on;
plot(nu{2}.axis, nu{2}.S2, 'r.-'); hold on;
%}

end

% =======================================================================

function [vec] = panelize_nan(vec, ref, assigned_val)

for j = 1:length(vec)
    
    if ( isnan(vec(j)) )
        
        if ( ~isnan(ref(j)) )
    
            vec(j) = assigned_val;
        
        end
        
    end
    
end

end    
    
% =======================================================================
% Interpolate the simulated data to compare with the experimental binning data

function [area_nu] = get_interp_data(exp_temp, simu_temp)
    
    area_nu = {};
    area_nu.axis = exp_temp(:,2);
    area_nu.E = exp_temp(:,3);
    area_nu.S = interp1(simu_temp.A2, simu_temp.Anu, area_nu.axis);
    
    % Remove bin with NaN and extrapolate
    
    temp = [area_nu.axis area_nu.S];
    temp( isnan(area_nu.S),: ) = [];    
    %area_nu.S2 = interp1(temp(:,1), temp(:,2), area_nu.axis, 'linear', 'extrap');
    area_nu.S2 = interp1(temp(:,1), temp(:,2), area_nu.axis);
    
end

% =======================================================================

function [output] = get_AF_bindata(dataS)

output = {};

F = dataS.traj.F;
G = dataS.traj.G;
param = dataS.traj.param;

DNA_F = param.Zs;
DNA_G = 1;

output.F.area = F.A;

output.F.act_RNAP = ones(length(F.t),1) * ( DNA_F /(DNA_F + param.K1 * param.Vmean) );

output.F.act_ribo = F.x ./(param.K2 * param.c * F.y + F.x);

output.G.area = G.A;

output.G.act_RNAP = DNA_G ./ (param.K1 * param.c * G.y + DNA_G);

output.G.act_ribo = G.x ./(param.K2 * param.c * G.y + G.x);

area_axis = 0.5*(1:24)';  % um^2
%output.F.bin_AF_RNAP = interp1(output.F.area, output.F.act_RNAP, area_axis, 'linear', 'extrap');
%output.G.bin_AF_RNAP = interp1(output.G.area, output.G.act_RNAP, area_axis, 'linear', 'extrap');
%output.F.bin_AF_ribo = interp1(output.F.area, output.F.act_ribo, area_axis, 'linear', 'extrap');
%output.G.bin_AF_ribo = interp1(output.G.area, output.G.act_ribo, area_axis, 'linear', 'extrap');

output.F.bin_AF_RNAP = interp1(output.F.area, output.F.act_RNAP, area_axis);
output.G.bin_AF_RNAP = interp1(output.G.area, output.G.act_RNAP, area_axis);
output.F.bin_AF_ribo = interp1(output.F.area, output.F.act_ribo, area_axis);
output.G.bin_AF_ribo = interp1(output.G.area, output.G.act_ribo, area_axis);
output.area_axis = area_axis;

end
