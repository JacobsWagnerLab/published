
function [simudata] = sub_simulate_single_param(param_temp, dataset_flag)

%============================================================
% Perform simulation
%simu = sub_DNA_model_V4(param_temp);
simu = sub_DNA_model_PQv2(param_temp);

%============================================================
% Compare input initial condition with steady-state sol

rec = NaN(1,6);  % [1] phiP [2] phi_SS [3] tau(exp) [4] tau(simu)
cell_tau = [40 100 175];

rec(1) = simu.F.phiP;
rec(2) = simu.F.phi;    
rec(3) = 1/rec(1);
rec(4) = 1/rec(2);
%rec(5) = cell_tau(dataset_flag);
%rec(6) = simu.F.tau;

% -----------------------------------------------------------

simudata = {};
simudata.traj = simu;
simudata.rec = rec;

end

%============================================================

