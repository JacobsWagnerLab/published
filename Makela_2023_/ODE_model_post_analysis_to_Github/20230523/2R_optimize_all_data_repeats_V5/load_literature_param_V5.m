
function [param] = load_literature_param_V5()

param = {};  % (1) M9GlyCAA (2) M9Gly (3) M9Ala

r1 = [1.99 0.89 0.60]* 10^(-3);   % 1/min
r2 = [20.0 8.21 3.70]* 10^(-3);   % 1/min
K1 = [1.40 2.99 3.03];      % 1/um^3
K2 = [622 1807 1467];       % 1/um^3
d  = [0.964 0.964 0.964];   % 1/min
Zs = [1.40 1.99 2.02];      % genome_copy/ um^3
c = [0.283 0.339 0.530]* 10^(-6);    % um^3 per number of protein 
cp = [0.469 0.616 1.061]* 10^(-6);   % um^2 per number of protein 
xini = [2.02 0.721 0.317]* 10^3;     % number of molecules per cell
yini = [3.82 1.77 0.95]* 10^6;       % number of molecules per cell
Vmean = [1.5 0.83 0.70];    % um^3 

% M9GlyCAA

k1 = 7.5* 10^4;   % 1 /(sec.M)
km1 = 1.1;        % 1 /sec
k2 = 0.012;       % 1 /sec
k3 = 0.0083;      % 1 /sec

for m = 1:3
    
    param{m}.r1 = r1(m);
    param{m}.r2 = r2(m);
    param{m}.K1 = K1(m);
    param{m}.K2 = K2(m);
    param{m}.d = d(m);
    param{m}.Zs = Zs(m);
    param{m}.c = c(m);
    param{m}.cp = cp(m);
    param{m}.xini = xini(m);
    param{m}.yini = yini(m);
    param{m}.Vmean = Vmean(m);    
        
end

param{1}.k1 = k1;
param{1}.km1 = km1;
param{1}.k2 = k2;
param{1}.k3 = k3;


end
