
function [param] = load_literature_param()

param = {};  % (1) M9GlyCAA (2) M9Gly (3) M9Ala

r1 = [0.568 0.616 0.548]* 10^(-3);
r2 = [20.0 8.21 3.70]* 10^(-3);
K1 = [1.0 1.8 2.14];        % 1/um^3
K2 = [622 1807 1467];       % 1/um^3
d  = [0.251 0.238 0.221];   % 1/min
c = [0.284 0.339 0.530]* 10^(-6);
cp = [0.537 0.616 1.061]* 10^(-6);
xini = [20.20 7.21 3.17]* 10^2;
yini = [3.82 1.77 0.95]* 10^6;
Zs = [1.5 1 1];
Vmean = [1.5 0.83 0.70];    % um^3 

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

end
