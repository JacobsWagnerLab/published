
function [param] = load_optimized_param_V4()

param = {};  % (1) M9GlyCAA (2) M9Gly (3) M9Ala

r1 = [1.22 0.64 0.46]* 10^(-3);
r2 = [22.4 15.3 13.6]* 10^(-3);
K1 = [1.17 3.02 2.78];      % 1/um^3
K2 = [797 1836 1513];       % 1/um^3
d  = [1.27 1.02 0.79];      % 1/min
Zs = [1.40 1.99 2.02];      % genome_copy/ um^3
c =  [0.283 0.15 0.3]* 10^(-6);
cp = [0.469 0.616 1.061]* 10^(-6);
xini = [2.02 0.721 0.317]* 10^3;
yini = [3.82 1.77 0.95]* 10^6;
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
