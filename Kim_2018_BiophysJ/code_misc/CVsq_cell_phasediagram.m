%% This script generates colormap of CVsq_cell-CVsq_gene as a function of CVsq_gene and f

clear;
x = 0:0.01:1; %f
y = 0:0.01:1.5; %CVsq_gene
for xi = 1:length(x)
    for yi = 1:length(y)
        %CVsq_cell-CVsq_DNA
        z(xi,yi)= -(x(xi)*y(yi)/(1+x(xi)))+ (x(xi)*(1-x(xi))/(1+x(xi))^2);        
    end;
end;
z = z';
[C,h] = contour(x,y,z,20);
clabel(C,h)
figure,
h = pcolor(x,y,z), shading interp; %grid off;
set(gca,'YDir','normal') %to flip y axis
Polarmap, colorbar, 
set(gca, 'Layer','top'); %this helps to show ticks and boundary box
xlabel('f'); ylabel('CV^2_D_N_A');



