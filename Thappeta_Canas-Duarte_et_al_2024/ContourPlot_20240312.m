%% Imports that data from appropiate analysis to generate countor plots

%for WT data
y_wt_m2r3 = [M2R3_WT_Cells.Norm_ConstrictionSiteDistance];
x_wt_m2r3 = [M2R3_WT_Cells.Norm_NucleoidCentroidDistance];

%for glg data
y_glg_m2r3 = [M2R3_glg_Cells.Norm_ConstrictionSiteDistance];
x_glg_m2r3 = [M2R3_glg_Cells.Norm_NucleoidCentroidDistance];

% Define grid
m2r3_xgrid=linspace(-0.15,0.15); %set the grid depending the range of values for the variables to be plotted
m2r3_ygrid=linspace(-0.15,0.15);
[m2r3_x1,m2r3_y1] = meshgrid(m2r3_xgrid, m2r3_ygrid);

%% Plotting 

m2r3_xi = [m2r3_x1(:) m2r3_y1(:)];
[m2r3_f_wt,m1r3_ep_wt]=ksdensity([x_wt_m2r3' y_wt_m2r3'],m2r3_xi); % remove the outputs to see a 3D plot of the distribution
[m2r3_f_glg,m1r3_ep_glg]=ksdensity([x_glg_m2r3' y_glg_m2r3'],m2r3_xi); % remove the outputs to see a 3D plot of the distribution

figure(8)
hold on
% format data in matrix for contourf and plot
m2r3_X_wt = reshape(m1r3_ep_wt(:,1),length(m2r3_xgrid),length(m2r3_ygrid));
m2r3_Y_wt = reshape(m1r3_ep_wt(:,2),length(m2r3_xgrid),length(m2r3_ygrid));
m2r3_Z_wt = reshape(m2r3_f_wt,length(m2r3_xgrid),length(m2r3_ygrid));
contour(m2r3_X_wt,m2r3_Y_wt,m2r3_Z_wt,10,'-k') %contour for not filled

m2r3_X_glg = reshape(m1r3_ep_glg(:,1),length(m2r3_xgrid),length(m2r3_ygrid));
m2r3_Y_glg = reshape(m1r3_ep_glg(:,2),length(m2r3_xgrid),length(m2r3_ygrid));
m2r3_Z_glg = reshape(m2r3_f_glg,length(m2r3_xgrid),length(m2r3_ygrid));
contour(m2r3_X_glg,m2r3_Y_glg,m2r3_Z_glg,10,'-m') %contour for not filled
ax=gca;
ax.XAxis.Exponent = 0;
xtickformat('%.4f')
xlabel('col1') %x-axis label
ylabel('col2') %y-axis label
hold off
