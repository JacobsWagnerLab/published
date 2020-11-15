% to make individual traces plots from MultiOutput structures
% parB + ParA distribution
% assumes that dimension of data might be varied

sset=1;
 ttrace=1;
 dt_step=1; % step for plooting, frames
                                      
% plotting parameters
%cols=['r';'g';'b'];
n_bins1=50; % number of bins for ParA distribution
fig_pos= [300,300,1000,300];
cols=[1,0,0; 0.7,0,0; 0,0,1; 0,0,0.7; 0.7,0.7,1];
 shades=0.5;
 col_A='r'; col_B='g'; col0='k';
mrk_A='x'; mrk_B='none';
 ls_A='none'; ls_B='-'; 
ms1=12; ms2=8; ms3=4;
 lw1=3; lw2=2; lw0=1;
fs1=18; fs2=16; fs3=14;
%  tit_add=['; trace#', num2str(ttrace), ' totA=', num2str(totA)];
 x_lab='Time, s'; y_lab='Position, \mum'; 
lgnd={'ParA','PC'};

% GET DATA
Data=DataX{sset}; %DataX4_sigma{sset}; 
 ParB=Data.ParB; 
 ParA_state=Data.ParA_state; 
 ParA=Data.ParA; 
 ParA0=Data.ParA_0;
params=Data.params; 
totA=params.totA;
 totB=params.totB;
dt=params.dt_out;
 t_fin=sum(params.t_fin);
l0=params.l0; w0=params.w0;
sigma_x=params.sigma_x;
if isfield (params, 'dim') 
  dim=params.dim; 
else
  dim=2;
end

tit_add=['; trace#', num2str(ttrace), ' sigma_x=', num2str(sigma_x)];

% prepare additional data 
TT=(0:dt_step*dt:t_fin)';
xmax=w0/2; ymax=l0/2;

% Get single trace data
ParB_1=ParB(0+1:dt_step:t_fin/dt+1,dim*(ttrace-1)*totB+1:dim:dim*(ttrace-1)*totB+dim*totB);
 ParA_1=ParA(0+1:dt_step:t_fin/dt+1,dim*(ttrace-1)*totA+1:dim:dim*(ttrace-1)*totA+dim*totA);
 ParA0_1=ParA0(0+1:dt_step:t_fin/dt+1,dim*(ttrace-1)*totA+1:dim:dim*(ttrace-1)*totA+dim*totA);


% #### PLOTTING #### 
% 1D trace(X-axis) with ParA density as kymograph
figure('Position',fig_pos); hold on
  %tit='Plowing replisome and HarP density'; %
  tit='Simulations example: ParB + ParA density';
 %col=cols(ttrace);
ind=ParB_1(:,1)<111; x_lim=[min(TT(ind)), max(TT(ind))+1];
% cols_X=winter(length(ParB_1(ind)));
% rnd_shift=dt*(rand(sum(ind),totA)-0.5);

tt0=min(TT(ind)); tt1=max(TT(ind));
XX=repmat(TT(ind),1,totA); dx=tt0-dt_step*dt/2:dt_step*dt:tt1+dt_step*dt/2;
 YY=ParA_1(ind,1:end); dy=-ymax:l0/n_bins1:ymax;
 edges{1}=dx; edges{2}=dy;
XX_YY = hist3([XX(:),YY(:)],'Edges',edges);
 XX_YY=XX_YY'; %/median(XX_YY(:));
 disp(sum(XX_YY(:)*(dx(2)-dx(1))*(dx(2)-dx(1)))); 
%imagesc(XX_YY'); ax=gca; colormap(ax,hot(256)); colorbar('location','EastOutside') ;
sc_lim=[min(XX_YY(:)),max(XX_YY(:))]; sc_lim=[sc_lim(1)+0.*diff(sc_lim),sc_lim(2)-0.*diff(sc_lim)];
imagesc(XX_YY,sc_lim); ax=gca; colormap(ax,hot(256)); colorbar('location','EastOutside') ;
%imagesc(XX_YY,sc_lim); ax=gca; colormap(ax,gray(256)); colorbar('location','EastOutside') ;
 plot(repmat(TT(ind)/(dt_step*dt)+1,1,totB),n_bins1*(ParB_1(ind,1:end)+l0/2)/l0,'LineStyle',ls_B,'LineWidth',lw1,'Marker',mrk_B,'MarkerSize',ms3,'Color',col_B,'MarkerFaceColor',col_B,'MarkerEdgeColor',col_B);
 
title([tit,tit_add]);
 set(gca,'FontSize',fs2); xlabel(x_lab,'fontsize',fs1); ylabel(y_lab,'fontsize',fs1);
%xlim(x_lim/dt_step); ylim([0,n_bins1])
xlim(x_lim/(dt*dt_step)); ylim([0,n_bins1])
  set(gca,'XTick',(1:t_fin/5/(dt_step*dt):t_fin/(dt_step*dt)+1)','XTickLabel',(0:t_fin/5:t_fin)','YTick',(0:n_bins1/4:n_bins1)','YTickLabel',(-l0/2:l0/4:l0/2)'); 
 
  
% 1D trace(X-axis) with equlibrium positions of ParA as a kymograph
figure('Position',fig_pos); hold on
 tit='Simulations example: ParB + equlibrium ParA positions';
 %col=cols(ttrace);
ind=ParB_1(:,1)<111; x_lim=[min(TT(ind)), max(TT(ind))+1];
% cols_X=winter(length(ParB_1(ind)));
% rnd_shift=dt*(rand(sum(ind),totA)-0.5);

tt0=min(TT(ind)); tt1=max(TT(ind));
XX=repmat(TT(ind),1,totA); dx=tt0-dt_step*dt/2:dt_step*dt:tt1+dt_step*dt/2;
 YY=ParA0_1(ind,1:end); dy=-ymax:l0/n_bins1:ymax;
 edges{1}=dx; edges{2}=dy;
XX_YY = hist3([XX(:),YY(:)],'Edges',edges);
 XX_YY=XX_YY'; %/median(XX_YY(:));
 disp(sum(XX_YY(:)*(dx(2)-dx(1))*(dx(2)-dx(1)))); 
imagesc(XX_YY); ax=gca; colormap(ax,hot(256)); colorbar('location','EastOutside') ;
 plot(repmat(TT(ind)/(dt_step*dt)+1,1,totB)+1,n_bins1*(ParB_1(ind,1:end)+l0/2)/l0,'LineStyle',ls_B,'LineWidth',lw1,'Marker',mrk_B,'MarkerSize',ms3,'Color',col_B,'MarkerFaceColor',col_B,'MarkerEdgeColor',col_B);
 
title([tit,tit_add]);
 set(gca,'FontSize',fs2); xlabel(x_lab,'fontsize',fs1); ylabel(y_lab,'fontsize',fs1);
%xlim(x_lim/dt_step); ylim([0,n_bins1])
xlim(x_lim/(dt*dt_step)); ylim([0,n_bins1])
% set(gca,'XTickLabel',(0:max(TT)/5:max(TT))'); set(gca,'YTickLabel',(-l0/2:l0/5:l0/2)'); 
 set(gca,'XTick',(1:t_fin/5/(dt_step*dt):t_fin/(dt_step*dt)+1)','XTickLabel',(0:t_fin/5:t_fin)','YTick',(0:n_bins1/4:n_bins1)','YTickLabel',(-l0/2:l0/4:l0/2)'); 
  
 
% Attempt to make 3D figure
figure; hold on
   cols3=cat(3,XX_YY.*ones(size(XX_YY))/max(XX_YY(:)),zeros(size(XX_YY)),zeros(size(XX_YY)));
  surf(XX_YY,cols3,'EdgeColor','none')
  ind_c=repmat(TT(ind)/dt+1,1,totB); ind_r=round(n_bins1*(ParB_1(ind,:)+l0/2)/l0+1);
   ind_lin=sub2ind(size(XX_YY),ind_r,ind_c);
  ZZ=XX_YY(ind_lin);
  ls_B='-'; mrk_B='none';
   %plot3(repmat(TT(ind)/dt,1,totB)+1,n_bins1*(ParB_1(ind,1:end)+l0/2)/l0,ZZ,'LineStyle',ls_B,'LineWidth',lw1,'Marker',mrk_B,'MarkerSize',ms3,'Color',col_B,'MarkerFaceColor',col_B,'MarkerEdgeColor',col_B);
   plot3(repmat(TT(ind)/dt+1,1,totB),round(n_bins1*(ParB_1(ind,:)+l0/2)/l0+1),ZZ,'LineStyle',ls_B,'LineWidth',lw1,'Marker',mrk_B,'MarkerSize',ms3,'Color',col_B,'MarkerFaceColor',col_B,'MarkerEdgeColor',col_B);

