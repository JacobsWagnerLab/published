function [Stats,Stats2,out]=plot_MultiplePlots_DataX_2_xD(DataX, set_to_plot, pplots, tit_0, cols, value_2_show)
% changed to accomodate data from simulations in different dimensions
% check all plots!!!
% to make various plots from simulated data by applying plot_MultiplePlots to different datasets in MultipleData :
%
% plot_MultiplePlots(MultiOutput, sset, pplots)
% - "snapshots" at given moment of time for multiple traces
% - etc

%t=2000;

Stats=[];Stats2=[];

t0=600; %500; % transition period
 t2=30; % stop-time t=t2+t0
dt=1;
%sset=1;
n_bins=100; n_bins_2=200; n_bins_3=5000;
n_sm=6000; % for smoothing (for example for the force plot)
 tau_sm=20; % for smoothing with Gaussian kernel, s
n_tr=2;
sim=6; % what simulation to show in trace examples
r0=0.2; % cutoff for steps correlation fitting
ind_max=900;
nn_max=5; % max iterations when step is changed, i.e. for displacement with variable time step n*dt

% plotting params
if nargin<6, value_2_show='totA'; end
if nargin <5, cols=['b';'g';'r';'c';'m';'k']; cols2=prism(n_tr);   end
if nargin<4, tit_0=' '; end%['Varied: ',MultiOutput.varied_params];


 ms=4; ms2=4; 
 lw1=3; lw2=2; lw3=1;
 fs1=18; fs2=16; fs3=11;
%colormap=cols;
 
 
disp(num2str(['n_runs','    dt_out','  t_fin','       totB','   totA','    khyd','     k_db' ],4)) 

out={};
lgnd_1={}; lgnd_4={};
for jj=1:length(set_to_plot)
sset=set_to_plot(jj);  
  
%txt_01 = num2str(params.(fldnms{kk}),3);
    
% getting data
Data=DataX{sset};
params=Data.params;
 n_runs=Data.n_runs;
 dt_out=params.dt_out;
  t_fin=sum(params.t_fin);
   t1=t_fin;   % cutoff time, i.e. to measure PCs distribution
 totB=params.totB;
  totA=params.totA;
 k_hyd=params.k_hyd;
 if isfield(params, 'k_db')
   k_db=params.k_db;
 elseif isfield(params, 'k_db0')
  k_db=[params.k_db0,params.k_db11,params.k_db12];
 else
   k_db=0;  
 end
 if isfield (params, 'dim') 
  dim=params.dim; 
else
  dim=2;
end
 sigma_x=params.sigma_x;
  sigma_y=params.sigma_y;
 k_sp=1/(sigma_x^2);  
 D_parB=params.D_parB;
 l0=params.l0;
 

ParB=Data.ParB;
if isfield (Data,'ParA')
  ParA=Data.ParA; 
   ParA_state=Data.ParA_state; 
   ParA0=Data.ParA_0;
end


% ParA_1=ParA(:,2*totA*(ttrace-1)+1:2*ttrace*totA);
%   ParA_1x=ParA_1(:,1:2:end);
%  ParA0_1=ParA0(:,2*totA*(ttrace-1)+1:2*ttrace*totA);
%   ParA0_1x=ParA0_1(:,1:2:end);
%  ParA_state_1=ParA_state(:,totA*(ttrace-1)+1:ttrace*totA);


lgnd_1{jj}=num2str(params.(value_2_show));

%  %  Some Statitistics
% X=ParB(t1/dt_out+1,1:2:end);
%  Xcell=reshape(X,totB,n_runs);
%  Xcell=Xcell';
%  Xcell=sort(Xcell,2); 
% no_left=find(Xcell(:,1)>0); no_right=find(Xcell(:,totB)<0); no_PC=length(no_left)+length(no_right);
%  only1_left=Xcell(:,2)>0; only1_right=Xcell(:,totB-1)<0; only1=sum(only1_left)+sum(only1_right); 
%  left=sum((Xcell<0),2); lr_ratio=left/totB; 
%  Stat.Fraction_wo_PC=no_PC/n_runs;
%   Stat.only1=only1/n_runs;
%  Stat.lr_ratio=[mean(lr_ratio),std(lr_ratio),min(lr_ratio),max(lr_ratio)];
% Stats(jj,:)=[no_PC, Stat.Fraction_wo_PC,only1,Stat.lr_ratio];
% 
% disp(num2str([n_runs, dt_out, t_fin, totB, totA, k_hyd, k_db, sigma_x, sigma_y],4))
% 
% %same thing but "integrated" between t0 and t1
% X_all=[];
% for tt=t0:dt_out:t1
%   X=ParB(int32(tt/dt_out+1),1:2:end);
%    Xcell2=reshape(X,totB,n_runs);
%    Xcell2=Xcell2';
%    Xcell2=sort(Xcell2,2);
%    X_all=[X_all;Xcell2];
% end
% no_left=find(X_all(:,1)>0); no_right=find(X_all(:,totB)<0); no_PC=length(no_left)+length(no_right);
%  only1_left=X_all(:,2)>0; only1_right=X_all(:,totB-1)<0; only1=sum(only1_left)+sum(only1_right); 
%  left=sum((X_all<0),2); lr_ratio=left/totB; 
%  Stat.Fraction_wo_PC_int=no_PC/n_runs/((t1-t0)/dt);
%   Stat.Fraction_only1_int=only1/n_runs/((t1-t0)/dt);
%  Stat.lr_ratio_int=[mean(lr_ratio),std(lr_ratio),min(lr_ratio),max(lr_ratio)];
% Stats2(jj,:)=[no_PC, Stat.Fraction_wo_PC_int,Stat.Fraction_only1_int,Stat.lr_ratio_int];


for kk=1:length(pplots)
switch pplots(kk)
  case 1 % FIG.1     % plotting all PC positions at t=t;
%     if exist('fig1','var'), figure(fig1); else fig1=figure; hold on; end
    figure; hold on
    ylab='Distance to pole, \mum';
     xlab='Trace #';
     tit_1=['PC positions at t' ' = ',num2str(t_fin), ' ',value_2_show,' = ', num2str(params.(value_2_show))];
     % auxilary params
    x0=0; %params.l0/2+0.5; % "pole" coordinate
    nB=totB;
     % sorting PCs by distance from the pole
    X=ParB(t1/dt_out+1,1:dim:end);
     Xcell=reshape(X,nB,n_runs);
     Xcell=Xcell';
     Xcell=sort(Xcell,2); 
     %plotting
       %plot(1:n_runs, x0-Xcell,'o','MarkerEdgeColor',cols,'MarkerFaceColor',cols,'MarkerSize',ms)
    plot(1:n_runs,x0*ones(1,n_runs), ':','Color',[0.5,0.5,0.5],'LineWidth',lw2)
     %plot(1:n_runs,x0-Xcell, 'o','MarkerSize',ms,'LineWidth',lw1)
     plot(1:n_runs,x0-Xcell, 'o','MarkerSize',ms,'LineWidth',lw1)
    ylim([0,params.l0+1]-0.5*(params.l0+1)); 
    xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
     ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
    set(gca,'FontSize', fs2)
    title(tit_1,'FontSize', fs3)
      
%      %  Some Statitistics
%     no_left=find(Xcell(:,1)>0); no_right=find(Xcell(:,totB)<0); no_PC=length(no_left)+length(no_right);
%      left=sum((Xcell<0),2); lr_ratio=left/totB; 
%     Stat.Fraction_wo_PC=no_PC/n_runs/2;
%      Stat.lr_ratio=[mean(lr_ratio),std(lr_ratio),min(lr_ratio),max(lr_ratio)];
% %     disp(['cells w/o PC',' Fraction w/o PC',' Mean    ',' Std   ','  min   ','  max ',' of (#PC in DaugtherCell)' ]);
% %      disp([no_PC, Stat.Fraction_wo_PC,Stat.lr_ratio]);
%     Stat(jj,:)=[no_PC, Stat.Fraction_wo_PC,Stat.lr_ratio];

  case 2 % FIG.2 Hist of relative positions
    if exist('fig2','var'), figure(fig2); else fig2=figure; hold on; end  
      %figure; hold on
     ylab='Frequency';
     xlab='Relative position'; 
     %getting counts
    X2=ParB(t0/dt_out+1:dt/dt_out:t1/dt_out+1,1:2:end);
     [HistY,HistX]=hist(X2(:)./l0,n_bins);
     HistYn=HistY/sum(HistY)/(HistX(2)-HistX(1));
     %plotting
    plot(HistX,HistYn, '-','Color',cols(jj,:),'LineWidth',lw2)
    xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
     ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
    set(gca,'FontSize', fs2)
    legend(lgnd_1, 'FontSize', fs3); title(tit_0,'FontSize', fs3)
  
  case 3 % FIG.3 Hist of abs (relative positions)
    if exist('fig3','var'), figure(fig3); else fig3=figure; hold on; end
        %figure; hold on
     ylab='Frequency';
     xlab='Position, \mum'; 
     %getting counts
    [HistY,HistX]=hist(abs(X2(:)),n_bins);
     HistYn=HistY/sum(HistY)/(HistX(2)-HistX(1));
     %plotting
    plot(HistX,HistYn, '-','Color',cols(jj,:),'LineWidth',lw2)
    xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
     ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
    set(gca,'FontSize', fs2)
    legend(lgnd_1, 'FontSize', fs3); title(tit_0,'FontSize', fs3)

  case 4 % FIG.4 Hist of inter-PC distance
    if exist('fig4','var'), figure(fig4); else fig4=figure; hold on; end
        %figure; hold on
     ylab='Frequency';
     xlab='Distance, \mum'; 
     %getting counts
    X3=diff(Xcell');
     [HistY,HistX]=hist(X3(:),n_bins); 
     HistYn=HistY/sum(HistY)/(HistX(2)-HistX(1));
     %plotting
    plot(HistX,HistYn, 'o-','Color',cols(jj,:),'LineWidth',lw2)
    xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
     ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
    set(gca,'FontSize', fs2)
    legend(lgnd_1, 'FontSize', fs3); title(tit_0,'FontSize', fs3)
  
  case 5 % FIG.5 PC-steps
    if exist('fig5','var'), figure(fig5); else fig5=figure; hold on; end
       %figure; hold on
     ylab='Frequency';
     xlab='Displacement, \mum'; 
     %getting counts
    dX=diff(ParB(1:dt/dt_out:end,:));
     [HistY,HistX]=hist(abs(dX(:)),n_bins_2); 
     HistYn=HistY/sum(HistY)/(HistX(2)-HistX(1));
    dX_gauss=exp(-(HistX.^2/(2*2*dt*D_parB)))/sqrt(pi*dt*D_parB);
     %plotting
    plot(HistX,HistYn, 'o-','Color',cols(jj,:),'LineWidth',lw2)
     plot(HistX,dX_gauss, '-','Color','k','LineWidth',lw2)
    xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
     ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
    set(gca,'FontSize', fs2)
    legend(lgnd_1, 'FontSize', fs3); title(tit_0,'FontSize', fs3)

  case 6 % FIG.6 Traces examples !! changed for varied dim
      %if exist('fig6','var'), figure(fig6); else fig6=figure; hold on; end  
    %figure; hold on
     fig_pos= [300,300,750,500]; figure('Position',fig_pos);
     ylab='Position, \mum';
     xlab='Time, s'; 
      tit_1=['Example: sim#' ' = ',num2str(sim), ' ',value_2_show,' = ', num2str(params.(value_2_show))];
     TT=0:dt:t_fin;
      XX=ParB(1:dt/dt_out:end,dim*totB*(sim-1)+1:dim:dim*totB*sim); 
     out.plot6_x=TT;
      out.plot6_y=XX;
     %plotting
   plot(TT,XX, '-','LineWidth',lw2)
    xlim([0,t_fin]); ylim([-params.l0/2-0.5,params.l0/2+0.5]);
%     plot(t0:dt_out*dt:t2,(ParB(t0/dt_out:dt:t2/dt_out,2*totB*(sim-1)+1:2:2*totB*sim)/l0), 'o','MarkerFaceColor',cols(jj,:),'MarkerEdgeColor',cols(jj,:),'LineWidth',lw2) 
%      xlim([t0,t2]); ylim([-0.5,+0.5]);
    xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
     ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
    set(gca,'FontSize', fs2)
    title(tit_1,'FontSize', fs3)
    
  case 7 % Fig.7: FFT
     if exist('fig7','var'), figure(fig7); else fig7=figure; hold on; end  
     ylab='Amplitude'; xlab='1/\nu, s'; 
     ParBcut=ParB(t0/dt_out+1:t1/dt_out+1,1:2:end);
      ParB_2=ParBcut-repmat(ParBcut(1,:),size(ParBcut,1),1);
      ParB_x=reshape(ParB_2(1:end,1:end),n_runs*((t_fin-t0)/dt_out+1)*totB,1);
     FF=fft(ParB_x,(t_fin-t0)/dt_out);
      w=1/dt_out;
     WW=w*linspace(0,1,(t_fin-t0)/dt_out);
     FF2=mean(abs(FF), 2);
     %semilogx(1./WW(1:floor(end/2)), abs(FF2(1:floor(end/2))), 'LineWidth', lw1, 'Color',col) 
     %loglog(1./WW(1:floor(end/2)), abs(FF2(1:floor(end/2))), 'LineWidth', lw1, 'Color',cols(jj)) 
     plot(1./WW(1:floor(end/2)), smooth(FF2(1:floor(end/2)), 9), 'LineWidth', lw1, 'Color',cols(jj,:)) 
       xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
        ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
       set(gca,'FontSize', fs2)
      legend(lgnd_1, 'FontSize', fs3); title(tit_0,'FontSize', fs3)
  
    case 8 % Fig.8: Next step correlation:  binned, vs abs(dx)
        if exist('fig8','var'), figure(fig8); else fig8=figure; hold on; end
         ylab='Next displacement, \mum'; xlab='Displacement, \mum'; 
        dX=diff(ParB(t0/dt_out+1:dt/dt_out:end,1:2:end));
         ind=1:size(dX,1)-1;
          XX1=abs(dX(ind,:));
           YY1=sign(dX(ind,:)).*dX(ind+1,:);
            disp(['total steps: ',num2str(length(XX1))]);
          [XX,YY]=bin2_fixN(XX1(:),YY1(:), n_bins_3); NN=ones(size(XX1)); tt=['1 bin=',num2str(n_bins_3),'steps '];  
          model='k*x';
           [coeff,fit1,YY1_fit]=oneFit(double(XX(XX<r0)),double(YY(XX<r0)),model,NN(XX<r0)); %abs(1./YY(XX<r0))
             ci = confint(fit1); n_pars=length(coeff);
          lgnd_2{2*(jj-1)+1} =lgnd_1{jj}; lgnd_2{2*jj} =['Lin fit, a=', num2str(coeff,2)];   
         plot(XX,YY, 'o','MarkerFaceColor',cols(jj,:),'MarkerEdgeColor',cols(jj,:),'LineWidth',lw2)
           plot(XX(XX<r0),YY1_fit, '-','Color','k','LineWidth',lw2)
        xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
         ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
        set(gca,'FontSize', fs2)
        legend(lgnd_2, 'FontSize', fs3); title(tit_0,'FontSize', fs3) 
          tit=[tt,' fit by ',model,' for dx<',num2str(r0),' dt=',num2str(dt)];
        title(tit,'FontSize', fs3)  
        
   case 9 % Fig.9: MSD-tau
     if exist('fig9','var'), figure(fig9); else fig9=figure; hold on; end  
     ylab='MSD, \mum'; xlab='Time, s'; 
     ParBcut=ParB(t0/dt_out+1:t1/dt_out+1,1:2:end);
      %ParB_2=ParBcut-repmat(ParBcut(1,:),size(ParBcut,1),1);
     %ind_max=min([size(ParBcut,1)-1,ind_max]);
%       ParB_x=reshape(ParB_2(1:end,1:end),n_runs*((t_fin-t0)/dt_out+1)*totB,1);
     for ii=0:ind_max
        ind1=1:size(ParBcut,1)-ii; ind2=ind1+ii;
        dX=ParBcut(ind2,:)-ParBcut(ind1,:);
        MSD(ii+1)=sum((dX(:)).^2)/sum(ones(size(dX(:))));
     end
    TT=dt_out*(0:ind_max)';
     plot(TT, MSD,'o','MarkerFaceColor',cols(jj,:),'MarkerEdgeColor',cols(jj,:),'LineWidth',lw2) 
       xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
        ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
       set(gca,'FontSize', fs2)
      legend(lgnd_1, 'FontSize', fs3); title(tit_0,'FontSize', fs3)   
      
   case 10 % FIG.10 Individual ACFs examples
      %if exist('fig6','var'), figure(fig6); else fig6=figure; hold on; end  
    figure; hold on
     ylab='ACF';
     xlab='Time, s'; 
      tit_1=['Example: sim#' ' = ',num2str(sim), ' ',value_2_show,' = ', num2str(params.(value_2_show))];
     %plotting
    ParB_cut1=ParB(t0/dt_out+1:t1/dt_out+1,2*totB*(sim-1)+1:2:2*totB*sim);
     cols_2=hsv(size(ParB_cut1,2));
    acf_tot=zeros(size(ParB_cut1(:,1))); 
    plot([0,t_fin-t0],[0,0],'-k')
    for ii=1:size(ParB_cut1,2); 
      [acf,tau,bounds]=autocorr(ParB_cut1(:,ii),size(ParB_cut1,1)-1); 
       acf_tot=acf_tot+acf;
    end 
     plot(tau*dt_out,acf_tot/(totB), '.-','LineWidth',lw2,'Color',cols(jj,:))
     for ii=1:size(ParB_cut1,2);
       [acf,tau,bounds]=autocorr(ParB_cut1(:,ii),size(ParB_cut1,1)-1); 
        acf_tot=acf_tot+acf;
         plot(tau*dt_out,acf,'-','LineWidth',lw2,'Color',cols_2(ii,:))
     end  
    xlim([0,t_fin-t0]); %ylim([-params.l0/2-0.5,params.l0/2+0.5]);
    xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
     ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
    set(gca,'FontSize', fs2)
     legend({'averaged' ,'individual'}, 'FontSize', fs3); title(tit_1,'FontSize', fs3)  
     
   case 11 % Fig.11: Averaged ACF
     if exist('fig11','var'), figure(fig11); else fig11=figure; hold on; end  
     ylab='ACF'; xlab='Time, s'; 
     ParBcut=ParB(t0/dt_out+1:t1/dt_out+1,1:2:end);
      acf_tot=zeros(size(ParBcut(:,1)));
     for ii=1:size(ParBcut,2);
       [acf,tau,bounds]=autocorr(ParBcut(:,ii),size(ParBcut,1)-1); 
        acf_tot=acf_tot+acf;
     end
    TT=dt_out*tau';
     plot(TT, acf_tot/(totB*n_runs),'.-','MarkerFaceColor',cols(jj,:),'MarkerEdgeColor',cols(jj,:),'Color',cols(jj,:),'LineWidth',lw2) 
       xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
        ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
       set(gca,'FontSize', fs2)
      legend(lgnd_1, 'FontSize', fs3); title(tit_0,'FontSize', fs3)   
      
   case 12 % FIG.12 Individual displacement ACFs examples
      %if exist('fig6','var'), figure(fig6); else fig6=figure; hold on; end  
    figure; hold on
     ylab='ACF';
     xlab='Time, s'; 
      tit_1=['Example: sim#' ' = ',num2str(sim), ' ',value_2_show,' = ', num2str(params.(value_2_show))];
     %plotting
    ParB_cut1=diff(ParB(t0/dt_out+1:t1/dt_out+1,2*totB*(sim-1)+1:2:2*totB*sim));
     cols_2=hsv(size(ParB_cut1,2));
    acf_tot=zeros(size(ParB_cut1(:,1))); 
    plot([0,t_fin-t0],[0,0],'-k')
    for ii=1:size(ParB_cut1,2); 
      [acf,tau,bounds]=autocorr(ParB_cut1(:,ii),size(ParB_cut1,1)-1); 
       acf_tot=acf_tot+acf;
    end 
     plot(tau*dt_out,acf_tot/(totB), '.-','LineWidth',lw2,'Color',cols(jj,:))
     for ii=1:size(ParB_cut1,2);
       [acf,tau,bounds]=autocorr(ParB_cut1(:,ii),size(ParB_cut1,1)-1); 
        acf_tot=acf_tot+acf;
         plot(tau*dt_out,acf,'-','LineWidth',lw2,'Color',cols_2(ii,:))
     end  
    xlim([0,t_fin-t0]); %ylim([-params.l0/2-0.5,params.l0/2+0.5]);
    xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
     ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
    set(gca,'FontSize', fs2)
     legend({'averaged' ,'individual'}, 'FontSize', fs3); title(tit_1,'FontSize', fs3)
        title(tit_1,'FontSize', fs3)
     
   case 13 % Fig.13: Averaged displacement ACF
     if exist('fig13','var'), figure(fig13); else fig13=figure; hold on; end  
     ylab='Displacement ACF'; xlab='Time, s'; 
     ParBcut=diff(ParB(t0/dt_out+1:t1/dt_out+1,1:dim:end));
      acf_tot=zeros(size(ParBcut(:,1)));
     for ii=1:size(ParBcut,2);
       [acf,tau,bounds]=autocorr(ParBcut(:,ii),size(ParBcut,1)-1); 
        acf_tot=acf_tot+acf;
     end
    TT=dt_out*tau';
     plot(TT, acf_tot/(totB*n_runs),'.-','MarkerFaceColor',cols(jj,:),'MarkerEdgeColor',cols(jj,:),'Color',cols(jj,:),'LineWidth',lw2) 
       xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
        ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
       set(gca,'FontSize', fs2)
      legend(lgnd_1, 'FontSize', fs3); title(tit_0,'FontSize', fs3)
      
   case 14 % Fig.14: Averaged displacement ACF vs delta_t
     figure; hold on
      ylab='ACF'; xlab='Time, s';
       tit_1=['dx(dt) with varied dt ',value_2_show,' = ', num2str(params.(value_2_show))];
      cols_2=hsv(nn_max);
     for nn=1:nn_max
       %ParBcut=diff_n(ParB(t0/dt_out+1:t/dt_out+1,1:2:end),nn);
        ParBcut=diff(ParB(t0/dt_out+1:nn:t1/dt_out+1,1:2:end));
        acf_tot=zeros(size(ParBcut(:,1)));
       for ii=1:size(ParBcut,2);
         [acf,tau,bounds]=autocorr(ParBcut(:,ii),size(ParBcut,1)-1); 
          acf_tot=acf_tot+acf;
       end
       TT=dt_out*nn'*tau';
       plot(TT, acf_tot/(totB*n_runs),'.-','MarkerFaceColor',cols_2(nn,:),'MarkerEdgeColor',cols_2(nn,:),'Color',cols_2(nn,:),'LineWidth',lw2)
       lgnd_3{nn}=['dt=',num2str(dt_out*nn)];
     end
     xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
      ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
       set(gca,'FontSize', fs2)
       legend(lgnd_3, 'FontSize', fs3); %title(tit_0,'FontSize', fs3)
       title(tit_1,'FontSize', fs3)
       
   case 15 % Fig.15: distance between "replicated" PCs, assumes tha they have ID# totB and totB-1 
     if exist('fig15','var'), figure(fig15); else fig15=figure; hold on; end  
     ylab='MSD, \mum'; xlab='Time, s'; 
     ParB_2=ParB(t0/dt_out+1:t_fin/dt_out+1,2*totB-1:2*totB:end); ParB_1=ParB(t0/dt_out+1:t_fin/dt_out+1,2*totB-3:2*totB:end);
      %dist_21_mean=mean(ParB_2-ParB_1,2);
      dist_21_mean=mean(abs(ParB_2-ParB_1),2);
%       dist_21_msd=mean((ParB_2-ParB_1).^2,2);
     TT=(t0:dt_out:t_fin)';
     plot(TT, dist_21_mean,'o-','MarkerFaceColor',cols(jj,:),'MarkerEdgeColor',cols(jj,:),'LineWidth',lw2) 
%       plot(TT, dist_21_msd,'x:','MarkerFaceColor',cols(jj,:),'MarkerEdgeColor',cols(jj,:),'LineWidth',lw2) 
       xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
        ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
       set(gca,'FontSize', fs2)
      legend(lgnd_1, 'FontSize', fs3); title(tit_0,'FontSize', fs3)      
    
   case 16 % FIG.16 Single Trace: number of PC-Bound ParA vs t
      %if exist('fig6','var'), figure(fig6); else fig6=figure; hold on; end  
    figure; hold on
     ylab='PC-bound ParA';
     xlab='Time, s'; 
      tit_1=['Example: sim#' ' = ',num2str(sim), ' ',value_2_show,' = ', num2str(params.(value_2_show))];    
    n_boundA=[];Force=[];
    %ParB_1=ParB(t0/dt_out+1:t_fin/dt_out+1,2*totB*(sim-1)+1);
    ParA0_1=ParA0(t0/dt_out+1:dt/dt_out:t_fin/dt_out+1,2*totA*(sim-1)+1:2:2*totA*sim);
     ParA_1=ParA(t0/dt_out+1:t_fin/dt_out+1,2*totA*(sim-1)+1:2:2*totA*sim);
     ParA_state_1=ParA_state(t0/dt_out+1:t_fin/dt_out+1,totA*(sim-1)+1:totA*sim);
     ParB_1=ParB(t0/dt_out+1:dt/dt_out:t_fin/dt_out+1,2*totB*(sim-1)+1:2:2*totB*sim);
    for tt=1:size(ParA_1,1)
      for bb=1:totB
        n_boundA(tt,bb)=sum(ParA_state_1(tt,:)==bb);
        ParBx1=ParB_1(tt,bb);
         boundAX=ParA0_1(tt,ParA_state_1(tt,:)==bb);
         n_boundA2(tt,bb)=sum(boundAX>=ParBx1)-sum(boundAX<=ParBx1);
%         indx=find(ParA_state_1(ii,:)>0);
%         Force(ii)=-sum(ParA_1x(ii,indx)-ParA0_1x(ii,indx));
      end
    end
    TT=(t0:dt_out:t_fin)';
    out.plot16_x=TT; out.plot17_y=n_boundA;
     %plotting
    plot(TT,n_boundA, '-','LineWidth',lw3)
     plot(TT,n_boundA2, '-','LineWidth',lw3,'Color',[0.,0.5,0.5])
    xlim([t0,t_fin]); %ylim([-params.l0/2-0.5,params.l0/2+0.5]);
     xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
      ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
    set(gca,'FontSize', fs2)
     title(tit_1,'FontSize', fs3)

   case 17 % FIG.17 Single Trace: Force on PC vs t
      %if exist('fig6','var'), figure(fig6); else fig6=figure; hold on; end  
    figure; hold on
     f_cutoff=8; % cutoff for "excessive" force, x-fold of median value
     ylab='Force on PC, a.u.';
     xlab='Time, s'; 
      tit_1=['Example: sim#' ' = ',num2str(sim), ' ',value_2_show,' = ', num2str(params.(value_2_show))]; 
      %lgnd={'F(t)',['average \tau=',num2str(n_sm*dt)],['Gauss-average \tau=',num2str(tau_sm)],['Gauss-average \tau=',num2str(tau_sm/2)]};
      lgnd={'Instantanous net force','Averaged net force', 'Averaged |Excessive force|','Averaged no Esxcessive force'};
    n_boundA=[]; force=[];
    ParB_1=ParB(t0/dt_out+1:dt/dt_out:t_fin/dt_out+1,2*totB*(sim-1)+1);
     ParA_1=ParA(t0/dt_out+1:dt/dt_out:t_fin/dt_out+1,2*totA*(sim-1)+1:2:2*totA*sim);
     ParA0_1=ParA0(t0/dt_out+1:dt/dt_out:t_fin/dt_out+1,2*totA*(sim-1)+1:2:2*totA*sim);
     ParA_state_1=ParA_state(t0/dt_out+1:dt/dt_out:t_fin/dt_out+1,totA*(sim-1)+1:totA*sim);
    F_all0=ParA_1(ParA_state_1(:)==1)-ParA0_1(ParA_state_1(:)==1);
         fm=median(abs(F_all0(:))); 
    for tt=1:size(ParA_1,1)
      for bb=1:totB
        n_boundA(tt,bb)=sum(ParA_state_1(tt,:)==bb);
        force(tt,bb)= -k_sp*sum(ParA_1(tt,ParA_state_1(tt,:)==bb)-ParA0_1(tt,ParA_state_1(tt,:)==bb));
         force1(tt,bb)= -k_sp^2*sum(ParA_1(tt,ParA_state_1(tt,:)==bb)-ParA0_1(tt,ParA_state_1(tt,:)==bb))*dt*D_parB/2;
         F_all=ParA_1(tt,ParA_state_1(tt,:)==bb)-ParA0_1(tt,ParA_state_1(tt,:)==bb);
%          fm=median(abs(F_all));        
         ind=(abs(F_all)<f_cutoff*fm);
         ind2=(abs(F_all)>=f_cutoff*fm);
        force2(tt,bb)= -k_sp*sum(F_all(ind));
         force2ex(tt,bb)= abs(k_sp*sum(F_all(ind2)));
         nf_2(tt,bb)= sum(ind2);
%         indx=find(ParA_state_1(ii,:)>0);
%         Force(ii)=-sum(ParA_1x(ii,indx)-ParA0_1x(ii,indx));
      end
    end
    for bb=1:totB
      %force_sm1(:,bb)=smooth(force(:,bb),n_sm); % force_sm(:,bb)=smooth(sign(force(:,bb)),n_sm);
      force_sm2(:,bb)=smooth_Gaussian(force2ex(:,bb),tau_sm/dt); % force_sm(:,bb)=smooth(sign(force(:,bb)),n_sm);
      force_sm3(:,bb)=smooth_Gaussian(force2(:,bb),tau_sm/dt); % force_sm(:,bb)=smooth(sign(force(:,bb)),n_sm);
      %force_sm3(:,bb)=smooth_Gaussian(nf_2(:,bb),tau_sm/dt/100); % force_sm(:,bb)=smooth(sign(force(:,bb)),n_sm);  
      force_sm1(:,bb)=smooth_Gaussian(force(:,bb),tau_sm/dt); % force_sm(:,bb)=smooth(sign(force(:,bb)),n_sm);
    end
    TT=(t0:dt:t_fin)';
    out.plot17_1_x=TT; out.plot17_1_y=force;
     out.plot17_2_x=TT; out.plot17_2_y=force2ex; %out.plot17_2_y=force_sm1;
     out.plot17_3_x=TT; out.plot17_3_y=force_sm2;
     out.plot17_4_x=TT; out.plot17_4_y=force_sm3;
     %plotting
    plot(TT,force, '-','LineWidth',lw3,'Color',cols(jj,:)); hold on
     plot(TT,force_sm1, '-','LineWidth',lw3,'Color',[0.6,0,0],'LineWidth',lw1);
     plot(TT,force_sm2, '-','LineWidth',lw3,'Color',[0,0.6,0],'LineWidth',lw1);
     plot(TT,force_sm3, '-','LineWidth',lw3,'Color',[0,0,0.6],'LineWidth',lw2);
    xlim([t0,t_fin]); %ylim([-params.l0/2-0.5,params.l0/2+0.5]);
     xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
      ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
    set(gca,'FontSize', fs2)
     title([tit_1, ' Gauss-Averaging tau=',num2str(tau_sm)],'FontSize', fs3)
     legend(lgnd, 'FontSize', fs3);
     text(0,0.5*max(force),['Exsessive force: |f| > ',num2str(f_cutoff),'*median'],'FontSize', fs3)
   
     
     
   case 18 % FIG.18 ParA distribution around PC: 1 simulation
     %if exist('fig18','var'), figure(fig18); else fig18=figure; hold on; end 
     figure; hold on
     ylab='ParA (ParB)'; xlab='Distance to PC, \mum'; 
      tit_1=['Example: sim#' ' = ',num2str(sim), ' ',value_2_show,' = ', num2str(params.(value_2_show))];
      lgnd_3={'All ParA','DNA-bound ParA','PC-bound ParA','ParB (20x)'}
     ParB_1=ParB(t0/dt_out+1:t_fin/dt_out+1,2*totB*(sim-1)+1:2:2*totB*sim);
      ParA_1=ParA(t0/dt_out+1:t_fin/dt_out+1,2*totA*(sim-1)+1:2:2*totA*sim);
      ParA_state_1=ParA_state(t0/dt_out+1:t_fin/dt_out+1,totA*(sim-1)+1:totA*sim);
     bins=(-l0/2+l0/n_bins/2:l0/n_bins:l0/2);
      aa_hist=zeros(size(bins))'; aa_hist0=aa_hist; aa_hist1=aa_hist;
      bb_hist=aa_hist;
     for bb=1:totB
       ParA_1B=ParA_1-repmat(ParB_1(:,bb),1,totA);
        aa_hist=aa_hist+histc(ParA_1B(:),bins);
         ind0=(ParA_state_1(:)==0); aa_hist0=aa_hist0+histc(ParA_1B(ind0),bins);
         ind1=(ParA_state_1(:)>0); aa_hist1=aa_hist1+histc(ParA_1B(ind1),bins); 
       ParB_1B=ParB_1-repmat(ParB_1(:,bb),1,totB);
        bb_hist=bb_hist+histc(ParB_1B(:),bins);
     end
     XX=bins(1:end-1)+diff(bins)/2;
      ind0=(ParA_state_1(:)==0); ind1=(ParA_state_1(:)>0); ind=(ParA_state_1(:)<0);
      n_tpoints=length(t0/dt_out+1:t_fin/dt_out+1);
      disp([sum(ind0),sum(ind1),sum(ind),sum(ind0)+sum(ind1)+sum(ind)]/length(ParA_state_1(:)))
       %plotting
     plot(XX,aa_hist(1:end-1)/n_tpoints, '-','LineWidth',lw2,'Color',cols(jj,:))
      plot(XX,aa_hist0(1:end-1)/n_tpoints, ':','LineWidth',lw2,'Color',cols(jj,:))
      plot(XX,aa_hist1(1:end-1)/n_tpoints, '--','LineWidth',lw2,'Color',cols(jj,:))
      plot(XX,bb_hist(1:end-1)/n_tpoints*20, ':','LineWidth',lw2,'Color','k')
       xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
       ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
      set(gca,'FontSize', fs2)
       title(tit_1,'FontSize', fs3)
       legend(lgnd_3, 'FontSize', fs3);
       
   case 19 % FIG.19 ParA distribution around PC: all simulations
     if exist('fig19','var'), figure(fig19); else fig19=figure; hold on; end  
     ylab='ParA' ; xlab='Distance to PC, \mum'; 
      tit_1=['All simulations averaged ',value_2_show,' = ', num2str(params.(value_2_show))];
      lgnd_3={'All ParA','DNA-bound ParA','PC-bound ParA','ParB (20x)'};
     bins=(-l0/2+l0/n_bins/2:l0/n_bins:l0/2);
      aa_hist=zeros(size(bins))'; aa_hist0=aa_hist; aa_hist1=aa_hist;
      bb_hist=aa_hist;
     for sim=1:n_runs 
       ParB_1=ParB(t0/dt_out+1:t_fin/dt_out+1,2*totB*(sim-1)+1:2:2*totB*sim);
        ParA_1=ParA(t0/dt_out+1:t_fin/dt_out+1,2*totA*(sim-1)+1:2:2*totA*sim);
        ParA_state_1=ParA_state(t0/dt_out+1:t_fin/dt_out+1,totA*(sim-1)+1:totA*sim);
       for bb=1:totB
         ParA_1B=ParA_1-repmat(ParB_1(:,bb),1,totA);
          aa_hist=aa_hist+histc(ParA_1B(:),bins);
           ind0=(ParA_state_1(:)==0); aa_hist0=aa_hist0+histc(ParA_1B(ind0),bins);
           ind1=(ParA_state_1(:)>0); aa_hist1=aa_hist1+histc(ParA_1B(ind1),bins); 
         ParB_1B=ParB_1-repmat(ParB_1(:,bb),1,totB);
          bb_hist=bb_hist+histc(ParB_1B(:),bins);
       end
     end
     XX=bins(1:end-1)+diff(bins)/2;
      ind0=(ParA_state_1(:)==0); ind1=(ParA_state_1(:)>0); ind=(ParA_state_1(:)<0);
      n_tpoints=length(t0/dt_out+1:t_fin/dt_out+1);
      disp([sum(ind0),sum(ind1),sum(ind),sum(ind0)+sum(ind1)+sum(ind)]/length(ParA_state_1(:)))
       %plotting
     plot(XX,aa_hist(1:end-1)/n_tpoints/n_runs, '-','LineWidth',lw2,'Color',cols(jj,:))
      plot(XX,aa_hist0(1:end-1)/n_tpoints/n_runs, ':','LineWidth',lw2,'Color',cols(jj,:))
      plot(XX,aa_hist1(1:end-1)/n_tpoints/n_runs, '--','LineWidth',lw2,'Color',cols(jj,:))
      plot(XX,bb_hist(1:end-1)/n_tpoints/n_runs*20, ':','LineWidth',lw2,'Color','k')
       xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
       ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
      set(gca,'FontSize', fs2)
       title(tit_1,'FontSize', fs3)
       legend(lgnd_3, 'FontSize', fs3);   
   
   case 20 % FIG.20 ParA distribution around PC relative to dx: all simulations  !! changed for varied dim
     if exist('fig20','var'), figure(fig20); else fig20=figure; hold on; end  
     ylab='ParA/totA'; xlab='Distance to PC, \mum'; 
      tit_1=['Relative to PC step ']; %,value_2_show,' = ', num2str(params.(value_2_show))];
      %lgnd_3={'All ParA','DNA-bound ParA','PC-bound ParA','ParB (20x)'};
      lgnd_4{2*jj-1}=[value_2_show,' = ', num2str(params.(value_2_show)),', DNA-bound'];
      lgnd_4{2*jj}=[value_2_show,' = ', num2str(params.(value_2_show)),', PC-bound'];
     bins=(-l0/2+l0/n_bins/2:l0/n_bins:l0/2);
      aa_hist=zeros(size(bins))'; aa_hist0=aa_hist; aa_hist1=aa_hist;
      bb_hist=aa_hist;
     for sim=1:n_runs 
       ParB_1=ParB(t0/dt_out+1:t_fin/dt_out+1,dim*totB*(sim-1)+1:dim:dim*totB*sim);
        ParA_1=ParA(t0/dt_out+1:t_fin/dt_out+1,dim*totA*(sim-1)+1:dim:dim*totA*sim);
        ParA_state_1=ParA_state(t0/dt_out+1:t_fin/dt_out+1,totA*(sim-1)+1:totA*sim);
         ParA_state_1cut=ParA_state_1(1+dt/dt_out:end,:);
       for bb=1:totB
         dx_ParB_1=diff_n(ParB_1(:,bb),dt/dt_out);  
         %ParA_1B=ParA_1(1:end-dt/dt_out,:)-repmat(ParB_1(1:end-dt/dt_out,bb),1,totA);
         ParA_1B=ParA_1(1+dt/dt_out:end,:)-repmat(ParB_1(1+dt/dt_out:end,bb),1,totA);
          ParA_1Bsign=ParA_1B.*repmat(sign(heaviside(dx_ParB_1)-0.2),1,totA);
          aa_hist=aa_hist+histc(ParA_1Bsign(:),bins);
           ind0=(ParA_state_1cut(:)==0); aa_hist0=aa_hist0+histc(ParA_1Bsign(ind0),bins);
           ind1=(ParA_state_1cut(:)>0); aa_hist1=aa_hist1+histc(ParA_1Bsign(ind1),bins); 
         ParB_1B=ParB_1(1:end-dt/dt_out,:)-repmat(ParB_1(1:end-dt/dt_out,bb),1,totB);
          ParB_1Bsign=ParB_1B.*repmat(sign(dx_ParB_1),1,totB);
          bb_hist=bb_hist+histc(ParB_1Bsign(:),bins);
       end
     end
     XX=bins(1:end-1)+diff(bins)/2;
      ind0=(ParA_state_1(:)==0); ind1=(ParA_state_1(:)>0); ind=(ParA_state_1(:)<0);
      n_tpoints=length(t0/dt_out+1:t_fin/dt_out+1);
      disp([sum(ind0),sum(ind1),sum(ind),sum(ind0)+sum(ind1)+sum(ind)]/length(ParA_state_1(:)))
       %plotting
      YY0=aa_hist0(1:end-1)/n_tpoints/totA/n_runs; 
       YY1=aa_hist1(1:end-1)/n_tpoints/totA/n_runs; 
      plot(XX,YY0, ':','LineWidth',lw2,'Color',cols(jj,:))
      plot(XX,YY1, '-','LineWidth',lw2,'Color',cols(jj,:))
      out.plot20_x=XX;
       out.plot20_y1=YY0;
       out.plot20_y2=YY1;
       xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
       ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
      set(gca,'FontSize', fs2)
       title(tit_1,'FontSize', fs3)        
       legend(lgnd_4, 'FontSize', fs3);
     
   case 21 % FIG.21 ParA amount(concentration) vs distance between adjacent PCs
     if exist('fig21','var'), figure(fig21); else fig21=figure; hold on; end     
     ylab='ParA'; xlab='Distance between PCs, \mum'; 
      tit_1=['All simulations ',value_2_show,' = ', num2str(params.(value_2_show))];
      lgnd_3={'inter-PC + PC-edge','PC-edge','only inter-PC'};
     dPC_01_all=[]; dPC_mid_all=[]; countsA_01_all=[]; countsA_mid_all=[];
     for sim=1:n_runs 
       ParB_1=ParB(t0/dt_out+1:t_fin/dt_out+1,2*totB*(sim-1)+1:2:2*totB*sim);
        ParA_1=ParA(t0/dt_out+1:t_fin/dt_out+1,2*totA*(sim-1)+1:2:2*totA*sim);
        ParA_state_1=ParA_state(t0/dt_out+1:t_fin/dt_out+1,totA*(sim-1)+1:totA*sim);
       ParB_1sort=sort(ParB_1,2);
       for tt=1:length(ParB_1)
         x_PC0=ParB_1sort(tt,:); x_PC=[-l0/2,x_PC0,l0/2];  
          dPC=diff(x_PC);
         bins=(x_PC);
         countsA=histc(ParA_1(tt,ParA_state_1(tt,:)==0),bins);
          dPC_01=[dPC(1),dPC(end)]; countsA_01=[countsA(1),countsA(end-1)];
          dPC_mid=dPC(2:end-1); countsA_mid=[countsA(2:end-2)];
         dPC_01_all=[dPC_01_all,dPC_01]; dPC_mid_all=[dPC_mid_all,dPC_mid];
          countsA_01_all=[countsA_01_all,countsA_01]; countsA_mid_all=[countsA_mid_all,countsA_mid];
       end
     end
     [XX,YY]=bin2_fixN([dPC_01_all,dPC_mid_all]',[countsA_01_all,countsA_mid_all]',n_bins_3);
      [XX_01,YY_01]=bin2_fixN(dPC_01_all',countsA_01_all',n_bins_3);
      [XXmid,YYmid]=bin2_fixN(dPC_mid_all',countsA_mid_all',n_bins_3);
       %plotting
% %      plot(XX,YY, 'o-','LineWidth',lw2,'Color',cols(jj,:))
% %       plot(XX_01,YY_01, 'x:','LineWidth',lw2,'Color',cols(jj,:))
% %       plot(XXmid,YYmid, '+--','LineWidth',lw2,'Color',cols(jj,:))
%      plot(XX,YY/totA, 'o-','LineWidth',lw2,'Color',cols(jj,:))
      plot(XX_01,YY_01/totA, 'x-','LineWidth',lw2,'Color',cols(jj,:))
      % to plot with equidistant binning
     [XX,YY]=bin2([dPC_01_all,dPC_mid_all]',[countsA_01_all,countsA_mid_all]',0:l0/2/n_bins:l0);
      [XX_01,YY_01]=bin2(dPC_01_all',countsA_01_all',0:l0/2/n_bins:l0);
      [XXmid,YYmid]=bin2(dPC_mid_all',countsA_mid_all',0:l0/2/n_bins:l0);
      plot(XX_01,YY_01/totA, 'o-','LineWidth',lw2,'Color',cols(jj,:))
%       plot(XXmid,YYmid/totA, '+--','LineWidth',lw2,'Color',cols(jj,:))
     xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
      ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
      set(gca,'FontSize', fs2)
       title(tit_1,'FontSize', fs3)
       legend(lgnd_3, 'FontSize', fs3); 
     % one more figure: ParA concnetration  
     if exist('fig21a','var'), figure(fig21a); else fig21a=figure; hold on; end
     ylab='ParA concentration';
     [XX,YY]=bin2_fixN([dPC_01_all,dPC_mid_all]',[countsA_01_all,countsA_mid_all]'./[dPC_01_all,dPC_mid_all]',n_bins_3);
      [XX_01,YY_01]=bin2_fixN(dPC_01_all',countsA_01_all'./dPC_01_all',n_bins_3);
      [XXmid,YYmid]=bin2_fixN(dPC_mid_all',countsA_mid_all'./dPC_mid_all',n_bins_3);
% %       plot(XX,YY, 'o-','LineWidth',lw2,'Color',cols(jj,:))
% %       plot(XX_01,YY_01, 'x:','LineWidth',lw2,'Color',cols(jj,:))
% %       plot(XXmid,YYmid, '+--','LineWidth',lw2,'Color',cols(jj,:))
%       plot(XX,YY/totA, 'o-','LineWidth',lw2,'Color',cols(jj,:))
      plot(XX_01,YY_01/totA, 'x:','LineWidth',lw2,'Color',cols(jj,:))
%       plot(XXmid,YYmid/totA, '+--','LineWidth',lw2,'Color',cols(jj,:))
      % to plot with equidistant binning
     [XX,YY]=bin2([dPC_01_all,dPC_mid_all]',[countsA_01_all,countsA_mid_all]',0:l0/2/n_bins:l0);
      [XX_01,YY_01]=bin2(dPC_01_all',countsA_01_all'./dPC_01_all',0:l0/2/n_bins:l0);
      [XXmid,YYmid]=bin2(dPC_mid_all',countsA_mid_all'./dPC_mid_all',0:l0/2/n_bins:l0);
      plot(XX_01,YY_01/totA, 'o-','LineWidth',lw2,'Color',cols(jj,:))
     xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
      ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
      set(gca,'FontSize', fs2)
       title(tit_1,'FontSize', fs3)
       %legend(lgnd_3, 'FontSize', fs3);
       legend(lgnd_1, 'FontSize', fs3); title(tit_0,'FontSize', fs3) 
                
   case 22 % Fig.22: Next step correlation vs inter-PC_distance
     if exist('fig22','var'), figure(fig22); else fig22=figure; hold on; end
      ylab='Next Step Correlation'; xlab='Distance to the closest PC, \mum';
      tit=['dt=',num2str(dt),' binned by inter-PC distance, 1 bin=',num2str(n_bins_3),'steps ']; 
     dx1_all=[]; dx2_all=[]; distB_all=[];
      XX=[]; YY=[]; 
     for sim=1:n_runs
       ParB_1=ParB(t0/dt_out+1:t_fin/dt_out+1,2*totB*(sim-1)+1:2:2*totB*sim);
        ParA_1=ParA(t0/dt_out+1:t_fin/dt_out+1,2*totA*(sim-1)+1:2:2*totA*sim);
        ParA_state_1=ParA_state(t0/dt_out+1:t_fin/dt_out+1,totA*(sim-1)+1:totA*sim);
       for bb=1:totB
         dX=diff_n(ParB_1(:,totB),dt/dt_out);
          ind=1:size(dX,1)-1;
         dx_ParB_1=dX(ind,:);
          dx_ParB_2=dX(ind+1,:);
          distB=abs(ParB_1(ind+1,:)-repmat(ParB_1(ind+1,bb),1,totB));
          distB=min(distB+dirac(distB),[],2); 
         dx1_all=[dx1_all;dx_ParB_1]; 
          dx2_all=[dx2_all;dx_ParB_2]; 
          distB_all=[distB_all;distB];
       end
     end  
     dBdX1dX2=sortrows([distB_all,dx1_all,dx2_all],1);
     for nn=1:n_bins_3:length(dBdX1dX2)-n_bins_3+1
       XX=[XX;mean(dBdX1dX2(nn:nn+n_bins_3-1,1))]; YY=[YY;corr(dBdX1dX2(nn:nn+n_bins_3-1,2),dBdX1dX2(nn:nn+n_bins_3-1,3))];  
     end
      if nn+n_bins_3-1<length(dBdX1dX2)
        XX=[XX;mean(dBdX1dX2(nn+n_bins_3:end,1))]; YY=[YY;corr(dBdX1dX2(nn+n_bins_3:end,2),dBdX1dX2(nn+n_bins_3:end,3))];
      end
     plot(XX,YY, 'o','MarkerFaceColor',cols(jj,:),'MarkerEdgeColor',cols(jj,:),'LineWidth',lw2)
      xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
       ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
      set(gca,'FontSize', fs2)
       legend(lgnd_1, 'FontSize', fs3); title(tit_0,'FontSize', fs3) 
        title(tit,'FontSize', fs3)  
        
   case 23 % Fig.23: ensemble-average MSD(abs(x))
     if exist('fig23','var'), figure(fig23); else fig23=figure; hold on; end  
     ylab='MSD, \mum^2'; xlab='Time, s'; 
     ParBcut=ParB(t0/dt_out+1:t1/dt_out+1,1:2:end);
     dX=abs(ParBcut(:,:)-repmat(ParBcut(1,:),size(ParBcut,1),1));
     MSD=mean(dX.^2,2);
     
    TT=(t0:dt_out:t1)';
     plot(TT, MSD,'o-','Color',cols(jj,:),'MarkerFaceColor',cols(jj,:),'MarkerEdgeColor',cols(jj,:),'LineWidth',lw2,'MarkerSize',ms2) 
       xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
        ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
       set(gca,'FontSize', fs2)
      legend(lgnd_1, 'FontSize', fs3); title(tit_0,'FontSize', fs3)  
      
   case 24 % FIG.24 ParA distribution around PC relative to dx: at given time points
     if exist('fig24','var'), figure(fig24); else fig24=figure; hold on; end     
     %tt_all0=[0,1,5,50,500]; % time-points at which  to plot profile 
     tt_all0=500; % time-points at which  to plot profile 
      tt_all=tt_all0/dt_out+1; % time-points at which  to plot profile 
     cols_w0=winter(length(tt_all)); cols_w(:,3)=2*(cols_w0(:,3)-0.5); cols_w(:,2)=0*cols_w0(:,2); cols_w(:,1)=cols_w0(:,2);  %winter(length(tt_all));
      for ii=1:length(tt_all), lgnd_5{ii}=['t = ',num2str(tt_all(ii)-1)]; end
       
     ylab='ParA/totA'; xlab='Distance to PC, \mum'; 
      tit_1=['Relative to PC step ']; %,value_2_show,' = ', num2str(params.(value_2_show))];
     lgnd_3={'DNA-bound ParA','PC-bound ParA'};
%       lgnd_4{2*jj-1}=[value_2_show,' = ', num2str(params.(value_2_show)),', DNA-bound'];
%       lgnd_4{2*jj}=[value_2_show,' = ', num2str(params.(value_2_show)),', PC-bound'];
     bins=(-l0/2+l0/n_bins/2:l0/n_bins:l0/2);
     for ii=1:length(tt_all)
         %figure; hold on  
         tt=tt_all(ii);
         aa_hist=zeros(size(bins))'; aa_hist0=aa_hist; aa_hist1=aa_hist;
%           bb_hist=aa_hist;
         for sim=1:n_runs 
           ParB_1=ParB(t0/dt_out+1:t_fin/dt_out+1,2*totB*(sim-1)+1:2:2*totB*sim);
            ParA_1=ParA(t0/dt_out+1:t_fin/dt_out+1,2*totA*(sim-1)+1:2:2*totA*sim);
            ParA_state_1=ParA_state(t0/dt_out+1:t_fin/dt_out+1,totA*(sim-1)+1:totA*sim);
             ParA_state_1cut=ParA_state_1(tt,:);
           for bb=1:totB
             dx_ParB_1=diff_n(ParB_1(:,bb),dt/dt_out);
             dx_ParB_1=dx_ParB_1(tt);
             %ParA_1B=ParA_1(1:end-dt/dt_out,:)-repmat(ParB_1(1:end-dt/dt_out,bb),1,totA);
             ParA_1B=ParA_1(tt,:)-repmat(ParB_1(tt,bb),1,totA);
              ParA_1Bsign=ParA_1B.*repmat(sign(heaviside(dx_ParB_1)-0.2),1,totA);
              aa_hist=aa_hist+histc(ParA_1Bsign(:),bins);
               ind0=(ParA_state_1cut(:)==0); aa_hist0=aa_hist0+histc(ParA_1Bsign(ind0),bins)';
               ind1=(ParA_state_1cut(:)>0); aa_hist1=aa_hist1+histc(ParA_1Bsign(ind1),bins)'; 
%              ParB_1B=ParB_1(tt,:)-repmat(ParB_1(tt,bb),1,totB);
%               ParB_1Bsign=ParB_1B.*repmat(sign(dx_ParB_1),1,totB);
%               bb_hist=bb_hist+histc(ParB_1Bsign(:),bins);
           end
         end
         XX=bins(1:end-1)+diff(bins)/2;
%           ind0=(ParA_state_1(:)==0); ind1=(ParA_state_1(:)>0); ind=(ParA_state_1(:)<0);
%           n_tpoints=length(t0/dt_out+1:t_fin/dt_out+1);
%           disp([sum(ind0),sum(ind1),sum(ind),sum(ind0)+sum(ind1)+sum(ind)]/length(ParA_state_1(:)))
           %plotting
          col=cols_w(ii,:);
           YY0=aa_hist0(1:end-1)/totA/n_runs;
           YY1=aa_hist1(1:end-1)/totA/n_runs;
          plot(XX,YY0, ':','LineWidth',lw2,'Color',col)
           plot(XX,YY1, '-','LineWidth',lw2,'Color',col)
          out.plot24_x{ii}=XX;
           out.plot24_y1{ii}=YY0;
           out.plot24_y2{ii}=YY1;
          %plot(XX,(aa_hist1(1:end-1)+aa_hist0(1:end-1))/totA/n_runs, '-','LineWidth',lw2,'Color',col)
          xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
           ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
          set(gca,'FontSize', fs2)
          title([tit_1,'  t = ' num2str(tt)],'FontSize', fs3)        
          %legend(lgnd_5, 'FontSize', fs3);
      end
       xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
       ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
      set(gca,'FontSize', fs2)
       title(tit_1,'FontSize', fs3)        
       legend(lgnd_5, 'FontSize', fs3);
      xlim([-1,1]); 
       
   case 25 % FIG.25 Average Force vs distance between adjacent PCs
     if exist('fig25','var'), figure(fig25); else fig25=figure; hold on; end     
      ylab='Force, a.u.'; xlab='Distance between PCs, \mum'; 
       tit_1=['All simulations ',value_2_show,' = ', num2str(params.(value_2_show))];
     dy_max=0.2; % to sort out data with PC "on one line"  
      lgnd={'All points',['only |\delta(Y)|<', num2str(dy_max)],['binned with ',num2str(n_bins_3), ' per bin']}; 
     % calculating forces and distances
     forces_all=[]; PC_dist_all=[]; PC_distY_all=[];
     for sim=1:n_runs
       ParA_1=ParA(t0/dt_out+1:t_fin/dt_out+1,2*totA*(sim-1)+1:2:2*totA*sim);
        ParA0_1=ParA0(t0/dt_out+1:t_fin/dt_out+1,2*totA*(sim-1)+1:2:2*totA*sim);
        ParA_state_1=ParA_state(t0/dt_out+1:t_fin/dt_out+1,totA*(sim-1)+1:totA*sim);
       ParB_1=ParB(t0/dt_out+1:t_fin/dt_out+1,2*totB*(sim-1)+1:2:2*totB*sim);
        ParB_1y=ParB(t0/dt_out+1:t_fin/dt_out+1,2*totB*(sim-1)+2:2:2*totB*sim);
        %ParB_1sort=sort(ParB_1,2);
       for tt=1:size(ParA_1,1)
         for bb=1:totB
%            x_PC0=ParB_1sort(tt,:); x_PC=[-l0/2,x_PC0,l0/2];  
%            dPC=diff(x_PC);
%            bins=(x_PC);
%            countsA=histc(ParA_1(tt,ParA_state_1(tt,:)==0),bins);
%             dPC_01=[dPC(1),dPC(end)]; countsA_01=[countsA(1),countsA(end-1)];
%           dPC_mid=dPC(2:end-1); countsA_mid=[countsA(2:end-2)];
%          dPC_01_all=[dPC_01_all,dPC_01]; dPC_mid_all=[dPC_mid_all,dPC_mid];   
             
          %n_boundA(tt,bb)=sum(ParA_state_1(tt,:)==bb);
          bb_but_one=setdiff(1:totB,bb);
          [dist,ind]=min(abs(ParB_1(tt,bb_but_one)-ParB_1(tt,bb)));
           PC_dist(tt,bb)=dist;
           PC_distY(tt,bb)=abs(ParB_1y(tt,bb_but_one(ind))-ParB_1y(tt,bb));
          force_1= -k_sp*sum(ParA_1(tt,ParA_state_1(tt,:)==bb)-ParA0_1(tt,ParA_state_1(tt,:)==bb));
           force(tt,bb)= force_1*sign(ParB_1(tt,bb_but_one(ind))-ParB_1(tt,bb));
%         indx=find(ParA_state_1(ii,:)>0);
%         Force(ii)=-sum(ParA_1x(ii,indx)-ParA0_1x(ii,indx));
         end
       end
       forces_all=[forces_all; force(:)];
        PC_dist_all=[PC_dist_all; PC_dist(:)]; 
        PC_distY_all=[PC_distY_all; PC_distY(:)];  
     end
      % to plot with equidistant binning
     [XX,YY]=bin2(PC_dist_all,forces_all,0:l0/2/n_bins:l0); 
      [XX1,YY1]=bin2(PC_dist_all(PC_distY_all<dy_max),forces_all(PC_distY_all<dy_max),0:l0/2/n_bins:l0); 
      % to plot with equi_number binning
     [XX2,YY2]=bin2_fixN(PC_dist_all,forces_all,n_bins_3);
      %plotting
     plot(XX,YY, '.-','LineWidth',lw3); 
      hold on
      plot(XX1,YY1, 'o-','LineWidth',lw3);
      plot(XX2,YY2, 'x-','LineWidth',lw3);
     xlim([0,1]); %ylim([-params.l0/2-0.5,params.l0/2+0.5]);
      xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
       ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
     set(gca,'FontSize', fs2)
      title(tit_1,'FontSize', fs3)
      legend(lgnd, 'FontSize', fs3);
      
%      figure; hold on
%      scatter(PC_dist_all,forces_all,9,PC_distY_all)

   case 26 % FIG.26 average 26ParA distribution: at given time points
     %tt_all0=[2990,3600,4200,4800,5400,5990]; % time-points at which  to plot profile 
     %tt_all0=[10,2990,3300,4000,5000,5990]; % 1-st cycle
     % tt_all0=[6010,8990,9300,10000,11000,11990]; % 2-nd cycle
     %tt_all0=[12010,14990,15300,16000,17000,17990]; % 3-d cycle
     %tt_all0=[18010,20990,21300,22000,23000,23990]; % 4-th cycle
     %tt_all0=[24010,26990,27300,28000,29000,29990]; % 5-th cycle
     %tt_all0=[30010,32990,33300,34000,35000,35990]; % 6-th cycle
     %tt_all0=[36010,38990,39300,40000,41000,41990]; % 7-th cycle
     %tt_all0=[42010,44990,45300,46000,47000,47990]; % 8-th cycle
     %tt_all0=[48010,50990,51300,52000,53000,53990]; % 9-th cycle
     %tt_all0=[54010,56990,57300,58000,59000,59990]; % 10-th cycle
     tt_all0=[0,6000,12000,18000,24000,30000,36000,42000,48000,54000,60000]; % start and end of replicatin for 10 cycles     
     tt_all=tt_all0/dt_out+1; % time-points at which  to plot profile 
     cols_w=winter(length(tt_all)); % jet(length(tt_all));
      %for ii=1:length(tt_all), lgnd_5{ii}=['t = ',num2str(tt_all0(ii))]; end
      for ii=1:length(tt_all), lgnd_5{ii}=['cycle #',num2str(ii-1)]; end
     ylab='ParA concentration'; xlab='Relative cell coordinate'; 
      tit_1=[tit_0,'  ', value_2_show,' = ', num2str(params.(value_2_show))];
     fig_pos= [200,200,300,600]; 
     %lgnd_3={'DNA-bound ParA','PC-bound ParA'};
%       lgnd_4{2*jj-1}=[value_2_show,' = ', num2str(params.(value_2_show)),', DNA-bound'];
%       lgnd_4{2*jj}=[value_2_show,' = ', num2str(params.(value_2_show)),', PC-bound'];
     bins=(-l0/2+l0/n_bins/2:l0/n_bins:l0/2);
     figure; hold on 
%      figure('Position',fig_pos); hold on
     ymax=0;
     n_tt=length(tt_all);
     for ii=1:n_tt
        tt=tt_all(ii);
        ParA_1n=ParA(tt,1:2:2*totA*n_runs);
         ParA_state_1n=ParA_state(tt,1:totA*n_runs);
             %ParA_state_1cut=ParA_state_1(tt,:);
         aa_hist=histc(ParA_1n(:),bins); aa_hist2=l0*aa_hist/totA/n_runs/(bins(2)-bins(1));  
        XX=bins(1:end-1)/l0+1/2;
         YY=aa_hist2(1:end-1);
%         subplot(n_tt,1,ii,'Position',[0.16,1-0.06-0.86/n_tt*ii,0.82,0.85/n_tt]) 
        %plot(XX, -ymax+YY-(ii-1)*0,'-','Color',cols_w(ii,:),'LineWidth',lw2)        
        plot(XX, YY-(ii-1)*0,'-','Color',cols_w(ii,:),'LineWidth',lw2)
%          text(0.7,2.3,['t = ',num2str(tt_all0(ii)) ' s'],'FontSize', fs3) ;
%           ylim([0,2.5]);
%           set(gca,'FontSize', fs2,'XTickLabel',[])
          if ii==1, title([tit_0,' ', value_2_show,' = ', num2str(params.(value_2_show))], 'FontSize', fs3); end
        out.plot26_x{ii}=XX;
         out.plot26_y1{ii}=YY;
         ymax=ymax+1.05*max(YY);
     end
          xlim([0,1]); %ylim([-params.l0/2-0.5,params.l0/2+0.5]);
      set(gca,'FontSize', fs2)
      xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
       ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
      legend(lgnd_5, 'FontSize', fs3);
      title(tit_1,'FontSize', fs3)
      
   case 27 % Fig.9: MSD-tau for ParA
     if exist('fig27','var'), figure(fig27); else fig27=figure; txt={}; lgnd_27={}; hold on; end  
     ylab='MSD_{\tau} for ParA, \mum'; xlab='Time, s'; 
    
     ParAcut=ParA(t0/dt_out+1:t1/dt_out+1,1:2:end);
     for ii=0:ind_max
        ind1=1:size(ParAcut,1)-ii; ind2=ind1+ii;
        dX=ParAcut(ind2,:)-ParAcut(ind1,:);
        MSD(ii+1)=sum((dX(:)).^2)/sum(ones(size(dX(:))));
     end
    TT=dt_out*(0:ind_max)';
     YY1=MSD(TT<t2)'; TT1=TT(TT<t2);
    [coeff,fit1,YYfit]=oneFit(TT1,YY1,'lin_growth'); disp(fit1)
     ci = confint(fit1); %TT11=TT2(2:end); 
     fit_tit='fit by a*t+b';
    coeff2(1)=coeff(2); coeff2(1)=coeff(2); %coeff=coeff1;
     ci2(:,1)=ci(:,2); ci2(:,2)=ci(:,1); ci=ci2;
    YYfit=coeff(1)*TT1+coeff(2);
    lgnd_27{2*jj-1}=lgnd_1{jj}; lgnd_27{2*jj}=[' D=',num2str(coeff(1)/2,2)];
    out.plot27(:,jj)=[coeff(1),ci(:,1)',coeff(2),ci(:,2)'];
     disp(num2str(coeff,3)); disp(num2str(ci2,3))
    plot(TT, MSD,'o','MarkerFaceColor',cols(jj,:),'MarkerEdgeColor',cols(jj,:),'LineWidth',lw2)
     plot(TT1,YYfit,':k','LineWidth',lw2)
       xlabel(xlab, 'FontSize', fs1); 
        ylabel(ylab, 'FontSize', fs1); 
       set(gca,'FontSize', fs2)
      legend(lgnd_27, 'FontSize', fs3); title(tit_0,'FontSize', fs3) 
      %text(0.6*x_lim(2),y_lim(1)+0.15*(y_lim(2)-y_lim(1)),txt,'FontSize', fs3)
    
      
end
end


end

% disp(['no PC | %',' | Stats {mean ','| std   ',' |  min   ',' | max }',' of (#PC in DaugtherCell)' ]);
% disp(['at t = ', num2str(t1)])
% disp(num2str(Stats,4));
% disp(['summed over t = ', num2str(t0),'-', num2str(t1)])
% disp(num2str(Stats2,4));
end