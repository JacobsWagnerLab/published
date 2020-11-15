function plot_MultiplePlots(MultiOutput, sset, pplots)

% to make various plots from simulated data:
% - "snapshots" at given moment of time for multiple traces
% - etc

%t=2000;
t0=500;
dt=10;
%sset=1;
n_bins=25; n_bins_2=200;
n_tr=1;

% plotting params
cols=['b','g','r','c','m','k']; cols2=prism(n_tr);
 ms=4;
 colormap=cols;
 lw1=3; lw2=2;
 fs1=18; fs2=16;

% getting data
Data=MultiOutput.Data{sset};
params=Data.params;
 n_runs=Data.n_runs;
 dt_out=params.dt_out;
  t_fin=params.t_fin;
 totB=params.totB;
  totA=params.totA;
 k_hyd=params.k_hyd;
  k_db=params.k_db;
 D_parB=params.D_parB;
 t=params.t_fin;
ParB=Data.ParB;

disp(['n_runs','  dt_out','    t_fin','    totB','   totA','    khyd','     k_db' ])
disp(num2str([n_runs, dt_out, t_fin, totB, totA, k_hyd, k_db ],3))


for kk=1:length(pplots)
switch pplots(kk)
  case 1 % FIG.1
     % plotting all PC positions at t=t;
    if exist('fig1'), figure(fig1); else fig1=figure; hold on; end
    ylab='Distance to pole, \mum';
     xlab='Trace #'; 
     % auxilary params
    x0=params.l0/2+0.5; % "pole" coordinate
    nB=totB;
     % sorting PCs by distance from the pole
    X=ParB(t/dt_out+1,1:2:end);
     Xcell=reshape(X,nB,n_runs);
     Xcell=Xcell';
     Xcell=sort(Xcell,2); 
     %plotting
       %plot(1:n_runs, x0-Xcell,'o','MarkerEdgeColor',cols,'MarkerFaceColor',cols,'MarkerSize',ms)
    plot(1:n_runs,x0*ones(1,n_runs), ':','Color',[0.5,0.5,0.5],'LineWidth',lw2)
     plot(1:n_runs,x0-Xcell, 'o','MarkerSize',ms,'LineWidth',lw1)
    ylim([0,params.l0+1]); 
    xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
     ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
    set(gca,'FontSize', fs2)
      
     %  Some Statitistics
    no_left=find(Xcell(:,1)>0); no_right=find(Xcell(:,totB)<0); no_PC=length(no_left)+length(no_right);
     left=sum((Xcell<0),2); lr_ratio=left/totB; 
    Stat.Fraction_wo_PC=no_PC/n_runs/2;
     Stat.lr_ratio=[mean(lr_ratio),std(lr_ratio),min(lr_ratio),max(lr_ratio)];
    disp(['cells w/o PC',' Fraction w/o PC',' Mean    ',' Std   ','  min   ','  max ',' of (#PC in DaugtherCell)' ]);
     disp([no_PC, Stat.Fraction_wo_PC,Stat.lr_ratio]);

  case 2 % FIG.2 Hist of relative positions
    figure; hold on
     ylab='Frequency';
     xlab='Position, \mum'; 
     %getting counts
    X2=ParB(t0/dt_out+1:t/dt_out+1,1:2:end);
     [HistY,HistX]=hist(X2(:),n_bins);
     HistYn=HistY/sum(HistY)/(HistX(2)-HistX(1));
     %plotting
    plot(HistX,HistYn, '-','Color',cols(1),'LineWidth',lw2)
    xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
     ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
    set(gca,'FontSize', fs2)
  
  case 3 % FIG.3 Hist of abs (relative positions)
    figure; hold on
     ylab='Frequency';
     xlab='Position, \mum'; 
     %getting counts
    [HistY,HistX]=hist(abs(X2(:)),n_bins);
     HistYn=HistY/sum(HistY)/(HistX(2)-HistX(1));
     %plotting
    plot(HistX,HistYn, '-','Color',cols(3),'LineWidth',lw2)
    xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
     ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
    set(gca,'FontSize', fs2)

  case 4 % FIG.4 Hist of inter-PC distance
    figure; hold on
     ylab='Frequency';
     xlab='Distance, \mum'; 
     %getting counts
    X3=diff(Xcell');
     [HistY,HistX]=hist(X3(:),n_bins); 
     HistYn=HistY/sum(HistY)/(HistX(2)-HistX(1));
     %plotting
    plot(HistX,HistYn, 'o-','Color',cols(2),'LineWidth',lw2)
    xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
     ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
    set(gca,'FontSize', fs2)
  
  case 5 % FIG.5 PC-steps
    figure; hold on
     ylab='Frequency';
     xlab='Displacement, \mum'; 
     %getting counts
    dX=diff(ParB(1:dt/dt_out:end,:));
     [HistY,HistX]=hist(abs(dX(:)),n_bins_2); 
     HistYn=HistY/sum(HistY)/(HistX(2)-HistX(1));
    dX_gauss=exp(-(HistX.^2/(2*2*dt*D_parB)))/sqrt(pi*dt*D_parB);
     %plotting
    plot(HistX,HistYn, 'o-','Color',cols(3),'LineWidth',lw2)
     plot(HistX,dX_gauss, '-','Color',cols(6),'LineWidth',lw2)
    xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
     ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
    set(gca,'FontSize', fs2)

  case 6 % FIG.6 Traces examples
    figure; hold on
     ylab='Position, \mum';
     xlab='Time, s'; 
     %plotting
    plot(0:dt_out:t_fin,ParB(:,1:2:2*totB), '-','LineWidth',lw2)
    xlim([0,t_fin]); ylim([-params.l0/2-0.5,params.l0/2+0.5]);
    xlabel(xlab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
     ylabel(ylab, 'FontSize', fs1); %ylabel(ylab, 'FontSize', fs1);
    set(gca,'FontSize', fs2)

end
end


end