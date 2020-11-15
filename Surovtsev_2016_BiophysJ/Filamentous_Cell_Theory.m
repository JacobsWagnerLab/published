% to make a showcase of theoretical distribution of ParA between PCs for a filamentous cell 
% PCs posiitons are chosen randomly
% Based on a theory described in Biophys.J. 2016 110 2790
% Briefly it considers 
% - diffusing ParA
% - each Partition Complex (PC) serve as a sink, from wheer ParA is
% distributed randomly over the interval (cell length)
% - there might be multiple PCs in the cell

l_cell=20; % cell length, um
 n_PC=8; % number rof PC per cell
n_cells=200; % umber of virtual cell to emulate
n_bins=50; % for median-averages calculation

% some experimental parameters
a_tot=1; % total ParA per cell;
k_hyd=0.03; % rate of ParA rebinding to DNA
 D=0.005; % ParA diffusion coefficent (apparent, on DNA surface)  
 % NOTE: only the k_hyd/D ratio matters 
 kD=k_hyd/D;

data_out=[]; %PC_positions=[];
n_ones=ones(n_PC-1,1);
for cell=1:n_cells
  % get random PC_positions  
  rand_PC=sort(l_cell*rand(1,n_PC-2));
   PC_positions=[0,rand_PC,l_cell];
  % preparing auxilary parameters  
  li=diff(PC_positions); % inter-PC distance  (from left to right)
   li3=li.^3;
   li3_sum=sum(li3);
  a0=a_tot/(1+li3_sum*kD/12/l_cell); % fraction of ParA bound to PC (or amount in a_tot units)
   ai=(kD/2)*a0*li3/l_cell/6; % fraction of total ParA in between PCs (or amount in a_tot units)
  
  % collecting output data organized in columns as: [cell#, PC posiition, interPC distance, ParA between PCs, ParA bound to PC] 
  data_out((cell-1)*(n_PC-1)+1:cell*(n_PC-1),:)=[cell*n_ones,PC_positions(2:end)',li',ai',a0*n_ones];
  
end

% to plot results
 fs1=18; fs2=16; fs3=14;
% getting median-averages 
[Xout,Yout]=bin2_fixN(data_out(:,3),data_out(:,4),n_bins,2);
figure; hold on
mrk='x'; ls='none'; lw1=2;ms1=7; col=[1,0,0.5];
 plot(data_out(:,3),data_out(:,4),'LineStyle',ls,'Marker',mrk,'MarkerSize',ms1,'MarkerFaceColor',col,'MarkerEdgeColor',col)
mrk='o'; ls='-'; lw1=2; ms1=7; col=col*0.5;
 plot(Xout,Yout,'LineStyle',ls,'LineWidth',lw1,'Color',col,'Marker',mrk,'MarkerSize',ms1,'MarkerFaceColor',col,'MarkerEdgeColor',col)
set(gca,'FontSize',fs2)
 xlabel('Inter-PC distance, \mum','FontSize',fs1)
  ylabel('ParA fraction between adjacent PCs','FontSize',fs1)
 legend({['n_{PC}=',num2str(n_PC)],'median-average'},'FontSize',fs3)

 
 
 