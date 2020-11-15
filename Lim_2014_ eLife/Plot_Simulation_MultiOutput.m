% to plot results of simulation from MultiOutput data structure 
% - repetaed simulations with varied parameters
%1
% change 111 to 0 for data prior to 11.21.2013 (line 83: ind=find(trajectory(:,1)==111,1,'first');)


% data output specifications
n_col=1; % maximum number of panels per row
switch_params=0; % 1, if want to switch param1 and param2 
dt_step=1; % how often to plot data, 1 - every time step, 2 every 2nd etc
aa=0.98; % decrease if some panels disappear

% plotting specifications

%col1=[0,0,1];
 shade1=0; %0.3; %col2=(1-shade1)*[1,1,1]+shade1*col1;
ls1='-'; ls2=':'; 
 lw1=3; lw2=2;
mrk1='none'; mrk2='+';
 ms1=12; ms2=5;
fsiz1=16; fsiz2=18; fsiz0=15;
y_max=MultiOutput.Data{1}.params.l0/2; y_min=-1*y_max;
 x_max=MultiOutput.Data{1}.params.t_fin; x_min=0;

% get params names and values
sset_vs_values=MultiOutput.sset_vs_values;
n1=2; n2=3;
params=strsplit(MultiOutput.varied_params,' & ');
 param_1=params{1}; param_2=params{2};
% was before, worked: 
%ii1=find(MultiOutput.sset_vs_values(:,3)==MultiOutput.sset_vs_values(1,3),2,'first');
ii1=find(MultiOutput.sset_vs_values(:,2)==MultiOutput.sset_vs_values(1,2),1,'last');
 values_2=MultiOutput.sset_vs_values(1:ii1,n2); 
 values_1=MultiOutput.sset_vs_values(1:ii1:end,n1); 
if switch_params
  values_10=values_1; param_10=param_1; 
   param_1=param_2; values_1=values_2; n1=3;
   param_2=param_10; values_2=values_10; n2=2;
end
nn_1=length(values_1); nn_2=length(values_2); n_row=ceil(nn_2/n_col); n_rowd=double(n_row); n_cold=double(n_col);
% change param names so they will be displayed correctly on plots
if ~isempty(strtok(param_2,'_'))
  [param_2,remain]=strtok(param_2,'_');
  param_2=[param_2,'_{',remain(2:end),'}'];
end
if ~isempty(strtok(param_1,'_'))
  [param_1,remain]=strtok(param_1,'_');
  param_1=[param_1,'_{',remain(2:end),'}'];
end

 scrsz = get(0,'ScreenSize');
 figure('Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(3)/3]);
 
 
 %cols=jet(nn_1);
 %cols=colorcube(nn_1);
 cols=hsv(nn_1);
 %cols=prism(nn_1);
 %cols=lines(nn_1);
 %cols=[1,0,0; 0,1,0; 0,0,1];
 ppanel=0;
 sset=0;
 txt1=['    ',param_1];
 for kk=1:nn_1, txt11{kk}=num2str(values_1(kk),2); end
 
 % plot multiple panels, one for each param2 value
 for ii=1:nn_2
   ppanel=ii;
    ppanel1=mod(ppanel,n_col); ppanel1=n_col*(1-sign(ppanel1-0.5))/2-sign(0.5-ppanel1)*ppanel1;
    ppanel2=ceil(ppanel/n_col);
   txt2=[param_2,' = ',num2str(values_2(ppanel),3)]; 
     
   
   
  % plot one panel
   subplot(n_row,n_col,ppanel);
   sets_to_plot=sset_vs_values(sset_vs_values(:,n2)==values_2(ii),:);
 for kk=1:nn_1
 sset=sets_to_plot(kk,1);
  ParB=MultiOutput.Data{sset}.ParB;
  dt=MultiOutput.Data{sset}.params.dt_out; t_fin=MultiOutput.Data{sset}.params.t_fin;
  
% DATA PREPARATION
 % to substistute zeros in trajectories with last non-zero values
 for jj=1:size(ParB,2)/2
   trajectory=ParB(:,2*jj-1:2*jj);
    tr_length=length(trajectory);
   ind=find(trajectory(:,1)==111,1,'first');
   if ~isempty(ind)&& ind>1
     trajectory(ind:end,:)=repmat(trajectory(ind-1,:),tr_length-ind+1,1);
   end
   ParB(:,2*jj-1:2*jj)=trajectory;
 end
 % get average trace and std 
 TT=0:dt*dt_step:t_fin; 
 TT2=TT+0.3*dt*dt_step*(-1+2*kk/nn_1);
 ParB_Xmean=mean(ParB(1:dt_step:end,1:2:end)');
  ParB_Xstd=std(ParB(1:dt_step:end,1:2:end)');
 
% PLOTTING
 col1=cols(kk,:);
 col2=(1-shade1)*[1,1,1]+shade1*col1;
  %plot([TT2;TT2],[ParB_Xmean-ParB_Xstd;ParB_Xmean+ParB_Xstd], 'LineStyle',ls2,'LineWidth',lw2,'Color',col2,'Marker',mrk2,'MarkerEdgeColor',col2,'MarkerFaceColor',col2,'MarkerSize',ms2);
  hold on
  plot(TT,ParB_Xmean, 'LineStyle',ls1,'LineWidth',lw1,'Color',col1,'Marker',mrk1,'MarkerEdgeColor',col1,'MarkerFaceColor',col1,'MarkerSize',ms1);

 end %  plot one panel

% fix panel appearance 
xlim([TT(1),TT(end)]); ylim([-y_max,y_max]);
 
% legend for colored lines
x01=0.80; dx01=0.02; % where to put lines for legend
 y01=-0.45; dy01=0.09; % where to put lines for legend
text(x01*x_max,(y01+dy01)*y_max,txt2,'FontSize',fsiz1);
 set(gca,'FontSize',fsiz1,'OuterPosition',[0.1+0.9/n_cold*(ppanel1-1),0.1+0.9-0.9/n_rowd*ppanel2,0.9*aa/n_cold,0.9*aa/n_rowd],'Position',[0.1+0.9/n_cold*(ppanel1-1),0.1+0.9-0.9/n_rowd*ppanel2,0.9*aa/n_cold,0.9*aa/n_rowd]); 
if ppanel1==1 && ppanel2==1
  for kk=1:nn_1
    col1=cols(kk,:);
    dy=0.12;
    plot([x01*x_max,(x01+dx01)*x_max],[(y01-dy01*kk)*y_max,(y01-dy01*kk)*y_max],'Color',col1,'LineStyle',ls2,'LineWidth',lw2);
    text((x01+2*dx01)*x_max,(y01-dy01*kk)*y_max,txt11{kk},'FontSize',fsiz0);
  end
  text((x01+1.5*dx01)*x_max,(y01-dy01*0)*y_max,txt1,'FontSize',fsiz1);
   
end

if ~(ppanel1==1), set(gca,'YTickLabel',{}); end   
if ~(ppanel2==n_row), set(gca,'XTickLabel',{}); end

if ppanel1==1 && ppanel2==1, ylabel('Cell coordinate, \mum'); title([param_1,' vs ',param_2]); end 

end % plot multiple panels one for each param2 value
xlabel('Time, s');
disp(MultiOutput.Data{1}.params)

