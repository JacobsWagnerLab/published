function [ParB,ParA,tt,params,ParB_dxymax,ParA_state,ParA_0] = diffusing_multiple_chained_ParB_oscillating_ParA_v4(params)
% UNDER CONSTRUCTION 
% based on diffusing_multiple_ParB_oscillating_ParA_v4_2014_04_23 (version='04.23.2014')
% parB are "chained", i.e connected by ideally flexible links with fixed
% length (i.e. distance between adjacent ParB <= l_max, there is no elastic, or other, force)
%
% version 09_10_2014:
% - some cosmetic changes added
% - note! Found an error in rand_distrib: some distributions were not
% correct (x=y)
%__________________________________________________________
%%


%% Parameters
 % change default values of parameters to ones provided in params
 script=mfilename; %disp(script)
 params = change_parameters(params);

dt=params.dt;  % Simulation time step
 dt_out=params.dt_out; % data output time step
 t_fin=params.t_fin; % end time of simulation

l0=params.l0;    % cell length
 w0=params.w0;  % cell width

x_goal=params.x_goal; % target ParB coordinate (when segregation is "done"), coordinates range: -0.5*l..0.5*l

totA=params.totA;  % number of ParA
 totB=params.totB;  % number of ParB

R_B=params.R_B;
R_B_link=params.R_B_link;
 R_A=params.R_A;
D_parA=params.D_parA; % D coefficient for DNA-bound ParA
 D_parB=params.D_parB; % D coefficient for ParB (in complex with ParS)
 
sigma_x=params.sigma_x; % Amplitude of chromosome oscillations, long axis (spatial SD of locus movement) 
 sigma_y=params.sigma_y; % Amplitude of chromosome oscillations, short axis (spatial SD of locus movement)

k_hyd=params.k_hyd; % rate of ParB-stimulated hydrolysis, s. Results in ParA dissociation from DNA
 k_db=params.k_db; % rate of ParA binding to DNA, s. Combined with ATP/ADP-exchange and ParA dimerization
 
x_B0=params.x_B0; y_B0=params.y_B0; % initial ParB coordinates
if length(x_B0)==1 && totB>0, x_B0=x_B0*ones(1,totB); y_B0=y_B0*ones(1,totB); end 
x_p=params.x_p; y_p=params.y_p; % coordinates to "park" free parA

dPdx_choice=params.dPdx_choice;
sigma_w=params.sigma_w; %width of normal distribution for ParA loclaization (migh not be used)

%% INTIALIZATION

% geometric parameters
R02=(R_B+R_A)^2;
R_B_link2=R_B_link^2;
l00=l0/2;    % cell boundaries:
w00=w0/2;    % [-l00..l00,-w00..w00]

% setting 3 states of ParA
indxA=1:totA;
 osciA=indxA;
 freeA=[];
 for bb=1:totB
   boundA{bb}=[];
   boundA_new{bb}=[];
 end
 boundA_all=cat(2,boundA{:});
 n_boundA=zeros(1,totB); % number of ParA bound to i-th ParB
% characteristic reaction times
tau_hyd=1/k_hyd;  % ParA ATP-hydrolysis/dissociation
 tau_free=1/k_db; % ParA dimerization/binding to DNA
% spring constants in kT units from Amplitude of chromosome oscillations
k_spA_x=1/(sigma_x^2); 
 k_spA_y=1/(sigma_y^2);

% useful multipliers 
DdtA11=2*dt*D_parA; DdtA12=dt*dt*D_parA; DdtA22=2*dt*dt*dt*D_parA/3; 
 DdtB11=2*dt*D_parB; DdtB12=dt*dt*D_parB; DdtB22=2*dt*dt*dt*D_parB/3;
MuA=[0,0]; SigmaA=[DdtA11,DdtA12; DdtA12,DdtA22];
 MuB=[0,0]; SigmaB=[DdtB11,DdtB12; DdtB12,DdtB22];
k_xA= D_parA*k_spA_x*dt; k_xA2=k_xA^2/2; % k_xA2 is used in 2-nd order aproximation
 k_yA= D_parA*k_spA_y*dt; k_yA2=k_yA^2/2;
k_xB= D_parB*k_spA_x*dt; k_xB2=k_xB^2/2; % k_xB2 is used in 2-nd order aproximation
 k_yB= D_parB*k_spA_y*dt; k_yB2=k_yB^2/2; 

% initial coordinates  and 
% equlibrium coordinates for parA 
x_B=x_B0;
 y_B=y_B0;
rand_XY=get_rand_XY(totA,dPdx_choice); 
x_A0=rand_XY(1,:);  %l0*rand(1,totA)-l00;
 y_A0=rand_XY(2,:); %w0*rand(1,totA)-w00;
x_A=x_A0; %+l0*(rand(size(x_A0))-0.5); perturbation for test purposes
 y_A=y_A0; %+w0*(rand(size(x_A0))-0.5); perturbation for test purposes

% data collectors
ParA=111*ones(ceil(t_fin/dt_out)+1,2*totA);         % ParA coordinates
 ParA_0=111*ones(ceil(t_fin/dt_out)+1,2*totA);      % ParA equlibrium positions
 ParA_state=222*ones(ceil(t_fin/dt_out)+1,totA);    % state of A
ParB=111*ones(ceil(t_fin/dt_out)+1,2*totB);         % ParB coordinates
 ParB_dxymax=111*ones(ceil(t_fin/dt_out)+1,2*totB); % ParB biggest steps, for controls
% initiating ParA/ParB coordinates and states 
ParA(1,1:2:end)=x_A;
 ParA(1,2:2:end)=y_A;
ParA_0(1,1:2:end)=x_A0;
 ParA_0(1,2:2:end)=y_A0;
ParA_state(1,:)=zeros(1,totA); % state of A
ParB(1,1:2:end)=x_B;
 ParB(1,2:2:end)=y_B; 
ParB_dxymax(1,1:2:end)=zeros(1,totB);
 ParB_dxymax(1,2:2:end)=zeros(1,totB);
deltaXb=zeros(1,totB);   % sum of deviations of ParA bound to i-th ParB 
 deltaYb=zeros(1,totB);


%% SIMULATION

ii0=1; % counter for data saving
%ff=0; hh=0; % counters for ParA hydrolysis and free parA rebinding events
%ddA=0; ddB=0; % counters for ParA and parB diffusion steps
T_free=[]; 
for bb=1:totB, T_hyd{bb}=[]; end

%pause on

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
for tt1=dt_out:dt_out:t_fin % "save cycle", save data at each step
  
  % pool of random reaction times   
  Tpool_hyd=exprnd(tau_hyd,[1,2*totA*ceil(dt_out/tau_hyd)]);  
   Tpool_free=exprnd(tau_free,[1,2*totA*ceil(dt_out/max([tau_hyd,tau_free]))]);
  % pool of random steps
  %Dpool_A=DdtA*randn(2,1*totA*ceil(dt_out/dt));
   %Dpool_B=DdtB*randn(1,2*ceil(dt_out/dt));
  Dpool_A=mvnrnd(MuA,SigmaA,2*totA*ceil(dt_out/dt))'; 
   Dpool_B=mvnrnd(MuB,SigmaB,2*ceil(dt_out/dt)*totB)';
 
  ddA=0; ddB=0; % counters for ParA and parB diffusion steps
  ff=0; hh=0; % counters for ParA hydrolysis and free parA rebinding events
  max_dx=zeros(1,totB); max_dy=max_dx; %collectors of bigest steps
  
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for tt2=dt:dt:dt_out % "silent" cycle, move without saving data
      
    tt=tt1+tt2; % current time
        
    % = update current values
    % - update ReactionTimes and to_parB-Distances 
    T_free=T_free-dt;
     T_hyd=cellfun(@(x) x-dt,T_hyd, 'UniformOutput', false);
     T_hyd_all=cat(2,T_hyd{:});
     %T_hyd=T_hyd-dt;
    %tot_boundA_new=0;
    %boundA_old=boundA;
    % boundA_all_old=boundA_all;
    % freeA_old=freeA;
    % osciA_old=osciA;
    
    for bb=randperm(totB)   % bb=1:totB
      dR2=(x_A(osciA)-x_B(bb)).^2+(y_A(osciA)-y_B(bb)).^2;
      boundA_new{bb}=osciA(dR2<=R02);
      osciA=osciA(dR2>R02);
      %tot_boundA_new=tot_boundA_new+length(boundA_new{bb});
    end
    
    % - update states of ParA
    %boundA_new=osciA(dR2<=R02); 
     freeA_new=boundA_all(T_hyd_all<=0); 
     osciA_new=freeA(T_free<=0);
    
    freeA=[freeA(T_free>0),freeA_new]; 
     osciA=[osciA,osciA_new];
     for bb=1:totB
       %[~,ii1,ii2]=intersect(boundA{bb},boundA_all);  
       boundA{bb}=[boundA{bb}(T_hyd{bb}>0),boundA_new{bb}]; 
       n_boundA(bb)=length(boundA{bb});
     end
     boundA_all=cat(2,boundA{:});
    
    if ~isempty(osciA_new) 
      rand_XY=get_rand_XY(length(osciA_new),dPdx_choice); 
      x_A(osciA_new)=rand_XY(1,:); x_A0(osciA_new)=rand_XY(1,:); 
       y_A(osciA_new)=rand_XY(2,:); y_A0(osciA_new)=rand_XY(2,:); 
    end 
    
    if ~(length(osciA)+length(boundA_all)+length(freeA)==length(unique([osciA,freeA,boundA_all])))
         qq=1;
     end   
    
    % - remove/add Reaction Times according changes in parA states
    for bb=1:totB
      T_hyd{bb}=T_hyd{bb}(T_hyd{bb}>0);
      d_indx=length(boundA_new{bb});
     if d_indx>0,  
       T_hyd{bb}=[T_hyd{bb}, Tpool_hyd(hh+1:hh+d_indx)]; 
       hh=hh+d_indx; 
    end
    end 
    T_free=T_free(T_free>0) ;
    d_indx=length(freeA_new);
     if d_indx>0,  
       T_free=[T_free, Tpool_free(ff+1:ff+d_indx)];  
       ff=ff+d_indx; 
     end  
          
      % = make moves
     % - ParB moves, if interparticle distance> max_length - apply string
     xB_old=x_B; yB_old=y_B;
     for bb=randperm(totB) %1:totB
       deltaXb=sum((x_A(boundA{bb})-x_A0(boundA{bb})));
        deltaYb=sum((y_A(boundA{bb})-y_A0(boundA{bb})));
     
       xB_new=x_B(bb)+Dpool_B(1,ddB+1)- k_xB*deltaXb-k_xB2*n_boundA(bb).*deltaXb-k_xB*n_boundA(bb).*Dpool_B(2,ddB+1);
        yB_new=y_B(bb)+Dpool_B(1,ddB+2)- k_yB*deltaYb-k_yB2*n_boundA(bb).*deltaYb-k_yB*n_boundA(bb).*Dpool_B(2,ddB+2);
       x_B(bb)=(sign(l00-abs(xB_new))+1).*xB_new/2+(sign(abs(xB_new)-l00)+1).*(l0-abs(xB_new)).*sign(xB_new)/2;
        y_B(bb)=(sign(w00-abs(yB_new))+1).*yB_new/2+(sign(abs(yB_new)-w00)+1).*(w0-abs(yB_new)).*sign(yB_new)/2;
       
       d_xB1=x_B(bb)-x_B(max(1,bb-1)); d_xB2=x_B(bb)-x_B(min(totB,bb+1));
        d_yB1=y_B(bb)-y_B(max(1,bb-1)); d_yB2=y_B(bb)-y_B(min(totB,bb+1));
       d_r2B1=d_xB1.^2+d_yB1.^2; d_r2B2=d_xB2.^2+d_yB2.^2;
       indx=find([d_r2B1,d_r2B2]-R_B_link2>0); %s_indx=sum(indx);
       if ~isempty(indx)
         for ii=randperm(length(indx))
           jj=indx(ii);  
           switch jj
             case 1  
               x_B(bb)=x_B(bb-1)+(x_B(bb)-x_B(bb-1))*R_B_link/sqrt(d_r2B1);
                y_B(bb)=y_B(bb-1)+(y_B(bb)-y_B(bb-1))*R_B_link/sqrt(d_r2B1);
             case 2
               x_B(bb)=x_B(bb+1)+(x_B(bb)-x_B(bb+1))*R_B_link/sqrt(d_r2B2);
                y_B(bb)=y_B(bb+1)+(y_B(bb)-y_B(bb+1))*R_B_link/sqrt(d_r2B2);
           end
         end
       end
        
       ddB=ddB+2;
     end     
     
     % - ParA moves
     % "park" free parA
     x_A(freeA)=x_p*ones(1,length(freeA));
      y_A(freeA)=x_A(freeA);
     % move bound ParA
     for bb=1:totB
       x_A(boundA{bb})=x_A(boundA{bb})+(x_B(bb)-xB_old(bb)); 
        y_A(boundA{bb})=y_A(boundA{bb})+(y_B(bb)-yB_old(bb));
     end 
     % move oscillating ParA
     d_indx=length(osciA);
     if d_indx>0
      deltaXo=(x_A(osciA)-x_A0(osciA));
       deltaYo=(y_A(osciA)-y_A0(osciA));
      x_new=x_A(osciA)+Dpool_A(1,ddA+1:ddA+d_indx)- k_xA*deltaXo-k_xA2*deltaXo-k_xA*Dpool_A(2,ddA+1:ddA+d_indx);
       y_new=y_A(osciA)+Dpool_A(1,ddA+d_indx+1:ddA+2*d_indx)- k_yA*deltaYo-k_yA2*deltaYo-k_yA*Dpool_A(2,ddA+d_indx+1:ddA+2*d_indx);
      x_A(osciA)=(sign(l00-abs(x_new))+1).*x_new/2+(sign(abs(x_new)-l00)+1).*(l0-abs(x_new)).*sign(x_new)/2;
       y_A(osciA)=(sign(w00-abs(y_new))+1).*y_new/2+(sign(abs(y_new)-w00)+1).*(w0-abs(y_new)).*sign(y_new)/2;
     end
     ddA=ddA+2*d_indx;
     
     % keep biggest ParB steps for control
     max_dx=max(abs(x_B-xB_old),max_dx);
      max_dy=max(abs(y_B-yB_old),max_dy);
       
  end % "silent" cycle
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   % Current Data collection
   ii0=ii0+1;
   ParA(ii0,1:2:end)=x_A;
    ParA(ii0,2:2:end)=y_A;
   ParA_0(ii0,1:2:end)=x_A0;
    ParA_0(ii0,2:2:end)=y_A0;
   ParB(ii0,1:2:end)=x_B;
    ParB(ii0,2:2:end)=y_B;
   ParB_dxymax(ii0,1:2:end)=max_dx;
    ParB_dxymax(ii0,2:2:end)=max_dy;
   %ParB(ii0,:)=[x_B,y_B,max_dx,max_dy];
   ParA_state(ii0,osciA)=0*ones(1,length(osciA));     % 0= ocsillating ParA
    ParA_state(ii0,freeA)=-1*ones(1,length(freeA));   % -1= free ParA
   for bb=1:totB                                      % ii= ParA bound to ii-th ParB  
     ParA_state(ii0,boundA{bb})=bb*ones(1,length(boundA{bb}));   
   end  
   
   params.rand_fun=rand_fun;
   %disp(num2str([tt, length(osciA),length(boundA_all),length(freeA), length(osciA)+length(boundA_all)+length(freeA),length(unique([osciA,freeA,boundA_all]))],3))
   
   % check if goal is reached by at least one ParB then return
   max_x_B=max(x_B);
   if max_x_B>x_goal
     return  
   end

end % "save" cycle
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

%% NESTED FUNCTIONS
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function rand_XY=get_rand_XY(nn,choice)
% returns n random XY pairs ([x1,x2,x3...; y1,y2,y3...]) using dPdx according to choice
  uni_rand=rand(2,nn); % uniform random numbers
      
  switch choice
    case 1 % X: linear gradient from 0 to 1 between the most right ParB location and new pole
           % Y: uniform distribution
      %x_B_max=max(x_B);
      x_B_max=max(x_B);
      rand_XY(1,:)=sqrt(uni_rand(1,:))*(l00-x_B_max)+x_B_max;
      rand_XY(2,:)=w0*uni_rand(2,:)-w00;
       rand_fun='X: linear gradient from 0 to 1; Y:uniform distribution';
    case 2 % X: linear gradient from 0 to 1 between current ParB location and new pole
           % Y: norm distribution with sigma_w
      x_B_max=max(x_B);
      rand_0=sigma_w*randn(1,2*nn*max([1,sigma_w/w0])); 
      rand_1=rand_0(abs(rand_0)<w00);
      rand_XY(1,:)=sqrt(uni_rand(1,:))*(l00-x_B_max)+x_B_max;
      rand_XY(2,:)=rand_1(1:nn);
       rand_fun='X: linear gradient from 0 to 1; Y:norm distribution';
    case 3 % X: uniform distribution between poles
           % Y: uniform distribution
      %x_B_max=max(x_B);     
      rand_XY(1,:)=l0*uni_rand(1,:)-l00;
      rand_XY(2,:)=w0*uni_rand(2,:)-w00;
       rand_fun='X: uniform distribution; Y:uniform distribution';
    case 4 % X: uniform distribution betwen poles
           % Y: norm distribution with sigma_w
      %x_B_max=max(x_B);
      rand_0=sigma_w*randn(1,2*nn*max([1,sigma_w/w0]));
      rand_1=rand_0(abs(rand_0)<w00);
      rand_XY(1,:)=l0*uni_rand(1,:)-l00;
      rand_XY(2,:)=rand_1(1:nn);
      rand_fun='X: uniform distribution; Y: norm distribution';
    case 5 % X: special distribution for explanation movies
           % Y: uniform distribution with wo/2 
       % totA=7 randomly around [0.25,0.4,0.55,0.7,0.7,0.85,0.85] relative positions
      c=0.05; % max deviations from target point
      rand_XY(1,:)=c*l0*(uni_rand(1,:)-0.5)+([0.25,0.4,0.55,0.68,0.77,0.85,0.9]*l0-l00);
      rand_XY(2,:)=w00*uni_rand(2,:)-w00/2;
       rand_fun='X: special for movies; Y:uniform distribution with w0/2';
    case 7 % X: linear gradient from y0 to y1 between current ParB location and new pole
           % Y: uniform distribution
      c=2; % ratio y1/y0
      x_B_max=max(x_B);
      rand_XY(1,:)=((sqrt(1+(c-1)*(c+1)*uni_rand(1,:))-1)/(c-1))*(l00-x_B_max)+x_B_max;
      rand_XY(2,:)=w0*uni_rand(2,:)-w00;
       rand_fun='X: linear gradient from from y0 to y1; Y:uniform distribution';
    otherwise
      disp('unknown choice... sorry')  
      return  
  end

end

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function params_out = change_parameters(params_in)
% declare default values of parameters and change them to values specified in params_in 
  % declare default parameters values
  params0.dt=0.001;  % Simulation time step
   params0.dt_out=1; % data output time step
   params0.t_fin=2000; % end time of simulation

  l_0=2.3; w_0=0.5; % cell length and width
   params0.l0=l_0; %l_0=2.5;    % cell length
   params0.w0=w_0; %w_0=0.5;  % cell width

  params0.x_goal=1.11; %0.45*l_0; % target ParB coordinate (when segregation is "done"), coordinates range: -0.5*l..0.5*l

  params0.totA=90;  % number of ParA
   params0.totB=1;  % number of ParB

  params0.R_B=0.05;
   params0.R_A=0.002;
  params0.R_B_link=10*3; 
  params0.D_parA=0.01; % D coefficient for DNA-bound ParA
   params0.D_parB=0.0001; % D coefficient for ParB (in complex with ParS)
 
  params0.sigma_x=2*0.06; % Amplitude of chromosome oscillations, long axis (spatial SD of locus movement) 
   params0.sigma_y=2*0.04; % Amplitude of chromosome oscillations, short axis (spatial SD of locus movement)

  params0.k_hyd=0.03; % rate of ParB-stimulated hydrolysis, s. Results in ParA dissociation from DNA
   params0.k_db=0.1; % rate of ParA binding to DNA, s. Combined with ATP/ADP-exchange and ParA dimerization
 
  params0.x_B0=-0.59;  params0.y_B0=0; %params0.x_B0=-0.45*l_0;  params0.y_B0=0; % initial ParB coordinates
   params0.x_p=-10;  params0.y_p=-10; % coordinates to "park" free parA

  params0.dPdx_choice=1; % choice of distrbution for ParA binidng to DNA
  params0.sigma_w=0.15; %width of normal distribution for ParA loclaization (migh not be used)
  
  params0.script=script; %'diffusing_multiple_chained_ParB_oscillating_ParA_v4';
  params0.model='DNA-relay, multi chained ParB';
   params0.version='09.10.2014';
  
  % check what is provided params_in and act change those in params
  fldnms_in = fieldnames(params_in);
   fldnms_0 = fieldnames(params0);
  params_out=params0; 
  for kk = 1:length(fldnms_in)
    %if isempty(nparam.(fni{k})), continue, end
    
    not_found=1;
    for jj = 1:length(fldnms_0)
      if (strcmpi(fldnms_in{kk},fldnms_0{jj}))
         params_out.(fldnms_0{jj}) = params_in.(fldnms_in{kk});
         not_found=0;
      end
    end
    if not_found
        disp(['Warning: parameter ',fldnms_in{kk},' not found in the list of parameters, check spelling...'])
    end
  end
  if params_out.dt_out<params_out.dt
     params_out.dt_out=params_out.dt; 
  end
  
end

%% END of ENDS

end
