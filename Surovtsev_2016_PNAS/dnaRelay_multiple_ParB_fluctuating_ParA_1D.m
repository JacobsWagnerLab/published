function [ParB,ParA,tt,params,ParA_state,ParA_0] = dnaRelay_multiple_ParB_fluctuating_ParA_1D(params)
% To simulate 1D dynamics of ParA dimers and ParB-rich cargo within DNA-relay model
% same as dnaRelay_multiple_ParB_fluctuating_ParA but 1D version
%   See Surovtsev et al. PNAS 2016 113 E7266 for details
%    http://www.pnas.org/content/113/46/E7268.full
%   This is a version used for the paper 
%
% This version of the DNA-relay model considers
%   (1) Multiple ParB-rich cargos, they can overlap without interaction
%   (2) ParB-cargos difuses freely with D_parB diffusion constant, if they are not engaged with DNA-bound ParA dimers.
%       Diffusion is implicit - modelled as a random binding with uniform probability (or choice of the distribution) throught the cell 
%   (3) ParB-cargos engaged with DNA-bound ParA dimers experience an
%        elastic force from each individual ParA-ParB interaction
%   (4) ParA dimers cycle through 3 states: 
%        (i)free ParA (called "free_A" in the code) diffusing freely though cytoplasm, 
%        (ii)DNA-bound ParA ("bound_A") fluctuates (with D_parA diffusion coefficeint) 
%            around an equilibrium position x0 and experience elastic forc F=k_sp(x_A-x0)  
%        (iii)ParB/DNA-bound ParA ("osci_A") moves together with the cargo
%   (5) Transitions (i)->(ii), (iii)->(i) are modelled as stochastic processes with exponential distribution of times
%        with mean times tau_free=1/k_db, tau_hyd=1/k_hyd, repsectively.
%        k_db and k_hyd are rate constants for "DNA-binding" and
%        "ATP-hydrolysis" by ParA
%   (6) Transition (ii)->(iii) occurs when cargo and ParA dimer, considered as spheers with radius R_B and R_A, repsectively, overlap
%   (7) If cargo(s) reaches target x-coordinate (x_goal), simulation stops (e.g., to model attchment of the cargo) 
%        To simulate cargos without attachement (e.g., plasmids) put x_goal>l0 (cell length)
%
% Simulations organized as a nested loops with "silent" internal loop that run simulations without saving coordinates, 
%   and external loop in which corrdinates are saved at each iteration.
% All simulations are within a rectangular container with reflective boundaries  
%
% INPUT:
%   params - parameters for simulation with multiple fields, for example:
%     params.dt - Simulation time step
%     params.dt_out - data output time step
%     params.t_fin - end time of simulation
%   Full list of parameters can be found in a nested function "change_parameters" below
%   Defaults values will be used for parameters value/fields not provided in input
%
% OUTPUT:
%   ParB - coordinates of ParB cargos vs time, i-th row is a snap shot of coordinates of all cargos (X1,X2...) at a given time,(i-1)*dt_out 
%   ParA - coordinates of ParA dimers vs time, i-th row is a snap shot of coordinates of all ParA dimers (X1,X2,...) at a given time,(i-1)*dt_out
%          all diffusing ParA ('free_A') have a fixed value, as their diffusion is not simulated explicitly and they do not interact with ParB untill they re-bind to DNA
%   tt - simulation time
%   params - parameters used for simulation, also contain versio/date info
%   ParA_state - state of the ParA dimer, i-th row is a list of states (S1,S2,S3,...) at a given time,(i-1)*dt_out, using following convention:
%           0=DNA-bound ParA, 
%           -1=free_ParA
%           ii=Cargo/DNA-bound ParA
%   ParA_0 - coordinates of equlibrium positions of ParA dimers, i-th row is a snap-shot of all ParA dimers (X1,X2,...) at a given time,(i-1)*dt_out
%   
% Ivan Surovtsev
% 2016.10.30 
%__________________________________________________________
%%


%% Parameters
%% Parameters
 % First, define all required parameters of simulations 
 % If values provided in the input "params" - change default values to ones provided in params
 % Note: parameters units should be self-consistent. For example, use seconds and micrometers as all time and length units  
 script=mfilename; % save code name for the output 
 params = change_parameters(params); % get all parameters values

% Split 'params' into individual parameters
 % simulation time parameters
dt=params.dt;  % Simulation time step
dt_out=params.dt_out; % data output time step
t_fin=params.t_fin; % end time of simulation
 % simulation box parameters
l0=params.l0;    % cell length
x_goal=params.x_goal; % target ParB coordinate (when segregation is "done"), coordinates range: -0.5*l..0.5*l

 % ParA and ParB-related parameters
totA=params.totA;  % number of ParA
totB=params.totB;  % number of ParB
R_B=params.R_B; % radius of ParB-rich cargo, considering a disk
R_A=params.R_A; % radius of ParA dimer, considering a disk
D_parB=params.D_parB; % Diffusion coefficient for ParB (in complex with ParS)
D_parA=params.D_parA; % Diffusion coefficient for DNA-bound ParA
 % elastisity of the chromosome
sigma_x=params.sigma_x; % Amplitude of chromosome oscillations, long axis (spatial SD of locus movement) 
 % rate of biochemical steps
k_hyd=params.k_hyd; % rate of ParB-stimulated hydrolysis. Results in ParA dissociation from DNA
k_db=params.k_db; % rate of ParA binding to DNA. Effectively combinec two steps: ATP/ADP-exchange and ParA dimerization

 % initial coordinates of the cargo(s) 
x_B0=params.x_B0; % get initial ParB coordinates
% if only one pair of coordinates provided by there are multiple cargos - replicate the values
if length(x_B0)==1 && totB>0, x_B0=x_B0*ones(1,totB); end 

% choice of the random function to determine coordinates of ParA when it rebinds to the DNA
% see "get_rand_XY" nested function at the end of the code for details
dPdx_choice=params.dPdx_choice;


% values assigned to the free ParA in output, relevant only for analysis of the output data 
x_p=params.x_p; % coordinates to "park" free parA

%% INTIALIZATION
rng('shuffle'); % to randomly reset RND-generator and change generator method from default

 % geometric parameters
R02=(R_B+R_A)^2; % "overlap" distance between ParA and parB 
% coordinates of the cell boundaries, with the origin in the center of the cell, [-l00..l00,-w00..w00]
l00=l0/2; % long-axis
 
 % setting 3 states of ParA
indxA=1:totA; % just an indexing array 1,2,3,...
osciA=indxA; % initially all ParA are DNA-bound
freeA=[];    % no freely diffusing ParA at t=0
% and no ParA bound to the cargo
% boundA is organized as as a cell array, with i-th cell containing the list of ParA bound to the i-th ParB-rich cargo 
for bb=1:totB
   boundA{bb}=[];
   boundA_new{bb}=[];
end
boundA_all=cat(2,boundA{:}); % to have a list of all ParB-bound ParAs
n_boundA=zeros(1,totB); % number of ParA bound to i-th ParB, none initially
 % characteristic reaction times, to use in generation of random numbers in stochastic reactions
tau_hyd=1/k_hyd;  % ParA ATP-hydrolysis/dissociation
tau_free=1/k_db; % ParA dimerization/binding to DNA
 % spring constants in kT units from Amplitude of chromosome oscillations, to use in force calculations
k_spA_x=1/(sigma_x^2); % along long-axis

% Below we define useful multipliers, they will be used at the calculations of displacements according to our BD-scheme
% Check SI in the paper for details, scheme is based on Branka AC and Heyes DM (1998) Phys.Rev.E 58: 2611–261
 % Some constants for random numbers generation 
DdtA11=2*dt*D_parA; DdtA12=dt*dt*D_parA; DdtA22=2*dt*dt*dt*D_parA/3; % for ParA
DdtB11=2*dt*D_parB; DdtB12=dt*dt*D_parB; DdtB22=2*dt*dt*dt*D_parB/3; % for ParB-cargo
 % Correlation matrices for generating correlated random numbers with constants from the aboove 
MuA=[0,0]; SigmaA=[DdtA11,DdtA12; DdtA12,DdtA22]; % for ParA
MuB=[0,0]; SigmaB=[DdtB11,DdtB12; DdtB12,DdtB22]; % for ParB-cargo
 % 1-st and 2-nd order constants, respectively, for force-related terms in displacement calculations 
k_xA= D_parA*k_spA_x*dt; k_xA2=k_xA^2/2; % for ParA, along long axis
k_xB= D_parB*k_spA_x*dt; k_xB2=k_xB^2/2; % for ParB, along long axis

% setting initial coordinates and equlibrium positions
% initial coordinates for ParB-cargo 
x_B=x_B0;
% equilibrium coordinates for ParA, randomly distributed
rand_X=get_rand_X(totA,dPdx_choice); 
x_A0=rand_X(1,:);  %l0*rand(1,totA)-l00;
% initial coordinates for ParA. We will start from equilibrium positions. Note: one can add random (reasonable) shift, it does not influence results 
x_A=x_A0; %+l0*(rand(size(x_A0))-0.5); perturbation for test purposes

% making data collectors for the output
ParA=111*ones(ceil(t_fin/dt_out)+1,1*totA);         % ParA coordinates
ParA_0=111*ones(ceil(t_fin/dt_out)+1,1*totA);      % ParA equlibrium positions
ParA_state=222*ones(ceil(t_fin/dt_out)+1,totA);    % state of ParA
ParB=111*ones(ceil(t_fin/dt_out)+1,1*totB);         % ParB coordinates
% saving initial values - ParA/ParB coordinates and states - at t=0;
ParA(1,1:1:end)=x_A;
ParA_0(1,1:1:end)=x_A0;
ParA_state(1,:)=zeros(1,totA); % state of A
ParB(1,1:1:end)=x_B;


%% SIMULATION

ii0=1; % counter for data saving

% initializing list of stochactis times for considered reactions(transitions)
T_free=[]; % list of stochactis times for free ParA -> DNA-bound ParA transition
for bb=1:totB, T_hyd{bb}=[]; end % list of stochactis times for ParB-cargo/DNA-bound ParA -> free ParA transition, by individual cargos

%pause on

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
for tt1=dt_out:dt_out:t_fin % "save cycle", save data at each step
  
  % Get pools of random times for reaction events
  % if k=0 take a time longer than the total simulation time
  % otherwise draw from exponential distribution with appropriate time constant
  %
  % for ATP-hydrolysis (ParB-cargo/DNA-bound ParA -> free ParA transition)
  switch sign(k_hyd) 
    case 0
      Tpool_hyd=10*t_fin*ones([1,2*totA]);  
    case 1
      Tpool_hyd=exprnd(tau_hyd,[1,2*totA*ceil(dt_out/tau_hyd)]);  
  end
  % for DNA re-binding (free ParA -> DNA-bound ParA transition)
  switch sign(k_db)
   case 0
     Tpool_free=10*t_fin*ones([1,2*totA]);  
   case 1
     Tpool_free=exprnd(tau_free,[1,2*totA*ceil(dt_out/max([tau_hyd,tau_free]))]);
  end
 
  % Get pools of random steps, i.e., Brownian dynamics steps
  Dpool_A=mvnrnd(MuA,SigmaA,1*totA*ceil(dt_out/dt))'; 
  Dpool_B=mvnrnd(MuB,SigmaB,1*ceil(dt_out/dt)*totB)';
 
  % counters for drawing numers from the pools 
  ddA=0; ddB=0; % counters for ParA and ParB diffusion steps
  ff=0; hh=0; % counters for ATP-hydrolysis and DNA-rebinding events
  
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for tt2=dt:dt:dt_out % "silent" cycle, move without saving data
      
    tt=tt1+tt2; % current time
        
    %~~~~~ UPDATE current values:
    % - update lists of Reaction Times since dt time past 
    T_free=T_free-dt; % free ParA -> DNA-bound ParA transition
    T_hyd=cellfun(@(x) x-dt,T_hyd, 'UniformOutput', false); % ParB-cargo/DNA-bound ParA -> free ParA transition
    T_hyd_all=cat(2,T_hyd{:}); % to have combined list of all ParB-cargo/DNA-bound ParA -> free ParA transition 
    
    % - update lists of to-ParB-cargo distances, going cargo by cargo but in random order 
    for bb=randperm(totB)   
      dR2=(x_A(osciA)-x_B(bb)).^2; % distance from all DNA-bound ParAs to a given cargo
      boundA_new{bb}=osciA(dR2<=R02); % if distance smaller than overlap, then transition DNA-bound ParA -> ParB-cargo/DNA-bound ParA
      osciA=osciA(dR2>R02); % only DNA-bound ParA non-overlapping with any cargo will keep their state
      %tot_boundA_new=tot_boundA_new+length(boundA_new{bb});
    end
    
    % - update states of ParA
    freeA_new=boundA_all(T_hyd_all<=0); % make new transitions: ParB-cargo/DNA-bound ParA -> free ParA
    osciA_new=freeA(T_free<=0); % make new transitions: free ParA -> DNA-bound ParA
    freeA=[freeA(T_free>0),freeA_new]; % combine "old" with "new" free ParA
    osciA=[osciA,osciA_new]; % combine "old" with "new" DNA-bound ParA
    for bb=1:totB  % combine "old" with "new" ParB-cargo/DNA-bound ParA, by individual cargo
      boundA{bb}=[boundA{bb}(T_hyd{bb}>0),boundA_new{bb}]; 
      n_boundA(bb)=length(boundA{bb});
    end
    boundA_all=cat(2,boundA{:}); % to have a list of all ParB-bound ParAs
    
     % - set coordinates for ParA newly bound to DNA 
    if ~isempty(osciA_new) 
      rand_X=get_rand_X(length(osciA_new),dPdx_choice); 
       % equlibrium coordinates
      x_A0(osciA_new)=rand_X(1,:);
       % current coordinates, will start from the equlibrium as they will quickly change
      x_A(osciA_new)=rand_X(1,:); 
    end 
       
    % - remove/add Reaction Times according changes in ParA states
     % ParB-cargo/DNA-bound ParA -> free ParA transition
    for bb=1:totB
      T_hyd{bb}=T_hyd{bb}(T_hyd{bb}>0);
      d_indx=length(boundA_new{bb}); % this is how many DNA-bound ParA are newly bound to the cargo
      if d_indx>0, % if there are some DNA-bound ParA are newly bound add their ATP-hydrolysis reaction times 
        T_hyd{bb}=[T_hyd{bb}, Tpool_hyd(hh+1:hh+d_indx)]; 
        hh=hh+d_indx; 
      end
    end
     % free ParA -> DNA-bound ParA transition
    T_free=T_free(T_free>0) ;
    d_indx=length(freeA_new);  % this is how many new free ParA were generated
     if d_indx>0,  % if there are some new free ParA add their DNA re-binding reaction times 
       T_free=[T_free, Tpool_free(ff+1:ff+d_indx)];  
       ff=ff+d_indx; 
     end  
    %~~~~~ END of "Update current values"
          
    % ~~~~~ Make BD moves
     % - ParB moves
     xB_old=x_B; 
     for bb=randperm(totB) 
       % distances for elastic force calculations 
       deltaXb=sum((x_A(boundA{bb})-x_A0(boundA{bb})));
       % new coordinates according to eq.S1 in the SI of the paper 
       xB_new=x_B(bb)+Dpool_B(1,ddB+1)- k_xB*deltaXb-k_xB2*n_boundA(bb).*deltaXb-k_xB*n_boundA(bb).*Dpool_B(2,ddB+1);
       % apply reflecting boundaries, if necessary
       x_B(bb)=(sign(l00-abs(xB_new))+1).*xB_new/2+(sign(abs(xB_new)-l00)+1).*(l0-abs(xB_new)).*sign(xB_new)/2;
       % shift the counter
       ddB=ddB+1;
     end     
     
     % - ParA moves
     % = "park" free parA
     x_A(freeA)=x_p*ones(1,length(freeA));
     % = move Parb-cargo/DNA-bound ParA using cargo displacements
     for bb=1:totB
       x_A(boundA{bb})=x_A(boundA{bb})+(x_B(bb)-xB_old(bb)); 
     end 
     % = move DNA-bound ParA
     d_indx=length(osciA);
     if d_indx>0
      % distances for elastic force calculations   
      deltaXo=(x_A(osciA)-x_A0(osciA));
      % new coordinates according to eq.S1 in the SI of the paper  
      x_new=x_A(osciA)+Dpool_A(1,ddA+1:ddA+d_indx)- k_xA*deltaXo-k_xA2*deltaXo-k_xA*Dpool_A(2,ddA+1:ddA+d_indx);
      % apply reflecting boundaries, if necessary
      x_A(osciA)=(sign(l00-abs(x_new))+1).*x_new/2+(sign(abs(x_new)-l00)+1).*(l0-abs(x_new)).*sign(x_new)/2;
     end
     % shift the counter 
     ddA=ddA+1*d_indx;
       
  end % "silent" cycle
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   % Current Data collection
   ii0=ii0+1;
   ParA(ii0,1:1:end)=x_A;
   ParA_0(ii0,1:1:end)=x_A0;
   ParB(ii0,1:1:end)=x_B;
   ParA_state(ii0,osciA)=0*ones(1,length(osciA));     % 0= ocsillating ParA
   ParA_state(ii0,freeA)=-1*ones(1,length(freeA));   % -1= free ParA
   for bb=1:totB                                      % ii= ParA bound to ii-th ParB  
     ParA_state(ii0,boundA{bb})=bb*ones(1,length(boundA{bb}));   
   end  
   
   params.rand_fun=rand_fun;% to save info on random function used
   
   % check if "finish line" is reached by at least one ParB then return
   max_x_B=max(x_B);
   if max_x_B>x_goal
     return  
   end

end % "save" cycle
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

%% NESTED FUNCTIONS
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function rand_X=get_rand_X(nn,choice)
% returns n random X sets ([x1,x2,x3...; ]) using dPdx according to choice
  uni_rand=rand(1,nn); % uniform random numbers
      
switch choice
    case 1 % X: linear gradient from 0 to 1 between the most right ParB location and new pole
      x_B_max=max(x_B);
      rand_X(1,:)=sqrt(uni_rand(1,:))*(l00-x_B_max)+x_B_max;
      rand_fun='X: linear gradient from 0 to 1';
    case 2 % X: linear gradient from 0 to 1 between current ParB location and new pole
      x_B_max=max(x_B);
      rand_X(1,:)=sqrt(uni_rand(1,:))*(l00-x_B_max)+x_B_max;
       rand_fun='X: linear gradient from 0 to 1';
    case 3 % X: uniform distribution between poles 
      rand_X(1,:)=l0*uni_rand(1,:)-l00;
       rand_fun='X: uniform distribution';
    case 4 % X: uniform distribution betwen poles
      rand_X(1,:)=l0*uni_rand(1,:)-l00;
      rand_fun='X: uniform distribution';
    case 21 % special case for initial coordinates: 
            % all with uniform distribution within the boundaries
            % except one pair has the same cooridinates
      xx=l0*(rand(1,nn-1)-0.5); xx=[xx,xx(end)];
      rand_fun='X : uniform distribution, one pair has the same [X]';   
       
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

  params0.x_goal=1.11; % target ParB coordinate (when segregation is "done"), set>l0 for no atachment case

  params0.totA=90;  % number of ParA
  params0.totB=1;  % number of ParB

  params0.R_B=0.05; % radius of ParB-cago disk
  params0.R_A=0.002; % radius of ParA dimer

  params0.D_parA=0.01; % D coefficient for DNA-bound ParA
  params0.D_parB=0.001; % D coefficient for ParB (in complex with ParS)
 
  params0.sigma_x=2*0.06; % Amplitude of chromosome oscillations, long axis (spatial SD of locus movement) 
  params0.sigma_y=2*0.04; % Amplitude of chromosome oscillations, short axis (spatial SD of locus movement)

  params0.k_hyd=0.03; % rate of ParB-stimulated hydrolysis, s. Results in ParA dissociation from DNA
  params0.k_db=0.1; % rate of ParA binding to DNA, s. Combined with ATP/ADP-exchange and ParA dimerization
 
  params0.x_B0=0.0;  % initial ParB-cargo coordinates
  params0.x_p=-100;  % coordinates to "park" free parA (should be outside of the simulation box)

  params0.dPdx_choice=3; % choice of distrbution for ParA binidng to DNA, 3=uniform destribution
  params0.sigma_w=0.15; %width of normal distribution for ParA loclaization (migh not be used)
  
  params0.dim=1; % Dimension of simulations;
  params0.script=script; % to keep track of what code was used;
  params0.model='DNA-relay, 1D multi ParB'; % model name
  params0.version='2016.08.09';
  
   % check what is provided as input parameters, i.e., in 'params_in', and substitute those in default values
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
