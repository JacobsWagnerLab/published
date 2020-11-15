function [X_A,X_R,TT,params_out]=simulate_Rebinding_byReplisome(X_A0,params)
% script to simulate redistribution of DNA-binding protein ("A") that is 
%  (1) irreversibly binds to DNA
%  (2) binds to DNA with uniform probability, i.e. everywhere with the same probability
%  (3) dissociates from the DNA upon encountering Replisome
% NOTE: (2) also implies that unobound A has enough time to distribute uniformly
%      over the DNA before binding occurs; and that there is no saturation
%      of binding
% Considering 
%  (1) only one replisome moving from ori to ter, 
%     presented as L1=[0,1] interval in genomic coordimates
%  (2) there is doubled DNA amount behind the replisome
%  (3) number of A doubles over cell cycle by constitutive synthesis
%     time of newly sythesized A is stochasctis and distributed uniformally
%     over cell cycle time
%  (4) replisome is a "point" (),i.e. size is=0
%
% INPUT:
%   X_A0 - initial coordinates of A
%   params - parameters for simulation with fields:
%   .dt - time step of simulation
%   .n_A - number of A to simulate
%   .ccT -  cell cycle time, G1, S, and G2/M, min
%   .dist - initial distribution of A: 'uniform', 'linear'
%
% OUTPUT:
%   X_A0 - coordinates of A vs time, 1 row is a snap shot of coordinates at a given time 
%          organized as a matrix (2*n_A by T/dt+1, T is cell cycle preiod) 
%          if A has not appeared then A coordinate=-999 (note that coordinates are not sorted)
%   X_R - coordinates of the replisome, -999 if before or after replication
%   TT - array of the time points at wich X_A0 and X_R were evalauted
%   params_out - parameters used for the simulations inside this function,
%                plus date/version/script_name info, just for record/checks
% Convention:
%  second copy of replicated chromosome has coordinates with "-" sign, i.e., 
%  L2 changes from 0 to [-1,0] interval during replication
%   
% Ivan Surovtsev
% 2016.09.08 

%=====================
% PARAMETERS

version='2016.09.08';
script=mfilename;

% ccT=[40,50,10];  % cell cycle time, G1, S, and G2/M, min
% n_A=100000;         % number of A to simulate
% dt=1;         % time step for simualtions, min 
  
%=====================
% SETTING PARAMS
ccT=params.ccT;  % cell cycle time, G1, S, and G2/M, min
 dt=params.dt;         % time step for simualtions, min 
 dist=params.dist; % what to use for initial distribution of A 
 
%=====================
% INITIALIZATION

t_g1=ccT(1); % breaking down cycle into individual times
 t_s=ccT(2);
 t_m=ccT(3);
 tau=sum(ccT); % whole cell cycle
v=1/t_s;      % velocity of replisome along the

% setting X_A0
if isempty(X_A0)
  if isfield(params,'n_A')
    n_A=params.n_A;         % number of A to simulate
  else  
    n_A=1000;
  end
  switch dist
    case 'uniform'  
      X_A0=rand(1,n_A);
    case 'linear' % descending linear gradient from 1 to 0
      X_A0=sqrt(rand(1,n_A))*(0-1)+1;  
  end
else
  n_A=length(X_A0);  
  X_A0=X_A0;  
end

TT=-1*ones(1+tau/dt,1);
X_R=-999*ones(1+tau/dt,1);
%  X_A0=rand(1,n_A);

ii=1;
TT(ii)=0;
X_A=[]; X_A(ii,:)=X_A0;
TA=tau*rand(n_A,1); %times of new A appearance 
 TA=sort(TA);
 
%=====================
% ------MAIN-----
% TIME EVOLUTION

% Before Replication
for  tt=dt:dt:t_g1
  
  ii=ii+1;  
   TT(ii)=tt;
   X_R(ii)=-777;
  
  % just adding newly sythesized A 
  X_Aii=X_A(ii-1,:);
  n_Anew=sum(( (TA>tt-dt) & (TA<=tt) ));
  if n_Anew>0
    X_Anew=rand(1,n_Anew);
    X_A=[X_A,-999*ones(ii-1,n_Anew)]; % adding -999 to "non-born" A
    X_A(ii,:)=[X_Aii,X_Anew];
  else
    X_A(ii,:)=X_Aii;   
  end
  
end

% Replication
X_R(ii)=0;
for  tt=t_g1+dt:dt:t_g1+t_s
  
  ii=ii+1;  
   TT(ii)=tt;
   X_R(ii)=X_R(ii-1)+v*dt;
      
    % just adding newly sythesized A 
  X_Aii=X_A(ii-1,:);
  n_Anew=sum(( (TA>tt-dt) & (TA<=tt) ));
  if n_Anew>0
    X_Anew=rand(1,n_Anew)*(1+X_R(ii))-X_R(ii);
    X_A=[X_A,-999*ones(ii-1,n_Anew)]; % adding -999 to "non-born" A
    X_Aii=[X_Aii,X_Anew];
    X_A(ii,:)=X_Aii;
  else
    X_A(ii,:)=X_Aii;   
  end
  
  % redistributing what was hit by the replisome
  ind=( X_Aii>X_R(ii-1) & X_Aii<=X_R(ii));
  n_Ared=sum(ind);
  if n_Ared>0
    X_Ared=rand(1,n_Ared)*(1+X_R(ii))-X_R(ii);
    X_Aii(ind)=X_Ared;
    X_A(ii,:)=X_Aii;  
  end
    
end

L2=-X_R(ii);
% After Replication
for  tt=t_g1+t_s+dt:dt:tau
  
  ii=ii+1;  
   TT(ii)=tt;
   X_R(ii)=-666;
  
  % just adding newly sythesized A 
  X_Aii=X_A(ii-1,:);
  n_Anew=sum(( (TA>tt-dt) & (TA<=tt) ));
  if n_Anew>0
    X_Anew=rand(1,n_Anew)*(1-L2)+L2;
    X_A=[X_A,-999*ones(ii-1,n_Anew)]; % adding -999 to "non-born" A
    X_A(ii,:)=[X_Aii,X_Anew];
  else
    X_A(ii,:)=X_Aii;   
  end
  
end

%=====================
% OUTPUT
params_out.dt=dt;
 params_out.ccT=ccT;
 params_out.n_A=n_A;
params.script=script;
 params.version=version;
 params.date=date;


end % END of the ends
