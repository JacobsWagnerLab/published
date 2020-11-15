function [Output,params]=generations_simulations_9(n_gen,X_exp)

% to simulate division cycles by generations
%
%
% changed internal variable to predifined ones...
% with correlation between Mother and Daughter growth rates
% 4: with Divsion ratio variations
% 4x: to include different growth by ST and SW cells
% 5: to include option with experimental dPd(DeltaL)
% 6ec: to play with parameters starting from values for E.coli
%      added a shortcut for negative value in distributions: take abs
%      for every distrib...
% 7: to inroduce possible correlation between dltaL and alpha
% 8: to add alpha-correlation between generations
% 9: saves used parameters in param structure

% n_gen = number of generatiobs to run

version='2014.07.10'; % version ID
 script=mfilename; disp(script)

n_cells_00=10e3; % number of cells at t=0;
 n_max=10e3; % maximum number of cells to follow
L0=5.5;   % average cell length at birth, only for initial distribution
 cv_L0=4.5/5.5; % CV of length at birth
 % in parameters below convention is: x(1) - for ST cell, x(2) - for SW cell 
DivR0=0.556;     % average division ratio (defined as ST/(ST+SW), i.e. ~0.54 for Caulobacter)
 cv_DivR=0.075; % CV of division ratio fluctuations 
DeltaL0(1)=1.81; %2.45; %2.8; %3.23; %3.14; % average length increment per cycle, ST cells
  DeltaL0(2)=DeltaL0(1)*DivR0/(1-DivR0) ; % average length increment per cycle, SW cells
 cv_DeltaL0(1)=0.125; %0.12; %0.15; %184; % CV of length increment per cycle
  cv_DeltaL0(2)=cv_DeltaL0(1); %cv_DeltaL0(1); % CV of length increment per cycle
%T0(1)=41; % average doubling time in min, ST cells
% T0(2)=T0(1); % average doubling time in min, SW cells
alpha0(1)=0.0255; % average growth rate in 1/min, ST cells
  alpha0(2)=0.0255; % average growth rate in 1/min, SW cells
 cv_alpha(1)=0.; %0.071; % CV of growth rates (stalk cells)
  cv_alpha(2)=cv_alpha(1); %cv_alpha(1); % CV of growth rates (swarmer cells)
 corr_alpha(1)=0.0; %0.44; %0.16; %0.16; % growth rate correlation between Mother and Daughter cells
  corr_alpha(2)=corr_alpha(1); % growth rate correlation between Mother and Daughter cells
%rng('shuffle','combRecursive'); % to randomly reset RND-generator and change generator method from default
rng('shuffle'); % to randomly reset RND-generator and change generator method from default
X_max=20; %maximal possible X from P(X) (i.e. maximum elongation), "fininte substitute for Infininty"
corr_choice=0; % to choose if alpha and DeltaL are correlated (0 - use various_distrib, 1 - various_distrib_2)
 corr_adl=0.0; % correlation between alpha and deltaL
dist_choice=[1,2,1]; % distribution choice, array with order [DivR,DeltaL,alpha]
 dist_choice0=[3,11]; % distribution choice for initial size and CCT distrib [Lb0,tau0]

 % List of available distributions
dist_list{1}{1}='Normal Distribution';   % for variuos_distrib
 dist_list{1}{2}='LogNormal Distribution';
 dist_list{1}{3}='Uniform Distribution';
 dist_list{1}{5}='Normal Distribution correlated with prior values';
 dist_list{1}{6}='Experimental Distribution';
 dist_list{1}{10}='Delta(0)';
 dist_list{1}{11}='Delta(x_mean)';
 dist_list{1}{20}='Special case: fixed number of specified values';
dist_list{2}{2}='Bi-variate: Normal + LogNormal Distributions'; % for variuos_distrib_2
 dist_list{2}{3}='Bi-variate: LogNormal + LogNormal Distributions';
 dist_list{2}{5}='Bi-variate: Normal correlated with prior values + LogNormal Distributions';

% initialization
Output=[]; params=[];
t0=0;
n_cells_0=n_cells_00;
%alpha0(1)=log(1/DivR0)/T0(1); % average specific growth rate, ST cells
% alpha0(2)=log(1/(1-DivR0))/T0(2); % average specific growth rate, SW cells

% Setting initial cells distribution
gen=1; % generation #
%Cells=single((1:n_cells_0))'; % "current cells", their ID number
Cells=(1:n_cells_0)'; % "current cells", their ID number
GGen=gen*ones(n_cells_0,1); % generation # for each cell
 % assuming that all starting cells are ST:
   %Lb=various_distrib(DeltaL0(2),cv_DeltaL0(2)/2,1,n_cells_0); % individual length at birth: LogNormal Distribution 
   %Lb=various_distrib(L0,cv_L0,3,n_cells_0); % individual length at birth: Uniform Distribution 
Lb=various_distrib(L0,cv_L0,dist_choice0(1),n_cells_0); % individual length at birth: Uniform Distribution 
  %DeltaL=various_distrib(DeltaL0(1),cv_DeltaL0(1),2,n_cells_0); % individual increment by the division
  %DeltaL=various_distrib(DeltaL0(1),cv_DeltaL0(1),6,n_cells_0); % individual increment by the division
if dist_choice(2)==5, dc_tmp=2; else dc_tmp=dist_choice(2); end  
switch corr_choice
  case 1
    AlphaDelta_1=various_distrib_2([alpha0(1),DeltaL0(1)],[cv_alpha(1),cv_DeltaL0(1)],corr_adl,dc_tmp,n_cells_0);
     DeltaL=AlphaDelta_1(:,2); alpha=AlphaDelta_1(:,1); % individual deltaL and growth rates
  case 0
    DeltaL=various_distrib(DeltaL0(1),cv_DeltaL0(1),dist_choice(2),n_cells_0); % individual increment by the division  
     alpha=various_distrib(alpha0(1),cv_alpha(1),dc_tmp,n_cells_0); % individual growth rates    
end
 Ld=Lb+DeltaL; % individual cell length at the division
DivR=various_distrib(DivR0,cv_DivR,dist_choice(1),n_cells_0); % individual division ratio at the division
  %DivR=various_distrib(DivR0,cv_DivR,3,n_cells_0); % individual division ratio at the division
  %  alpha=various_distrib(alpha0(1),cv_alpha(1),11,n_cells_0); % individual growth rates
   
TT=log(1+DeltaL./Lb)./alpha;
   %TTb=t0*ones(n_cells_0,1); % indivdiual birth times
   %TTb=-TT.*various_distrib(0.5,1,3,n_cells_0); % indivdiual birth times
 TTb=various_distrib(t0,t0,dist_choice0(2),n_cells_0); % indivdiual birth times
 TTd=TTb+TT; % indivdiual division times
MotherID=0*ones(n_cells_0,1); % ID of Mother cells
 ST_SW=ones(n_cells_0,1); % Marker of ST cells
 
if length(Cells)>n_max
  rand_indx=randperm(length(Cells),n_max);
  %rand_indx=randperm(length(Cells)); 
  ii=sort(rand_indx); %ii=ii(1:n_max);
else
  ii=1:length(Cells);  
end
ii_2(1:2:2*length(ii))=ii; ii_2(2:2:2*length(ii))=ii;
Output=[Output;single([Cells(ii),GGen(ii),TTb(ii),TTd(ii),Lb(ii),Ld(ii),TT(ii),DivR(ii),MotherID(ii),ST_SW(ii),alpha(ii)])];

for gen=2:n_gen
Cells_old=Cells(ii);
 Cells=(max(Cells_old)+1:max(Cells_old)+2*length(Cells_old))';
GGen=gen*ones(length(Cells),1); % generation # for each cell 
% Setting new generation of cells 
  Lb(1:2:2*length(Cells_old),1)=Ld(ii).*DivR(ii); % individual length at birth
   Lb(2:2:2*length(Cells_old),1)=Ld(ii).*(1-DivR(ii)); % individual length at birth
    %DeltaL(1:2*length(Cells_old),1)=various_distrib(DeltaL0,cv_DeltaL0,2,2*length(Cells_old)); % individual increment by the division
    %AlphaDelta_1=various_distrib_2([alpha0(1),DeltaL0(1)],[cv_alpha(1),cv_DeltaL0(1)],corr_adl,3,length(Cells_old));
    % AlphaDelta_2=various_distrib_2([alpha0(2),DeltaL0(2)],[cv_alpha(2),cv_DeltaL0(2)],corr_adl,3,length(Cells_old));
    %alpha(ii_2)
  DivR(1:1:2*length(Cells_old),1)=various_distrib(DivR0,cv_DivR,dist_choice(1),2*length(Cells_old)); % individual division ratio at the division
  alpha_now=alpha(ii); %ST_SW_now=ST_SW(ii);
  switch corr_choice
    case 1
      AlphaDelta_1=various_distrib_2([alpha0(1),DeltaL0(1)],[cv_alpha(1),cv_DeltaL0(1)],corr_adl,dist_choice(2),length(Cells_old),alpha_now,corr_alpha(1));
       AlphaDelta_2=various_distrib_2([alpha0(2),DeltaL0(2)],[cv_alpha(2),cv_DeltaL0(2)],corr_adl,dist_choice(2),length(Cells_old),alpha_now,corr_alpha(2)); 
      DeltaL(1:2:2*length(Cells_old),1)=AlphaDelta_1(:,2); % individual increment by the division
       DeltaL(2:2:2*length(Cells_old),1)=AlphaDelta_2(:,2); % individual increment by the division
      alpha(1:2:2*length(Cells_old),1)=AlphaDelta_1(:,1); % individual growth rates
       alpha(2:2:2*length(Cells_old),1)=AlphaDelta_2(:,1); % individual growth rates
    case 0
     DeltaL(1:2:2*length(Cells_old),1)=various_distrib(DeltaL0(1),cv_DeltaL0(1),dist_choice(2),length(Cells_old)); % individual increment by the division
      DeltaL(2:2:2*length(Cells_old),1)=various_distrib(DeltaL0(2),cv_DeltaL0(2),2,length(Cells_old)); % individual increment by the division       
        %DeltaL(1:2:2*length(Cells_old),1)=various_distrib(DeltaL0(1),cv_DeltaL0(1),6,length(Cells_old)); % individual increment by the division
     % DeltaL(2:2:2*length(Cells_old),1)=various_distrib(DeltaL0(2),cv_DeltaL0(2),6,length(Cells_old)); % individual increment by the division
     %alpha(1:2*length(Cells_old),1)=alpha0*(1+cv_alpha*various_distrib(12*length(Cells_old))); % individual growth rates
     %alpha(1:2*length(Cells_old),1)=various_distrib(alpha0,cv_alpha,5,2*length(Cells_old)); % individual growth rates
     alpha(1:2:2*length(Cells_old),1)=various_distrib(alpha0(1),cv_alpha(1),dist_choice(3),length(Cells_old),alpha_now,corr_alpha(1)); % individual growth rates
      alpha(2:2:2*length(Cells_old),1)=various_distrib(alpha0(2),cv_alpha(2),dist_choice(3),length(Cells_old),alpha_now,corr_alpha(2)); % individual growth rates
  end
     %DivR(1:2*length(Cells_old),1)=various_distrib(DivR0,cv_DivR,1,2*length(Cells_old)); % individual division ratio at the division
       %DivR(1:1:2*length(Cells_old),1)=various_distrib(DivR0,cv_DivR,3,2*length(Cells_old)); % individual division ratio at the division
  MotherID(1:2:2*length(Cells_old),1)=Cells_old; % ID of Mother cells
   MotherID(2:2:2*length(Cells_old),1)=Cells_old; % ID of Mother cells
  ST_SW(1:2:2*length(Cells_old),1)=ones(length(Cells_old),1); % Marker of ST cells   
   ST_SW(2:2:2*length(Cells_old),1)=zeros(length(Cells_old),1); % Marker of ST cells   
  TTb(1:2:2*length(Cells_old),1)=TTd(ii);
   TTb(2:2:2*length(Cells_old),1)=TTd(ii);
  Ld=Lb+DeltaL; % individual cell length at the division
    %DivR=DivR0*ones(length(Cells),1); % individual division ratio at the division
  TT=log(1+DeltaL./Lb)./alpha;
  TTd=TTb+TT;
  
if length(Cells)>n_max
     %rand_cells=randi(length(Cells),n_max,1);
  rand_indx=randperm(length(Cells),n_max);
     %rand_indx=randperm(length(Cells)); 
  ii=sort(rand_indx); %ii=ii(1:n_max);
else
  ii=1:length(Cells);  
end
ii_2(1:2:2*length(ii))=ii; ii_2(2:2:2*length(ii))=ii;

Output=[Output;single([Cells(ii),GGen(ii),TTb(ii),TTd(ii),Lb(ii),Ld(ii),TT(ii),DivR(ii),MotherID(ii),ST_SW(ii),alpha(ii)])];

  
end

%jj=10;
%n_cells=length(Output(:,1)>0);
%Output=Output(1:n_cells,:);

% assembling parameters output.
params.script=script;
 params.version=version;
params.dist_choice=dist_choice;
 params.distrib_DivR=dist_list{1}{dist_choice(1)};
 params.distrib_deltaL=[dist_list{corr_choice+1}{dist_choice(2)}];
 params.distrib_alpha=dist_list{corr_choice+1}{dist_choice(3)};
params.initial_dist_choice=dist_choice0;
 params.distrib_Lb=dist_list{1}{dist_choice0(1)};
 params.distrib_CCT=dist_list{1}{dist_choice0(2)};
params.n_cells_0=n_cells_0; % number of cells at t=0;
 params.n_max=n_max; % maximum number of cells to follow
params.L0=L0;   % average cell length at birth, only for initial distribution
 params.cv_L0=cv_L0; % CV of length at birth
params.DivR0=DivR0;     % average division ratio (defined as ST/(ST+SW), i.e. ~0.54 for Caulobacter)
 params.cv_DivR=cv_DivR; % CV of division ratio fluctuations 
params.DeltaL0=DeltaL0; %2.45; %2.8; %3.23; %3.14; % average length increment per cycle, ST cells
 params.cv_DeltaL0=cv_DeltaL0; %0.12; %0.15; %184; % CV of length increment per cycle
params.alpha0=alpha0; % average growth rate in 1/min, ST cells
 params.cv_alpha=cv_alpha; %0.071; % CV of growth rates (stalk cells)
 params.corr_alpha=corr_alpha; %0.16; %0.16; % growth rate correlation between Mother and Daughter cells
params.CorrAlphaDeltaL=[corr_choice,corr_adl]; % correlation between alpha and deltaL 
 %params.CorrAlphaDeltaL=corr_adl; % correlation between alpha and deltaL
params.MaxExpDist=X_max; %maximal possible X from P(X) (i.e. maximum elongation), "fininte substitute for Infininty"


disp([num2str(n_gen),' generations simulated. Output Format (by columns):'])
disp(['   1    | ','    2      | ',' 3  | ',' 4  | ',' 5  | ',' 6  | ','7 | ','    8    | ','    9     | ',' 10 | ',' 11  ']);
disp(['Cell ID | ','generation | ','T_b | ','T_d | ','L_b | ','L_d | ','T | ','DivRatio | ','Mother ID | ','ST? | ','alpha']);


% DONE


%% NESTED FUNCTIONS
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

% choice of distributions
function distrib_out=various_distrib(x_mean,x_cv,cchoice,nn,x_0,x_corr)
  % x0 is previous sampling for correlated variables (i.e. growth rates)  
  % to choose distribution
  switch cchoice
    case 1 % Normal Distribution
      distrib = single(randn(nn,1));
      distrib_out=x_mean*(1+x_cv*distrib); 
    case 2 % LogNormal Distribution
      mu=log(x_mean^2/sqrt(x_mean^2+(x_cv*x_mean)^2)); sigma=sqrt(log(1+(x_cv*x_mean)^2/x_mean^2));   
      distrib_out = single(random('Lognormal',mu,sigma,nn,1));
    case 3 % Uniform Distribution between [L0-SD,L0+SD];
      distrib = single(rand(nn,1));
      distrib_out = x_mean*(1-x_cv)+2*x_cv*x_mean*distrib;
    case 5 % Normal distribution correlated with input x_0
        % GR_1 = mean (Growth rates) [GR_mother, GR_STdaughter, GR_SWdaughter]
        % GR_2 = std (Growth rates) [GR_mother, GR_STdaughter, GR_SWdaughter,]
        % GR_3 = [corr(GR_mother, GR_STdaughter), corr(GR_mother, GR_SWdaughter)]
      distrib_out = (x_mean+x_corr*(x_0-x_mean))+x_cv*x_mean*sqrt(1-x_corr^2)*randn(nn,1);
    case 6 % using experimental distribution
      XXe=sort(X_exp); XXe=[0;XXe;X_max];
       YYe=ones(length(XXe)-1,1); YYe=[0;YYe]; YYe=cumsum(YYe)/sum(YYe);
      RRand=rand(nn,1);
      distrib_out= interp1(YYe,XXe,RRand,'pchip');  
    case 10 % delta(0) Distribution (i.e. all the same 0 value)
      distrib_out = single(zeros(nn,1));
    case 11 % delta(x_mean) Distribution (i.e. all the same x_mean value)
      distrib_out=x_mean*single(ones(nn,1));  
    case 20 % Special case: fixed number of cells withgiven L0
      distrib_out = single([1;2;3;4;5;6;7;8;9;10]);
    otherwise % do delta-distribution(1) (i.e. same value)
      distrib_out=x_mean*single(ones(nn,1));
  end    
  
  distrib_out=abs(distrib_out);
   
end

function distrib_out=various_distrib_2(x_mean,x_cv,xx_corr,cchoice,nn,x_0,x0_corr)
  % returns 2 sets of correlated variables  
  % x0 is previous sampling for correlated variables (i.e. growth rates)  
  % choice = to choose distribution
  switch cchoice % Currently only Gaussian capula's
    case 2 % Normal + LogNormal Distributions
      mu_1=x_mean(1); sigma_1=x_cv(1)*x_mean(1);  
      mu_2=log(x_mean(2)^2/sqrt(x_mean(2)^2+(x_cv(2)*x_mean(2))^2)); 
       sigma_2=sqrt(log(1+(x_cv(2)*x_mean(2))^2/x_mean(2)^2));
      Mu=[mu_1,mu_2];
       Sigma=[sigma_1^2,xx_corr*sigma_1*sigma_2; xx_corr*sigma_1*sigma_2,sigma_2^2];
      ZZ=mvnrnd(Mu,Sigma,nn);
      UU(:,1)=normcdf(ZZ(:,1),mu_1,sigma_1);
       UU(:,2)=normcdf(ZZ(:,2),mu_2,sigma_2);
      XX(:,1)=ZZ(:,1);
       XX(:,2)=logninv(UU(:,2),mu_2,sigma_2); 
      distrib_out = single(XX);
    case 3 % LogNormal + LogNormal Distributions
      mu_1=log(x_mean(1)^2/sqrt(x_mean(1)^2+(x_cv(1)*x_mean(1))^2)); 
       sigma_1=sqrt(log(1+(x_cv(1)*x_mean(1))^2/x_mean(1)^2));
      mu_2=log(x_mean(2)^2/sqrt(x_mean(2)^2+(x_cv(2)*x_mean(2))^2)); 
       sigma_2=sqrt(log(1+(x_cv(2)*x_mean(2))^2/x_mean(2)^2));
      Mu=[mu_1,mu_2];
       Sigma=[sigma_1^2,xx_corr*sigma_1*sigma_2; xx_corr*sigma_1*sigma_2,sigma_2^2];
      ZZ=mvnrnd(Mu,Sigma,nn);
      UU(:,1)=normcdf(ZZ(:,1),mu_1,sigma_1);
       UU(:,2)=normcdf(ZZ(:,2),mu_2,sigma_2);
      XX(:,1)=logninv(UU(:,1),mu_1,sigma_1); 
       XX(:,2)=logninv(UU(:,2),mu_2,sigma_2); 
      distrib_out = single(XX);  
    case 5 % Normal correlated with previous values input x_0 and LogNormal
      mu_1=x_mean(1); sigma_1=x_cv(1)*x_mean(1);  
      mu_2=log(x_mean(2)^2/sqrt(x_mean(2)^2+(x_cv(2)*x_mean(2))^2)); 
       sigma_2=sqrt(log(1+(x_cv(2)*x_mean(2))^2/x_mean(2)^2));
      Mu=[mu_1,mu_2];
       Sigma=[sigma_1^2,xx_corr*sigma_1*sigma_2; xx_corr*sigma_1*sigma_2,sigma_2^2];
      ZZ=mvnrnd(Mu,Sigma,nn);
      UU(:,1)=normcdf(ZZ(:,1),mu_1,sigma_1);
       UU(:,2)=normcdf(ZZ(:,2),mu_2,sigma_2);
      XX(:,1)=norminv(UU(:,1),mu_1+x0_corr*(x_0-mu_1),sigma_1*sqrt(1-x0_corr^2));
       XX(:,2)=logninv(UU(:,2),mu_2,sigma_2); 
      distrib_out = single(XX);
    otherwise % do uncorrelated normal
      Mu=[x_mean(1),x_mean(2)];
       Sigma=[x_cv(1)*x_mean(1),0; 0,x_cv(2)*x_mean(2)];
      distrib_out = single(mvnrnd(Mu,Sigma,nn));
      disp('Warning: Requested unknown bivariate distribution, used bi-variate uncorrelated Normal')
  end    
  
  %distrib_out=abs(distrib_out);
   
end


%% END of ENDS
end
