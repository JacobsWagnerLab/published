function ans = TASEPmodeling_noEl_bursty_par(promoter,parRun,fOn,mRNALL,proteinLL)
%{
-About-
This function runs TASEP modeling of transcription, translation and mRNA degradation based on parfor function
Bursty promoter case
No elongation! 

-Inputs-
promoter: name of the promoter (this will be included in the output file name)
parRun: integer, indicating # of files in a series (only to name output file name) 
fOn: the fraction of time the promoter is ON. Value is between 0 and 1
mRNALL: mRNA lifetime;
proteinLL: protein lifetime. If 0, do not consider protein degradation.

-varargin-

-Outputs-
output: matlab data file containing RNAP loading number per template,
headway between RNAPs at initiation, 
mRNA number, protein number, mRNA lifetime   

-Example-
TASEPmodeling_noEL_bursty_par('P1',1,0.25,90,0)
   
-Supplementary-

-Keywords-
TASEP modeling, gene expression, parfor, no elongation, elongation free
model

-Dependencies-
masterscript: see masterscript for how to use this function

-References-

-Author-
Sangjin Kim, 2017 September 30
%}

%--------------------------------------------------------------------------

%% Parameters
% All changeable parameters are listed here. Feel free to modify depending
% on the purpose
dx = 1; % lattice site = 1 nt
N_totalloci = 100; %eventually 1000
simTime = 0:1:40*60; % simulation time in seconds

% Input parameters for ribosome initiation
kRiboLoading = 0.2; %(1/sec); 

% Bursty initiation parameters
kOn = 0.007; % sec-1 from So et al. Nat Genetics (2011)
kLoading = 0.45; % sec-1 from So et al. Nat Genetics (2011)
% fOn = 0.25; %this is input of this function 
kOff = (1-fOn)*kOn/fOn;

% Parameters for analysis
sampleWindow = [15*60,30*60]; % simulation time window when the analysis happened
loadingSample = zeros(1,N_totalloci); % # of loading per loci during sampling

% mRNA life-time
binLT = 0:5:1000;
lifeTimeHist1 = zeros(length(binLT),N_totalloci);
lifeTimeAvg = zeros(1,N_totalloci);

%for FISH
probe1 = 1;
fishTime = 0:60:30*60; %max(simTime);
fishSignal1 = zeros(length(fishTime),N_totalloci); % #mRNA per loci

% for steady-state (SS) protein amount
proteinSS =  zeros(length(fishTime),N_totalloci); % for #protein per loci

%for distribution of RNAP interval (t-headway)
binT = 0:5:200;
tDiffStartHist = zeros(length(binT),N_totalloci);

parallelobject = parpool(12);
tic 
parfor lociID = 1:N_totalloci %each loci is independent
    
    %---Promoter ON/OFF
    % Initialize for each loci 
    onTime = 0; %onTime of each cycle
    offTime = 0; %offTime of each cycle

    j = 0;
    % Flip a coin to decide if promoter is ON or OFF at the start of simulation 
    firstOn = rand <= fOn;
    if firstOn % it was middle of ON at t = min(simTime)
        onTime(1) = min(simTime) - exprnd(1/kOff)*rand;
    else % it was middle of OFF at min(simTime)
        onTime(1) = min(simTime) + exprnd(1/kOn)*rand;
    end;
    offTime(1) = onTime(1) + exprnd(1/kOff);
    
    while offTime(j+1) < max(simTime)
        j = j+1;
        onTime(j+1) = offTime(j) + exprnd(1/kOn);
        offTime(j+1) = onTime(j+1) + exprnd(1/kOff); 
    end;
    Numonoffcycle = j+1; %or length(onTime)
    
    %---RNAP loading attempts
    % For each RNAP, there are...
    % [1] loading time = time point of loading "attempt"
    loadTime = [];
        
    j = 0;
    for i = 1:Numonoffcycle
        t = onTime(i) + exprnd(1/kLoading);
        while t < offTime(i)
            j = j + 1; %loaded RNAP on a DNA
            loadTime(j) = t;
            t = t + exprnd(1/kLoading);
        end;
    end;
    
    %---mRNA degradation and FISH
    lifeTime = exprnd(mRNALL,1,length(loadTime));
    if isempty(loadTime)
        lifeTime = [];
    end;
    decayTime = loadTime + lifeTime;
    
    %---Translation
    % Initialize parameters
    RiboloadTime = [];
    h = 0;
    for polID = 1:length(loadTime)
        ribot = loadTime(polID);
        ribot = ribot + exprnd(1/kRiboLoading);
        while ribot < decayTime(polID)
            h = h+ 1; 
            RiboloadTime(h) = ribot;
            ribot = ribot + exprnd(1/kRiboLoading);
        end;          
    end; 
    
    %---Ribo Analysis
    % Check proteins accumulated per DNA during t = 1200-1800 (10 min)
    % Check proteins at each time point
    if isempty(RiboloadTime)
        proteinLoci2(lociID) = 0;
        proteinSS(:,lociID) = 0;
    else
        proteinLoci2(lociID) = length(find(RiboloadTime>1200 & RiboloadTime<=1800));
        proteindecayTime = RiboloadTime + exprnd(proteinLL,1,length(RiboloadTime));
        % for protein SS (FISH)
        for tt = 1:length(fishTime)
            made1 = length(find(RiboloadTime < fishTime(tt)));
            decayed1 = length(find(proteindecayTime <fishTime(tt)));
            loc = zeros(length(fishTime),1); loc(tt) = 1;
            proteinSS(:,lociID) = proteinSS(:,lociID) + (made1-decayed1)*loc;
        end
    end;
    
    %-- mRNA number per simulation
    for tt = 1:length(fishTime)
       made1 =  length(find(loadTime<fishTime(tt)))*1;
       decayed1 = length(find(decayTime<fishTime(tt)))*1;

       loc = zeros(length(fishTime),1); loc(tt) = 1;
       fishSignal1(:,lociID) = fishSignal1(:,lociID) + (made1-decayed1)/length(probe1)*loc;
    end;
    
    %---Analysis during sampleWindow time
    i = 0; 
    lastPol = length(loadTime);
    lifeTimeSite = [0];
    tStart = []; tEnd = [];
    if isempty(lastPol)
        loadingSample(lociID) = 0;
    else
        for polID = 1: lastPol
            if loadTime(polID)>=sampleWindow(1) && loadTime(polID)<sampleWindow(2)     
                i = i+1; % traj index in a lociID
                % number of successful RNAP initiation per loci
                loadingSample(lociID) = loadingSample(lociID) + 1;
                % 5'-end mRNA life time
                lifeTimeSite(i,1) = lifeTime(1,polID);
                % For the interval between RNAPs upon loading
                tStart(i) = loadTime(polID); % or exitTime(1,polID)
            end;
        end;
        lifeTimeHist1(:,lociID) = hist(lifeTimeSite(:,1),binLT)*100/size(lifeTimeSite,1)';
        lifeTimeAvg(:,lociID) = mean(lifeTimeSite(:,1));
        tDiffStartHist(:,lociID) = hist(diff(tStart),binT)*100/max((length(tStart)-1),1);
    end;
end;
delete(parallelobject);
toc

fileName = strcat('NoEL-',promoter,'-',sprintf('%01.0f',parRun),'par.mat');
save(fileName,'loadingSample','fishSignal1',...
    'lifeTimeHist1','lifeTimeAvg',...
    'tDiffStartHist','proteinLoci2','proteinSS',...
    '-v7.3');
