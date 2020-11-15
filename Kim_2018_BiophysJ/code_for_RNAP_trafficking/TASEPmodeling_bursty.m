function lociList = TASEPmodeling_bursty(pauseProfile,promoter,totalLoci,fOn,avgSpeed,pauseSite, pauseHeight)
%{
-About-
This function runs TASEP modeling of transcription only (no translation or
mRNA degradation).
This function iterates by regular for loop (not by parfor), hence slower
than _par version.
Bursty promoter case

-Inputs-
pauseProfile: name of the pause profile (added for consistency with _par version
function)
promoter: name of the promoter (added for consistency with _par version
function)
totalLoci: (integer) number of iterations (this will be the size of
output)(see "Note" under Example)
fOn: The fraction of time the promoter is ON. Value is between 0 and 1
avgSpeed: average speed of RNAP 
pauseSite: location of the pause (nt)
pauseHeight: duration of the pause (sec)

-varargin-

-Outputs-
lociList: RNAP trajectories, time points of promoter ON/OFF, time points of RNAP loading events,
times and locations of RNAP-RNAP collisions

-Example-
TASEPmodeling_bursty('flat','P1',10,0.25,30,0,0)


Note: to generate figures of RNAP trajectories, it is okay to use low value
for totalLoci
to obtain headway distribution, use high value for totalLoci

-Supplementary-

-Keywords-
TASEP modeling, gene expression

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
geneLength = 3075; % nt 
polWidth = 35;  % nt
N_totalloci = totalLoci;
simTime = 0:1:40*60; % simulation time in seconds

% Bursty initiation parameters
kOn = 0.007; % sec-1 from So et al. Nat Genetics (2011)
kLoading = 0.45; % sec-1 from So et al. Nat Genetics (2011)
% fOn this is input of this function 
kOff = (1-fOn)*kOn/fOn;

%% Define RNAP dwelltime profile 
avgDwelltime = dx/avgSpeed; %sec per nucleotide
if strcmp(pauseProfile, 'flat')
    specificDwelltime = avgDwelltime * ones(geneLength/dx,1); % no pause, flat profile
    runCondition = strcat(pauseProfile,'-NO-');
elseif strcmp(pauseProfile, 'OnepauseAbs')
    specificDwelltime = avgDwelltime * ones(geneLength/dx,1); % no pause, flat profile
    specificDwelltime(pauseSite) = pauseHeight;
    runCondition = strcat(pauseProfile,sprintf('%01.0f',pauseSite),'x',sprintf('%01.0f',pauseHeight),'-NO-');
elseif strcmp(pauseProfile, 'MultipauseAbs')
    specificDwelltime = avgDwelltime * ones(geneLength/dx,1); % no pause, flat profile
    specificDwelltime(pauseSite) = pauseHeight;
    runCondition = strcat(pauseProfile,sprintf('%01.0f',length(pauseSite)),'x',sprintf('%01.0f',pauseHeight),'-NO-');
end;


for lociID = 1:N_totalloci %each loci is independent
    
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
    % [2] exitTime = stepping time at each nt of an RNAP
    % [3] polStatus = status of RNAP at each nt
    % 0 not yet loaded or elongating (0.5 for collision); 3 loading inhibited;
    loadTime = [];
    cycleIndex = [];
    exitTime = -inf(geneLength/dx+polWidth/dx,1); 
    polStatus = zeros(geneLength/dx+polWidth/dx,1);
    
    j = 0;
    for i = 1:Numonoffcycle
        t = onTime(i) + exprnd(1/kLoading);
        while t < offTime(i)
            j = j + 1; %loaded RNAP on a DNA
            loadTime(j) = t;
            cycleIndex(j) = i; %cycle number for each RNAP
            if j == 1
                % for the first RNAP, determine stepping through entire length 
                % because there is none ahead.
                exitTime(1:geneLength/dx,j) = t + cumsum(exprnd(specificDwelltime));
            else
                exitTime(1:geneLength/dx,j) = zeros(geneLength/dx,1);
            end;
            polStatus(1:geneLength/dx,j) = zeros(geneLength/dx,1);
            t = t + exprnd(1/kLoading);
        end;
    end;
    
    %---Elongation/ RNAP translocation
    % Iterate over RNAP (p = polID)-> iterate over x = 1:geneLength
    for p = 2:length(loadTime)
        % check which pol is ahead on the template
        ahead = 0;
        for i = p:-1:2
            if polStatus(1,i-1)<1
                ahead = i-1; 
                break;
            end;
        end;
        
        % check if the loading attempt (loadTime) is hindered by RNAP ahead
        % loadTime = arrivalTime at x = 1 = (exitTime at x = 0)
        if loadTime(p) <= exitTime(polWidth/dx,ahead) 
            %loading inhibited
            polStatus(:,p) = 3*ones(size(polStatus,1),1); % Loading failed
            exitTime(1:geneLength/dx,p) = -inf(geneLength/dx,1); % reset (no loading)
        else
            % translocate the RNAP 
            for x = 1:geneLength
                % arrival time at x = exitTime at x-1
                if x == 1
                    arrival = loadTime(p);
                else
                    arrival = exitTime(x-1,p);
                end;

                % exitTime at x = arrival time at x + dwell time at x (=shift)
                shift = exprnd(specificDwelltime(x));
                exitTime(x,p) = arrival + shift;

                % check for collision with RNAP ahead
                % collision -> trailing RNAP polStatus = 0.5; leading RNAP
                % polStatus = 0.2
                % collision -> wait until leading RNAP steps
                if arrival + shift <=  exitTime(x+polWidth/dx,ahead) 
                    polStatus(x,p) = polStatus(x,p)+0.5; 
                    polStatus(x+polWidth/dx,ahead) = polStatus(x+polWidth/dx,ahead)+ 0.2;
                    exitTime(x,p) = exitTime(x+polWidth/dx,ahead) + shift;
                end;
            end;
        end;
    end;
    
    lociList(lociID).ontime = onTime;
    lociList(lociID).offtime = offTime;
    lociList(lociID).loadtime = loadTime;
    lociList(lociID).cycleindex = cycleIndex;
    lociList(lociID).polexittime = exitTime;
    lociList(lociID).polstatus = polStatus;
end;
    
    
