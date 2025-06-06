function lociList = TASEPmodeling_nonbursty(pauseProfile,promoter,totalLoci,kLoading,avgSpeed,pauseSite, pauseDuration, pauseProb)
%{
-About-
This function runs TASEP modeling of transcription only (no translation or
mRNA degradation).
This function iterates by regular for loop (not by parfor), hence slower
than _par version.
Nonbursty promoter case

-Inputs-
pauseProfile: name of the pause profile (added for consistency with _par version
function)
promoter: name of the promoter (added for consistency with _par version
function)
totalLoci: (integer) number of iterations (this will be the size of
output)(see "Note" under Example)
kLoading: RNAP loading rate (transcription initiation rate) in 1/sec
avgSpeed: average speed of RNAP
pauseSite: location of the pause (nt)
pauseDuration: duration of the pause (sec)
pauseProb: probability of pausing (%)

-varargin-

-Outputs-
lociList: RNAP trajectories, RNAP loading time points, time and location of
RNAP-RNAP collisions

-Example-
TASEPmodeling_nonbursty('flat','P1nb',5,1/15,30,0,0,0)

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
% Nonbursty initiation parameters: kLoading (input) in 1/sec

%% Define RNAP dwelltime profile 
avgDwelltime = dx/avgSpeed; %sec per nucleotide
if strcmp(pauseProfile, 'flat')
    specificDwelltime = avgDwelltime * ones(geneLength/dx,1); % no pause, flat profile
    runCondition = strcat(pauseProfile,'-NO-');
elseif strcmp(pauseProfile, 'OnepauseAbs')
    specificDwelltime = avgDwelltime * ones(geneLength/dx,1); % no pause, flat profile
    specificDwelltime(pauseSite) = pauseDuration;
    runCondition = strcat(pauseProfile,sprintf('%01.0f',pauseSite),'x',sprintf('%01.0f',pauseDuration),'xp',sprintf('%01.0f',pauseProb),'-NO-');
elseif strcmp(pauseProfile, 'MultipauseAbs')
    specificDwelltime = avgDwelltime * ones(geneLength/dx,1); % no pause, flat profile
    specificDwelltime(pauseSite) = pauseDuration;
    runCondition = strcat(pauseProfile,sprintf('%01.0f',length(pauseSite)),'x',sprintf('%01.0f',pauseDuration(1)),'x',sprintf('%01.0f',pauseDuration(2)),'xp',sprintf('%01.0f',pauseProb(1)),'-NO-');
end;
 
for lociID = 1:N_totalloci %each loci is independent
    %---RNAP loading attempts
    % For each RNAP, there are...
    % [1] loading time = time point of loading "attempt"
    % [2] exitTime = stepping time at each nt of an RNAP
    % [3] polStatus = status of RNAP at each nt
    % 0 not yet loaded or elongating (0.5 for collision); 3 loading inhibited;
    % This is the section that only differs from bursty mode of run
    loadTime = [];
    pauseRNAP = [];
    exitTime = -inf(geneLength/dx+polWidth/dx,1); 
    polStatus = zeros(geneLength/dx+polWidth/dx,1);
    
    j = 0;
    t = simTime(1) + exprnd(1/kLoading)*rand;
    % because previous loading (exprnd(1/kLoading) before can be before
    % simTime(1)
    while t < simTime(end)
        j = j+ 1; %loaded RNAP on a locus
        loadTime(j) = t;
        pauseRNAP(j) = rand<= (pauseProb/100); % 1= yes pause 0 = no pause
        if j == 1
            if pauseRNAP(j) == 0
                specificDwelltime1 = avgDwelltime * ones(geneLength/dx,1); % no pause, flat profile
            else
                specificDwelltime1 = specificDwelltime;
            end;
            exitTime(1:geneLength/dx,j) = t + cumsum(exprnd(specificDwelltime1));
        else
            exitTime(1:geneLength/dx,j) = zeros(geneLength/dx,1);
        end;
        polStatus(1:geneLength/dx,j) = zeros(geneLength/dx,1);
        t = t + exprnd(1/kLoading);
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
            if pauseRNAP(p) == 0 % no pause
                specificDwelltime1 = avgDwelltime * ones(geneLength/dx,1); % no pause, flat profile
            else
                specificDwelltime1 = specificDwelltime; %original input, with pause
            end;
            for x = 1:geneLength
                % arrival time at x = exitTime at x-1
                if x == 1
                    arrival = loadTime(p);
                else
                    arrival = exitTime(x-1,p);
                end;

                % exitTime at x = arrival time at x + dwell time at x (=shift)
                shift = exprnd(specificDwelltime1(x));
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
    
    lociList(lociID).loadtime = loadTime;
    lociList(lociID).polexittime = exitTime;
    lociList(lociID).polstatus = polStatus;
end;
    
    
