function ansAll = TASEPmodeling_par_analysis_genedosage(pauseProfile,promoter,totalparRun,pauseSite, pauseHeight,f)

%{
-About-
This function analyzes results of TASEPmodeling_(non)bursty_par in
case of the cells in the population contain 1 or 2 genes per cell.
This code combine results from separate parfor runs and analyze the result on the cell population level.
Fraction of cells with 2 gene copies (f) can be changed to see the effect
of gene dosage (e.g. the effect of gene location along the chromosome).

-Inputs-
pauseProfile: name of the pause profile (this will be included in the output file name)
promoter: name of the promoter (this will be included in the output file name)
totalparRun: total par file number (the output files should be numbered as
1, 2, 3, ... totalparRun. If only out output file, use 1)
pauseSite: location of the pause (nt)
pauseHeight: duration of the pause (sec)
f: fraction_cellswith2DNA (between 0 or 1) use 0.578 for lacZ on the E.coli
chromosome

-varargin-

-Outputs-
ansAll: Distribution of mRNAs and
proteins per cell (mean, fano factor, CV2), duration with no new protein
production

-Example-
   
-Supplementary-

-Keywords-
TASEP modeling, gene expression

-Dependencies-
see masterscript

-References-

-Author-
Sangjin Kim, 2017 September 30
%}

%--------------------------------------------------------------------------

%% This code combine results from separate parfor runs and analyze 
%% mRNA and protein distribution for gene dosage problems

fraction_cellswith2DNA = f;

if strcmp(pauseProfile, 'flat')
    runCondition = strcat(pauseProfile,'-NO-');
elseif strcmp(pauseProfile, 'OnepauseAbs')
    runCondition = strcat(pauseProfile,sprintf('%01.0f',pauseSite),'x',sprintf('%01.0f',pauseHeight),'-NO-');
elseif strcmp(pauseProfile, 'MultipauseAbs')
    runCondition = strcat(pauseProfile,specific,'-NO-');
end;

fishSignal1 = []; fishSignal2 = []; fishSignal3 = []; 
proteinLoci2 = []; rEndStampA = [];
% Load all par runs and combine together
for i = 1:totalparRun
    fileName = strcat(runCondition,promoter,'-',sprintf('%01.0f',i),'par.mat'); 
    tmp = load(fileName);
        
    fishSignal1 = [fishSignal1 tmp.fishSignal1];
    fishSignal2 = [fishSignal2 tmp.fishSignal2];
    fishSignal3 = [fishSignal3 tmp.fishSignal3];
    
    rEndStampA = [rEndStampA tmp.rEndStampA];
    proteinLoci2 = [proteinLoci2 tmp.proteinLoci2];
end;

%% mRNA per cell (chromosomal lacZ gene dosage, if f = 0.578)

% We will re-distribute 1000 DNA into cells, such that 60 % cells have 2
% DNA and rest have 1 DNA
cellswith1DNA = round(100*totalparRun/(1+2*fraction_cellswith2DNA/(1-fraction_cellswith2DNA)));
cellswith2DNA = round(cellswith1DNA*fraction_cellswith2DNA/(1-fraction_cellswith2DNA));

mRNAstat1 = []; mRNAstat2 = []; mRNAstat3 = [];
% Repeat for bootstrapping
for shufflei = 1:3000
    clear cellfish*;
    lociIndex = randperm(size(fishSignal1,2));
    for cellnum = 1:cellswith2DNA
        cellfishSignal1(:,cellnum) = fishSignal1(:,lociIndex(cellnum*2-1))+fishSignal1(:,lociIndex(cellnum*2));
        cellfishSignal2(:,cellnum) = fishSignal2(:,lociIndex(cellnum*2-1))+fishSignal2(:,lociIndex(cellnum*2));
        cellfishSignal3(:,cellnum) = fishSignal3(:,lociIndex(cellnum*2-1))+fishSignal3(:,lociIndex(cellnum*2));
    end;
    tmp = cellnum*2;
    for cellnum = (cellswith2DNA+1):(cellswith2DNA+cellswith1DNA)
        cellfishSignal1(:,cellnum) = fishSignal1(:,lociIndex(cellnum-cellswith2DNA+tmp));
        cellfishSignal2(:,cellnum) = fishSignal2(:,lociIndex(cellnum-cellswith2DNA+tmp));
    end;
    
    % Calculate mean, fano factor, CV, CV^2 in each sampling
    % Signal1 = 5'end, Signal2 = 3'end, Signal3 = 72tiling
    ssTime = 30; %steady-state time index of the fish time
    sscellfishSignal1 = []; sscellfishSignal2 = []; sscellfishSignal3 = []; 
    for i = ssTime
        sscellfishSignal1 = [sscellfishSignal1, cellfishSignal1(i,:)]; 
        sscellfishSignal2 = [sscellfishSignal2, cellfishSignal2(i,:)]; 
        sscellfishSignal3 = [sscellfishSignal3, cellfishSignal3(i,:)]; 
    end;

    mRNAstat1(shufflei,:) = [mean(sscellfishSignal1), var(sscellfishSignal1)/mean(sscellfishSignal1),...
        std(sscellfishSignal1)/mean(sscellfishSignal1),(std(sscellfishSignal1)/mean(sscellfishSignal1))^2];
    mRNAstat2(shufflei,:) = [mean(sscellfishSignal2), var(sscellfishSignal2)/mean(sscellfishSignal2),...
        std(sscellfishSignal2)/mean(sscellfishSignal2),(std(sscellfishSignal2)/mean(sscellfishSignal2))^2];
    mRNAstat3(shufflei,:) = [mean(sscellfishSignal3), var(sscellfishSignal3)/mean(sscellfishSignal3),...
        std(sscellfishSignal3)/mean(sscellfishSignal3),(std(sscellfishSignal3)/mean(sscellfishSignal3))^2];
end;

% ans10: mRNA mean, fano, CV, CV&2 and then ste for each
ans10(1,:) = [mean(mRNAstat1,1),std(mRNAstat1,0,1)];
ans10(2,:) = [mean(mRNAstat2,1),std(mRNAstat2,0,1)];
ans10(3,:) = [mean(mRNAstat3,1),std(mRNAstat3,0,1)];
ansAll.mRNAstat = ans10;

% ans11: mRNA distribution based on last sampling result
binN = 0:1:50;
ans11(:,1) = binN';
ans11(:,2) = hist(sscellfishSignal1', binN')*100/length(sscellfishSignal1);
sscellfishSignal1_1DNA = []; sscellfishSignal1_2DNA = [];
for i = ssTime
    sscellfishSignal1_1DNA = [sscellfishSignal1_1DNA, cellfishSignal1(i,(cellswith2DNA+1):(cellswith2DNA+cellswith1DNA))]; 
    sscellfishSignal1_2DNA = [sscellfishSignal1_2DNA, cellfishSignal1(i,1:cellswith2DNA)]; 
end;
%distribution among cells with 2 DNA, normalized by total cells
ans11(:,3) = hist(sscellfishSignal1_2DNA, binN')*100/length(sscellfishSignal1);
%distribution among cells with 1 DNA, normalized by total cells
ans11(:,4) = hist(sscellfishSignal1_1DNA', binN')*100/length(sscellfishSignal1);
% repeat to obtain above two lines to obtain distributions of Signal2 and Signal3 
ansAll.mRNAdist = ans11;

%%  Ribosome per cell (using chromosomal lacZ gene dosage)
for shufflei = 1:3000 
    lociIndex = randperm(size(rEndStampA,2)); 
    clear cellProtein2;
    % cellswith 2 DNA 
    for cellnum = 1:cellswith2DNA
        cellRibostamp(:,cellnum) = rEndStampA(:,lociIndex(cellnum*2-1))+rEndStampA(:,lociIndex(cellnum*2));
        cellProtein2(cellnum) = proteinLoci2(lociIndex(cellnum*2-1))+ proteinLoci2(lociIndex(cellnum*2));
    end;
    % cells with 1 DNA
    tmp = cellnum*2;
    for cellnum = (cellswith2DNA+1):(cellswith2DNA+cellswith1DNA)
        cellRibostamp(:,cellnum) = rEndStampA(:,lociIndex(cellnum-cellswith2DNA+tmp));
        cellProtein2(cellnum) = proteinLoci2(lociIndex(cellnum-cellswith2DNA+tmp));
    end;
    % Calculate mean, fano factor, CV, CV^2 in each sampling
    proteinStat(shufflei,1) = mean(cellProtein2); %mean of protein # per cell
    proteinStat(shufflei,2) = var(cellProtein2)/mean(cellProtein2); %fano
    proteinStat(shufflei,3) = std(cellProtein2)/mean(cellProtein2); %CV
    proteinStat(shufflei,4) = (proteinStat(shufflei,3))^2; %CV^2
end;
% protein number (mean+/-ste; fano +/- ste)
proteinLacZ = [mean(proteinStat(:,1)),std(proteinStat(:,1)),mean(proteinStat(:,2)),std(proteinStat(:,2)),...
    mean(proteinStat(:,3)),std(proteinStat(:,3)),mean(proteinStat(:,4)),std(proteinStat(:,4))];
ansAll.proteinstat = proteinLacZ;

% For distribution
ans14(:,1) = (0:100:2000)'; %binP
ans14(:,2) = hist(cellProtein2',ans14(:,1))*100/length(cellProtein2); % for all cells
ans14(:,3) = hist(cellProtein2(1:cellswith2DNA)',ans14(:,1))*100/length(cellProtein2); % for cells with 2 DNA
ans14(:,4) = hist(cellProtein2(cellswith2DNA+1:end)',ans14(:,1))*100/length(cellProtein2); %for cells with 1 DNA
ansAll.proteindist = ans14;

% Statistics of instantaneous protei production rate from rEndStampA(t,DNA)
% E.g., # of time points with no new protein production
X = cellRibostamp(2:end-1,:);
k = 0; zeroDall = [];
for i = 1:size(X,2)
    % time point with 0 will get 1 time points >0 will get 0
    zeroTs = X(:,i) == 0;
    if isempty(find(X(:,i) ==0)) % no zeros 
        zeroDall = [zeroDall 0];
    else
        h = 0; k = 0; zeroD = 0;
        for j = 1:length(zeroTs)
            if zeroTs(j)==1
                h = h+1; %count 1's 
                if (j<length(zeroTs) && zeroTs(j+1) == 0)
                    k = k+1;
                    zeroD(k) = h;
                end;
            else
                h = 0;
            end;
        end;
        if zeroD == 0
            zeroDall = [zeroDall 0];
        else
            zeroDall = [zeroDall zeroD(2:end)];
        end;
    end;
end;
ans20 = BootstrapMeanNoise(zeroDall,3000);
ansAll.durationzeroprotein = ans20(1:2); %mean and ste
