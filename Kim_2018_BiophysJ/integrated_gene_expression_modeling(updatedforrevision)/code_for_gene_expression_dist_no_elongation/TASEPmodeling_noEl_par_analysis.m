function ansAll = TASEPmodeling_noEl_par_analysis(promoter,totalparRun)

%{
-About-
This function analyzes results of TASEPmodeling_noEl_(non)bursty_par in
case of the cells in the population contain only 1 gene copy per cell
This code combine results from separate parfor runs and analyze 
This generates data that is plotted in the paper

-Inputs-
promoter: name of the promoter (this will be included in the output file name)
totalparRun: total par file number

-varargin-

-Outputs-
ansAll: Transcription initiation rate (or intervals), headway
distribtuion, Distribution of mRNAs and
proteins per cell (mean, fano factor, CV2), 

-Example-
TASEPmodeling_noEl_par_analysis('P1',10)
   
-Supplementary-

-Keywords-
TASEP modeling, gene expression, no elongation, elongation free model

-Dependencies-
see masterscript

-References-

-Author-
Sangjin Kim, 2017 September 30
%}

%--------------------------------------------------------------------------

runCondition = 'NoEL-';

loadingSample = []; 
lifeTimeHist1 = []; lifeTimeAvg = []; 
fishSignal1 = []; 
tDiffStartHist = []; 
proteinLoci2 = [];
proteinSS = [];

% Load all par runs and combine together
for i = 1:totalparRun
    fileName = strcat(runCondition,promoter,'-',sprintf('%01.0f',i),'par.mat'); 
    tmp = load(fileName);
    loadingSample = [loadingSample, tmp.loadingSample];
    
    tDiffStartHist = [tDiffStartHist, tmp.tDiffStartHist];
    
    lifeTimeHist1 = [lifeTimeHist1 tmp.lifeTimeHist1]; 
    lifeTimeAvg = [lifeTimeAvg tmp.lifeTimeAvg];
    
    fishSignal1 = [fishSignal1 tmp.fishSignal1];
    
	proteinLoci2 = [proteinLoci2 tmp.proteinLoci2];
    proteinSS = [proteinSS tmp.proteinSS];

end;

sampleWindow = [15*60,30*60];

%% ans1 loadingSample = Number of successful initiation per DNA during sampleWindow
effectiveLoading = mean(loadingSample), std(loadingSample)
effectiveLoadingInterval=(sampleWindow(2)-sampleWindow(1))/effectiveLoading % in seconds
ans1 = BootstrapMeanNoise(loadingSample/(15*60),3000); %This is for effective loading rate (sec-1)
ans1(2,:) = BootstrapMeanNoise((15*60)./loadingSample,3000); %This is for effective loading interval (sec)
ansAll.tsxInitiationrate = ans1(:,1:2);

%% ans3 tDiffStartHist = temporal separation between RNAPs at the start 
binT = 0:5:200;
ans3(:,1) = binT';
ans3(:,2) = mean(tDiffStartHist,2);
ans3(end+1,:) = [205,0];
ansAll.RNAPheadway = ans3;

%% ans4 mRNA life-time
%% DNA with no RNAP initiation during sampleWindow will show mRNA life-time = 0
%% and hence should be eliminated from the analysis
tmp2 = find(loadingSample>0);
binLT = 0:5:1000;
ans4 = [];
ans4(:,1) = binLT';
ans4(:,2) = mean(lifeTimeHist1(:,tmp2),2); %life-time of 5'-end mRNA
ans4(:,3) = hist(lifeTimeAvg(1,tmp2),binLT)*100/length(tmp2); %life-time of 5'-end mRNA
ans4stat = BootstrapMeanNoise(lifeTimeAvg(1,tmp2),3000);
ansAll.mRNAlifetime = ans4;
ansAll.mRNAlifetimestat = ans4stat(:,1:2);

%% ans5/6 mRNA FISH signal
%% probe1 = 5'-end; probe2 = 3'-end; probe3 = 72 tiling probes
ssTime = size(fishSignal1,1); %steady-state time index of the fish time
ssfish1 = []; 
for i = ssTime
    ssfish1 = [ssfish1, fishSignal1(i,:)]; 
end;

% For mRNA distribution at steady-state
ans5 = [];
ans5(:,1) = (0:1:50)'; %binN';
ans5(:,2) = hist(ssfish1', ans5(:,1))*100/length(ssfish1);
ansAll.mRNAdist = ans5;

% For mean, fano, CV, CV^2
ans6(1,:) = BootstrapMeanNoise(ssfish1',3000);
ansAll.mRNAstat = ans6;

%% ans7/8 = Total protein number per DNA during 10 min
% For distribution
ans7(:,1) = (0:100:2000)'; %binP
ans7(:,2) = hist(proteinLoci2',ans7(:,1))*100/length(proteinLoci2);
ansAll.proteindist = ans7;

% For mean and fano (mean+/-ste; fano +/- ste; CV +/-ste, CV^2 +/- ste)
ans8 = BootstrapMeanNoise(proteinLoci2,3000);
ansAll.proteinstat = ans8;


%% examine steady-state protein levels
fishTime = 0:60:(size(fishSignal1,1)-1)*60;
figure, yyaxis left
xlabel('Time (sec)'); ylabel('mRNA per loci')
plot(fishTime, mean(fishSignal1,2), 'b-'); hold on; 

yyaxis right
ylabel('Protein per loci'); plot(fishTime, mean(proteinSS,2), 'r-'); hold off;
yyaxis left 
ylabel('mRNA per loci')

ssTime = size(fishSignal1,1); %steady-state time index of the fish time
ssProtein = [];
for i = ssTime
    ssProtein = [ssProtein, proteinSS(i,:)]; 
end;
ans9(1,:) = BootstrapMeanNoise(ssProtein',3000);
ansAll.steadystateprotein = ans9;




