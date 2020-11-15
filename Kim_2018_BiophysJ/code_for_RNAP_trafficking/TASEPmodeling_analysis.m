function headwayList = TASEPmodeling_analysis(lociList)

%{
-About-
This function uses results of TASEPmodeling_(non)bursty for plotting
RNAP trajectories and analyzing RNAP headways

-Inputs-
lociList: results from TASEPmodeling_(non)bursty, contatining RNAP
trajectories from each DNA template

-varargin-

-Outputs-
output1: Figure, plot of RNAP trajectories
headwayList: headway and d_headway (delta headway, headway change) 

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

%% To plot, RNAP trajectories from example trajs
geneLength = 3075;
xLattice = 1:geneLength; xLattice = xLattice';
for lociID = 1:1 %<- change this as you want
    exitTime = lociList(lociID).polexittime;    
    figure,
    for polID = 1:size(exitTime,2)
        if exitTime(1,polID)>0
            plot(exitTime(1:geneLength,polID),xLattice); hold on;
            axis([1200,1800,0,3075]);
            xticks([1200,1500,1800]);
            yticks([0 1000 2000 3000]);
        end;
    end;
    xlabel('Time (sec)'); ylabel('Position (nt)');
end;


%% RNAP headway
sampleWindow = [15*60,30*60]; 
N_totalloci = size(lociList,2); 
geneLength = 3075;
dtStart = []; dtEnd = [];
for lociID = 1:N_totalloci
    exitTime = lociList(lociID).polexittime ;
    polStatus = lociList(lociID).polstatus;
    loadTime = lociList(lociID).loadtime;
    lastPol = max(find(exitTime(1,:)>0));
    tStart = []; tEnd= [];
    for polID = 1: lastPol
        if exitTime(1,polID) >=0 %realTraj
            if exitTime(1,polID)>=sampleWindow(1) && exitTime(1,polID)<sampleWindow(2)
                tStart = [tStart loadTime(1,polID)];
                tEnd = [tEnd exitTime(geneLength,polID)];
            end;
        end;
    end;
    dtStart = [dtStart diff(tStart)]; %headway between RNAPs at initiation (at the promoter)
    dtEnd = [dtEnd diff(tEnd)]; %headway between RNAPs at termination (at the end of gene)
end;
dHeadway = dtEnd - dtStart; %headway change
headwayList = [dtStart;dtEnd;dHeadway]';