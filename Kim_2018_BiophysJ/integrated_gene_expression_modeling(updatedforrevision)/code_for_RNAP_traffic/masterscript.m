%{
-About-
This script contains command lines to run TASEPmodeling_(non)bursty and to analyze the run results

-Inputs-
see TASEPmodeling_(non)bursty

-varargin-

-Outputs-
see TASEPmodeling_(non)bursty
see TASEPmodeling_analysis

-Example-
masterscript
   
-Supplementary-

-Keywords-
TASEP modeling, gene expression

-Dependencies-

-References-

-Author-
Sangjin Kim, 2017 September 30
%}

% To calculate RNAP traffic from a bursty promoter with fOn = 0.25
lociList = TASEPmodeling_bursty('flat','P1',1,0.25,30,0,0,0);
headway = TASEPmodeling_analysis(lociList);

% To calculate RNAP traffic from a nonbursty promoter with kLoading =
% 1/15 sec-1;
lociList = TASEPmodeling_nonbursty('flat','P1nb',1,1/15,30,0,0,0);
headway = TASEPmodeling_analysis(lociList);

% To calculate RNAP traffic from a bursty promoter with fOn = 0.25
% with pause at 1500 nt, for 10 sec and probability = 80%
lociList = TASEPmodeling_bursty('OnepauseAbs','P1',1,0.25,30,1500,10,80);
headway = TASEPmodeling_analysis(lociList);

lociList = TASEPmodeling_nonbursty('OnepauseAbs','P1nb',1,1/15,30,1500,10,80);
headway = TASEPmodeling_analysis(lociList);
