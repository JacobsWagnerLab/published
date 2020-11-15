%{
-About-
This script contains command lines to run
TASEPmodeling_noEl_bursty_par or TASEPmodeling_noEl_nonbursty_par and to analyze the
run results

-Inputs-
conditions for promoter activity (bursty, nonbursty, kinetic parameters), 
see TASEPmodeling_noEl_par

-varargin-

-Outputs-
output1: results of the TASEPmodeling_noEl_par runs will be saved as independent
matlab data files. See the last lines of TASEPmodeling_noEl_(non)bursty_par.m
for what is contained in those data files
output2: results of TASEPmodeling_noEl_par_analysis will be in the matlab workspace. See
TASEPmodeling_noEl_par_analysis for the details

-Example-
masterscript
   
-Supplementary-

-Keywords-
TASEP modeling, gene expression, no elongation, elongation free

-Dependencies-

-References-

-Author-
Sangjin Kim, 2017 September 30
%}

%--------------------------------------------------------------------------
% This "masterscript" is to run the TASEP model (NO ELONGATION) in a parallel fashion 
% This parfor-based code is advantageous for running 10000 simulations in a timely manner
% This parfor-based code only records statistics from each simulation and does not save raw trajectories
%% -----------------------------------------------------------------------------------------------------
%% Example scenario: 
% A bursty promoter, average RNAP loading interval ~15sec, no pause
close all; clear;
promoter='P1'; %arbitrary name, but will be included the output file name
fOn = 0.25; %fOn constant for promoter
for i = 1:10 %each run does 100 iterations
     TASEPmodeling_noEl_bursty_par(promoter,i,fOn);
end;
%% Result files are saved in the folder.

%% For analysis of the par results:
%% per DNA template or per simulations
TASEPmodeling_noEl_par_analysis('flat','P1',10,0, 0)


close all; clear;
promoter='P1nb'; %arbitrary name, but will be included the output file name
kLoading = 1/15;
for i = 1:10 %each run does 100 iterations
     TASEPmodeling_noEl_nonbursty_par(promoter,i,kLoading);
end;

TASEPmodeling_noEl_par_analysis('flat','P1nb',10,0, 0)
