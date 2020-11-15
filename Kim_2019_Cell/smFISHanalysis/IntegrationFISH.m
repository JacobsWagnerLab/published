function intP = IntegrationFISH(rawFISHdata,mRNAlifetime)
%{
-About-
This function calculates the total mRNAs made from pulsed induction. The calculation follows eq 3 in SI of Kim, et al. LacZ paper.
It is to estimate premature termination between two probe regions.

-Inputs-
rawFISHdata:  three column data. First column = time points (minutes),
second column = Z5 signal, and third column = Z3 signal.
mRNAlifetime:  mean mRNA lifetimes (minutes), obtained by fitting an exponential decay function to the last 4-6 data points in the pulsed induction curve (i.e. t ~ 6-12 min)

-varargin-

-Outputs-
output:     total lacZ mRNAs (Z5 and Z3) made until each time point (first
column of rawFISHdata)

-Example-
rawFISHdata =
[0,0.0134487000000000,0.0181838000000000;1,1.79951000000000,0.0261977000000000;2,4.73399000000000,0.122877000000000;3,3.79987000000000,0.287993000000000;4,3.08398000000000,1.61866000000000;5,1.32696000000000,1.80267000000000;6,0.976509000000000,1.54928000000000;7,0.587405000000000,0.963588000000000;8,0.288642000000000,0.430795000000000;9,0.206146000000000,0.236055000000000;10,0.172012000000000,0.194890000000000;12,0.104088000000000,0.0946338000000000];
intP = IntegrationFISH(rawFISHdata, [1.52, 1.66]);

-Supplementary-

-Keywords-
lacZ FISH, signal integration, premature termination

-Dependencies-
smlacZFISHWorkflow

-References-

-Author-
Sangjin Kim, 2018 July 18

%}
t = rawFISHdata(:,1); S = rawFISHdata(:,2:3);
kc = 1./mRNAlifetime; %degradation rate (1/min)

%% get integration of P
for k = 1:size(S,2) %k is index for Z5 or Z3 or Y5
    intP(1,k) = 0;
    for i = 2:size(S,1)
        S_tmp = S(1:i,k); t_tmp = t(1:i,1);
        intP(i,k) = kc(k)*trapz(t_tmp,S_tmp)+S(i,k)-S(1,k); %equation (3) in SI
    end;
end;

figure,
plot(t,S(:,1),'ro-'); hold on;
plot(t,S(:,2),'bo-'); 
xlabel('Time (min)'); ylabel('mRNA numbers per cell');
plot(t,intP(:,1),'r-'); 
plot(t,intP(:,2),'b-'); 

