function Xans = NETseqAnalysis(kind)
%{
-About-
This function analyze genome-wide gene expression properties from wig data.
 
-Inputs-
kind:  kind of genome-wide data. See examples
genome: reference genome. See examples

-varargin-

-Outputs-
Xans:  list of gene expression properties for each gene in the reference
gene.

-Example-
Xans = NETseqAnalysis('LarsonOld');
   
-Supplementary-
I supplied Larson's NETseq data files in the same folder (publication).

-Keywords-
NETseq, RNA seq, genome

-Dependencies-

-References-

-Author-
Sangjin Kim, 2016 October 13
%}

% Obtain gene location information from genebank
% Revisit MG if want to know some other information about genes
% getgenbank might not work on certain MATLAB versions.

% MG = getgenbank('U00096'); % newer
MG = getgenbank('U00096.2'); % used by Larson 
for i = 1:size(MG.CDS,2)
    if length(MG.CDS(i).indices) == 2
        geneLoc(i,1:2) = MG.CDS(i).indices;
    end;
end;
% save('MG1655geneLoc.mat','geneLoc'); 
% save and use later. skip the above for loop in that case.

geneLength = abs(geneLoc(:,1)-geneLoc(:,2))+1;
% general property test
for i = 1:size(geneLoc,1)
    x = NETseqreads(kind,geneLoc(i,:)); %get [coordinate, reads]
    readLength(i) = size(x,1);
    readLengthpercent(i) = size(x,1)/geneLength(i,1)*100;
    readSum(i) = sum(x(:,2));
    % this is mean read per base
    readSumnorm(i) = readSum(i)/geneLength(i); %assuming read = 0 at unread base
    readSumnorm2(i) = readSum(i)/size(x,1); %normalized to read length
end;
Xans = [readLength',readLengthpercent',readSum'];
end


function x = NETseqreads(dataType, geneLoc)
if geneLoc(1,1)>geneLoc(1,2)
    direction = '-';
else
    direction = '+';
end;
        
if strcmp(dataType, 'LarsonOLD')
    X = importdata(strcat('GSM1367304_Ecoli_WT.',direction,'.wig'));
elseif strcmp(dataType, 'LarsonBCMOLD')
    X = importdata(strcat('GSM1367308_Ecoli_WT_BCM.',direction,'.wig'));
end;
        
Y = X.data; clear X;
if geneLoc(1,1)>geneLoc(1,2)
    geneStart = find(Y(:,1) - geneLoc(1,2) >=0,1);
    geneEnd = find(Y(:,1) - geneLoc(1,1) <=0,1,'last');
else
    geneStart = find(Y(:,1) - geneLoc(1,1) >=0,1);
    geneEnd = find(Y(:,1) - geneLoc(1,2) <=0,1,'last');
end;
x = Y(geneStart:geneEnd,1:2);
end
