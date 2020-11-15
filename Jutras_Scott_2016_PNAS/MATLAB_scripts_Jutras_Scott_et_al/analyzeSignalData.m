%{
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Many data transformations performed in the script or its dependencies make
particular assumptions regarding the optical setup that images were
acquired on (in particular, the camera bit depth, camera pixel size and
diffraction limit; this list may not be exhaustive). If you are using this
script for your own analysis, there is no guarantee that it will perform as
expected!

REQUIREMENTS
This script requires that the cell array "signals" has been assembled from 
read_cellList.m
Additionally, there many code dependencies which have been packaged with
this script. Keep them in the same folder or accessible in your path.

FUNCTION
This script groups the data found in signals based on similarity and weight the 
contribution of each signal to its group according to its similarity to the bin 
center. 

Each cell in the signals cell array contains line scans measured from the cell 
centerline of cytoplasmic gfp, dna, hada and the inverted phase contrast signal 
from single biological cells. Features are extracted from the gfp and dna 
signals and assembled based on similarity. Prior to similarity assembly, all 
signals in each biological cell are shifted so that the central peak in the hada
signal occurs exactly in the middle. If there is no central hada peak (n = 22 
cells in the current dataset), no alignment based on hada is performed.

Brad Parry, Christine Jacobs-Wagner lab; 2016 April
%}

%load signals if stored elsewhere
% load('cellList_signals.mat')

%all output will be printed if set to true, if not desired, set to 'false'
printOutput = true;
%the number of bins the data should be segregated into
Nbins = 10;
%anonymous function for normalization. In the context of this script, this
%function is used to normalize signals from the cellList. Prior to doing
%so, you must make sure that the spacing between all signal indices is
%identical. Otherwise, this summation normalization will give poor results.
signalNorm = @(x) x / sum(x(:));

%
%   Identify signal features
%
%determine how constricted the dna and gfp signals are using a window of 9
%pixels
sWin = 9;
[gfpConstriction] = signalConstriction(signals,'gfp',sWin);
[dnaConstriction] = signalConstriction(signals,'dna',sWin);
%identify hada peaks near the cell center
peaks = findPeaks(signals, 'hada', false);
%get an array of cell lengths
cellLengths = cellfun(@(x) x.cellLength, signals);

%Construct a smoothed representation of the data which will be used to
%construct bins
[x, y] = gaussianWeightX(gfpConstriction, dnaConstriction, .3);
%measure the cumulative path length of x,y
z = [0,cumsum((diff(x).^2 + diff(y).^2).^(1/2))];
%construct Nbins on line x,y that are eqidistant in that space
bins = 0:(max(z)/(Nbins-1)):max(z);
binInds = [];
for k = bins
    [~,ix] = min( abs(z-k) );
    binInds(end+1) = ix(1);
end
Xbin = x(binInds);
Ybin = y(binInds);

%Display the bins and color point weights
cm = parula(301);
sigma = .15;
normalizationFn = 1/(2*pi*sigma^2);
for k = 1:length(Xbin)
    dx = Xbin(k) - gfpConstriction;
    dy = Ybin(k) - dnaConstriction;
    weights = normalizationFn * exp(-(dx.^2 + dy.^2)/(2*sigma^2));
    
    %impose a threshold by cutting off weights under half max of the
    %probability feature
    ix1 = find(weights > (1/(2*pi*sigma^2))/2);
    ix0 = setdiff(1:length(weights), ix1);
    
    figure
    hold on
    plot(x,y,'r')
    plot(gfpConstriction(ix0), dnaConstriction(ix0),'o','color',[.7,.7,.7])
    for w = ix1
        q = (weights(w) - (1/(2*pi*sigma^2))/2) / ((1/(2*pi*sigma^2))/2);
        q = floor(q*300) + 1;
        plot(gfpConstriction(w),dnaConstriction(w),'o','color',cm(q,:))
    end

     set(gca,'fontsize',16)
     
    if printOutput
        tx = ['bin',num2str(length(Xbin)-k+1)];
        print(tx,'-dpdf')
    end
end

gfp = cell(1,length(Xbin));
phs = cell(1,length(Xbin));
dna = cell(1,length(Xbin));
hada = cell(1,length(Xbin));
muL = []; %average cell length
W = cell(1,length(Xbin)); %signal weights

for k = 1:length(Xbin)
    
    dx = Xbin(k) - gfpConstriction;
    dy = Ybin(k) - dnaConstriction;
    weights = 1/(2*pi*sigma^2) * exp(-(dx.^2 + dy.^2)/(2*sigma^2));
    ix = find(weights > (1/(2*pi*sigma^2))/2);
    W{k} = weights(ix)';

    for C = ix
        if isnan(peaks(C))
            %if a peak was not found for the HADA signal, append the result
            phs{k}(end+1,:) = pchip_interp(signalNorm( signals{C}.phs )');
            gfp{k}(end+1,:) = pchip_interp(signalNorm( signals{C}.gfp )');
            dna{k}(end+1,:) = pchip_interp(signalNorm( signals{C}.dna )');
            hada{k}(end+1,:) = pchip_interp(signalNorm( signals{C}.hada )');
            continue
        end
        %if a peak was found for the HADA signal, shift all signals and
        %append
        phs{k}(end+1,:) = center_interp_signal(signalNorm( signals{C}.phs )', peaks(C));
        gfp{k}(end+1,:) = center_interp_signal(signalNorm( signals{C}.gfp )', peaks(C));
        dna{k}(end+1,:) = center_interp_signal(signalNorm( signals{C}.dna )', peaks(C));
        hada{k}(end+1,:) = center_interp_signal(signalNorm( signals{C}.hada )', peaks(C));
    end
    
    phs{k}(isnan(phs{k})) = 0;
    gfp{k}(isnan(gfp{k})) = 0;
    dna{k}(isnan(dna{k})) = 0;
    hada{k}(isnan(hada{k})) = 0;
   
    muL(end+1) = mean( cellLengths(ix) ) * 0.0642;
end

%anonymous functions to apply weights to signals and calculate weighted
%averages
Wmean = @(x,w) sum(x.*repmat(w,[1,size(x,2)])) ./ sum( (x~=0).*repmat(w,[1,size(x,2)]) );
mn = @(x,w) Wmean(x,w) / sum(Wmean(x,w));

%plot weighted averages for each bin
for k = 1:length(gfp)
    figure
    hold on
    plot(mn(phs{k},W{k}),'color',[.7,.7,.7])
    plot(mn(gfp{k},W{k}),'g')
    plot(mn(dna{k},W{k}),'b')
    plot(mn(hada{k},W{k}),'k')
    set(gca,'fontsize',16)
    tx = ['\langleCell length\rangle = ',num2str(muL(k)),' \mum'];
    xlabel(tx)
    ylabel('Signal, AU')
    
    if printOutput
        tx = ['binnedData',num2str(length(gfp)-k+1)];
        print(tx,'-dpdf')
    end
end