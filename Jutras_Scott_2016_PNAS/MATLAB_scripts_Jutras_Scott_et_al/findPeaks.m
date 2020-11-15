function peaks = findPeaks(signals, fieldname, doFit)
%This script is dependent on on output from read_cellLists.m to generate
%the signals cell array first.
%finds peak locations in the signals array for the signal specified by
%fieldname. if dFit = true, the peak will be fit with a biased gaussian. if
%doFit = false, the peak location will be returned as its center of mass.
%
%findPeaks begins by looking for a peak within searchR of the signal
%center. Any peaks found must be the highest amplitude point within a
%radius of searchRlocal. Initial smoothing is performed to reduce the
%effect of noise on the initial guess for a peak. If no peak is found, non
%is inserted into the peaks output array
%
%Brad Parry, Christine Jacobs-Wagner lab; April 2016

% fieldname = 'hada';
% doFit = true

%these parameter values are dependent on the optical set up, and in
%particular, the pixel size of your camera:
searchR = 50; %initial guess radius from cell center
searchRlocal = 20; %peak must be brightest point in this radius
sWin = 19; %smoooooth the data first with a conv of this size.

%Biased Gaussian that will be used for fitting. The reason for bias is that
%signals are typically present on an uneven background. If background is
%even, m goes to 0 and the first term drops out.
gaussFn = fittype('m*x+A+K/(sigma*sqrt(2*pi))*exp(-(1/2)*((x-u)/sigma).^2)',...
    'independent',{'x'},'dependent','Y');
peaks = [];
for C = 1:length(signals)
    
    %pull the signal at fieldname and smooth it
    S = signals{C}.(fieldname);
    %will focus on the center, so keep it the same size to get around
    %indexing challenges
    Sconv = conv(S, ones(1,sWin),'same');

    %grab the central region to search for a peak
    inds = round(length(S)/2) + (-searchR:searchR);
    inds(inds<1) = [];
    inds(inds>length(Sconv)) = [];
    
    %lets make sure that this is indeed a peak
    [~, ix] = max(Sconv(inds));
    ix = ix + inds(1) - 1;
    inds = ix + (-searchRlocal:searchRlocal);
    ind_center = inds(searchRlocal+1);
    inds(inds<1) = [];
    inds(inds>length(Sconv)) = [];
    
    if sum(Sconv(ind_center) > Sconv(inds)) ~= (length(inds) - 1)
        %then we don't have a peak!!!!
        peaks(end+1) = nan;
        tx = ['no peak found ',num2str(C)];
        disp(tx)
        continue
    end
    
    %if we are here, it is a peak by our criteria.
    peakWindow = Sconv(inds);
    %extract the "peak-iest" part of the peak. we'll just fit that
    [ix0, ix1] = trimPeak(peakWindow);
    
    Y = peakWindow(ix0:ix1);
    xFit = 1:length(Y);
    xFit = xFit(:);
    Y = Y(:);

    if doFit
        [~,mx] = max(Y);
        sigma = abs( (xFit(mx) - max( xFit([1,length(xFit)]) )) / 2 );
        [sfit, ~] = fit(xFit,Y,gaussFn, 'startpoint', [min(xFit),1,0,sigma, xFit(mx)]);

        %the fit can be evaluated:
        %yhat = gaussFn(sfit.A,sfit.K,sfit.m,sfit.sigma, sfit.u, xFit);

        peaks(end+1) =  sfit.u  + inds(1) + ix0 - 2;
    else
        %compute centroid
        peaks(end+1) = sum(xFit.*Y) / sum(Y) + inds(1) + ix0 - 2;
    end

end
end

function [ix0, ix1] = trimPeak(peakWindow)

    %trim the peak to the narrowest region with derivatives of the greatest
    %magnitude. 
    [ix0] = trimPeak_(peakWindow);
    [ix1] = trimPeak_(peakWindow(length(peakWindow):-1:1));
    ix1 = length(peakWindow) - ix1 + 1;
    
end

function [ix0] = trimPeak_(peakWindow)
%This sub-function discretely calculates the first derivative of
%peakWindow and then searches for the final negative derivative -- that is,
%the location corresponding to the base of the peak. If none is found, the
%smallest derivative is used instead as the location.

    d1 = diff(peakWindow);
    [~, ix] = max(peakWindow);

    %this is the fragment of derivatives up to the peak
    [~,mx] = max(d1(1:(ix-1)));

    ix0 = find(d1(1:mx) < 0);
    if isempty(ix0)
        [~,ix0] = min(d1(1:mx));
    end
    %add back one because we were indexing from the derivative
    ix0 = ix0(end) + 1;
    
end
