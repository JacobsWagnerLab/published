%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function sigmas = calculate_width(cellList,d_indx,scope_sigma,signal,pkThresh,lengthThresh,pix2mu)
%@author: Molly Scott
%@date: March 4, 2016
%==========================================================================
%************************Output**********************:
%sigmas:                matrix containing the values for the sigma of the fit to the
%                       peaks of fluorescence intensity (in microns)
%************************Input**********************:
%cellList:              output from Oufti
%d_indx:                the range of pixel values around the peak of fluorescence that
%                       will be considered for fitting
%scope_sigma:           experimentally determined PSF for your scope
%signal:                'signal1' or 'signal2'
%pkThresh:              minimum fluorescence value above the median to be considered
%                       a 'peak' (a way to verify that you are fitting real peaks)
%lengthThresh:          minimum length of cells to be evaluated. Here, this is
%                       used in the E. coli dataset to get a first pass thresholding for cells
%                       nearing the end of their cell cycle.
%pix2mu:                pixel to micron conversion
%==========================================================================
%This function allows you to calculate the sigma for the width of a fit of
%a convolved rectangular function to peaks of fluorescence intensity. It
%relies on the function fitErf.m to call the fitting function.
%-------------------------------------------------------------------------- 
%--------------------------------------------------------------------------
function sigmas = calculate_width(cellList,d_indx,scope_sigma,signal,pkThresh,lengthThresh,pix2mu)

pause on; % this lets us scan through each figure output by pressing a key
 figure;

 count = 0;
%  d_indx = 15;
%  sigma = 1.593;
 kk = 0;
 
for frame = 1:length(cellList.meshData)
    for cell = 1:length(cellList.meshData{frame})
        if ~isempty(cellList.meshData{frame}{cell}.length)
            if ~isempty(cellList.meshData{frame}{cell}.(signal))
                if cellList.meshData{frame}{cell}.length > lengthThresh %first threshold by min cell length
                 cellLength = cellList.meshData{frame}{cell}.lengthvector;
                 sig = cellList.meshData{frame}{cell}.(signal); 
                 %now, verify that there's some type of peak in
                 %fluorescence
                 if max(sig) > median(sig) + pkThresh
                    [vm,im] = max(sig)
                    if im > d_indx
                        indx1=im-d_indx; indx2=im+d_indx; % find the range of pixels around the peak
                        if indx2 < length(sig)
                            signal_cut=sig(indx1:indx2); % let's cut the signal to our region of interest +/- 50 pixels
                            XX=(1:length(signal_cut))'/scope_sigma/sqrt(2); % the x axis for the fit (the cell length axis) automatically multiplies our input by the conversion factor necessary to account for the "matlab erf" function (basically, our derived function was close, but had different limits of integration-- this factor accounts for that change)


                            [coeff,fit1,YYfit]=fitErf(XX,signal_cut,'2_erf'); % apply the fit for XX = the cell length for our shortened region, signal_cut = signal intensity values contained within that region of the cell length, 2_erf = our fitting error function
                            disp(fit1) % show me the coefficients for each fit (a,w,b,x0) as I scroll through each output figure (will be updated each scroll)
                            ci = confint(fit1); % determine the confidence interval for each coefficient
                            kk=kk+1; % start a count so that each loop I start on a different row, in this case
                            aa(kk,:)=[coeff(1),ci(1,1),ci(2,1)]; % store an array with a value plus its 95% confidence interval
                            ww(kk,:)=[coeff(2),ci(1,2),ci(2,2)]; % store an array with w value plus its 95% confidence interval
                            bb(kk,:)=[coeff(3),ci(1,3),ci(2,3)]; % store an array with b value plus its 95% confidence interval
                            xx(kk,:)=[coeff(4),ci(1,4),ci(2,4)]; % store an array with x0 value plus its 95% confidence interval
                            cl(kk)=cellList.meshData{frame}{cell}.length; % store an array with the cell's length
                            fluo1(kk)=sum(cellList.meshData{frame}{cell}.signal1); % store an array with the cell's total fluorescence
                            fframe(kk)=frame; % tell me the frame for the signal
                            ccell(kk)=cell; % tell me the cell signal
                            
                            plot(XX,signal_cut,'ob'); hold on
                            plot(XX,YYfit,'-r'); % plot the real values of signal int with blue circles, the fit with red lines
                            title(['frame=',num2str(frame),' ','cell=',num2str(cell),', ','peak number=',num2str(kk)])  % title by converting the number to a string 
                            xlabel('Cell length (pixels)','FontSize',18)
                            ylabel('Fluorescence (A.U.)','FontSize',18)
                            legend('Signal','Fit')
                            pause
                            hold off
                        end
                    end
                 end
                end
            end
        end
    end
end
sigmas = ww(:,1) * pix2mu;
end