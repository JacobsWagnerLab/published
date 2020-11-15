%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function cells = peakFinder_cellCheck(cellList,signal,fracOfMedian,distBtwnPks)
%@author: Molly Scott
%@date: March 4, 2016
%==========================================================================
%************************Output**********************:
%cells:                 a cell array where each cell contains the fields:
%       -pks:           signal values for peaks
%       -locs:          pixel location along cell for peaks
%       -length:        cell length (pixels)
%       -lengthvector:  pixel values for each step of the cell length
%       -normLength:    normalized lengthvector
%       -signal:        signal data for each step of the cell length
%************************Input**********************:
%cellList:              Output from Oufti
%signal:                'signal1' or 'signal2', etc.
%fracOfMedian:          how much above the median signal in the cell your peak must be to
%                       be counted (percent)
%distBtwnPks:           distance between peaks to be differentiated as peaks
%==========================================================================
%This function takes the input cellList (from Oufti) and determines, on a
%per cell basis, the peaks of signal fluorescence (i.e. zones) based on
%two parameters: (1) a minimum fluorescence relative to the median
%fluorescence of the cell and (2) a minimum distance between consecutive
%points of high fluorescence. The user scans through plots of the peaks of
%fluorescence as well as the peaks that have been called. The output is a
%cell array containing information about the location of peaks of
%fluorescence intensity in the cell (in pixels) as well as the fluorescence
%values of these peaks (in A.U.).
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 

function cells = peakFinder_cellCheck(cellList,signal,fracOfMedian,distBtwnPks)

%first, count how many cells we have with good data
cells = 0;
for ii = 1:length(cellList.meshData)
    for jj = 1:length(cellList.meshData{ii})
        if ~isempty(cellList.meshData{ii}{jj})
            if isfield(cellList.meshData{ii}{jj},signal)
                if ~isempty(cellList.meshData{ii}{jj}.(signal)) && length(cellList.meshData{ii}{jj}.(signal)) > 3
                    if max(cellList.meshData{ii}{jj}.(signal)) - min(cellList.meshData{ii}{jj}.(signal)) > 0.1
                        cells = cells + 1;
                    end
                end
            end
        end
    end
end

%create cell array for all cells
cells = cell(cells,1);
count = 0;

pause on
figure;

%fill cell array with data from each cell
for ii = 1:length(cellList.meshData)
    for jj = 1:length(cellList.meshData{ii})
        if ~isempty(cellList.meshData{ii}{jj})
            if isfield(cellList.meshData{ii}{jj},signal)
                if ~isempty(cellList.meshData{ii}{jj}.(signal)) && length(cellList.meshData{ii}{jj}.(signal)) > 3
                    if max(cellList.meshData{ii}{jj}.(signal)) - min(cellList.meshData{ii}{jj}.(signal)) > 0.001
                        if cellList.meshData{ii}{jj}.length > distBtwnPks
                            count = count + 1;
                            signalValue = smooth(double(cellList.meshData{ii}{jj}.(signal)),12); 
                            [pks,locs] = findpeaks(signalValue,'MINPEAKHEIGHT',median(cellList.meshData{ii}{jj}.(signal))+fracOfMedian*median(cellList.meshData{ii}{jj}.(signal)),'minpeakdistance',distBtwnPks);
                            for ps = 1:length(pks) %let's double check the peaks and make sure they're real local maxima
                                k = find(signalValue == pks(ps),1);
                                if (k > 20) && (length(signalValue) > k+20);
                                    sgn_range = signalValue(k-20:k+20);
                                    if pks(ps) < median(sgn_range) + std(signalValue)
                                        pks(ps) = NaN;
                                        locs(ps) =NaN;
                                    end
                                end
                            end
                            locs = locs(~isnan(locs)); %remove false positive peaks
                            pks = pks(~isnan(pks));
                            %fill in cell array
                            cellLength = cellList.meshData{ii}{jj}.lengthvector;
                            normLength = (cellLength - min(cellLength))/(max(cellLength) - min(cellLength));
                            cells{count}.pks = pks;
                            cells{count}.locs = locs;
                            cells{count}.length = cellList.meshData{ii}{jj}.length;
                            cells{count}.lengthvector = cellLength;
                            cells{count}.normLength = normLength;
                            cells{count}.signal = cellList.meshData{ii}{jj}.(signal);

                            %plot the result: peaks and signal values for each
                            %cell. Scan thru cells by clicking any key
                            if length(cellList.meshData{ii}{jj}.(signal)) == length(cellLength)
                                plot(cellLength,cellList.meshData{ii}{jj}.(signal),'b'); hold on
                                plot(cells{count}.lengthvector(locs),cells{count}.pks,'k^','markerfacecolor','r')
                                xlabel('Cell length (pixels)','FontSize',18)
                                ylabel('Fluorescence (A.U.)','FontSize',18)
                                legend('Signal','Peaks')
                            end
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