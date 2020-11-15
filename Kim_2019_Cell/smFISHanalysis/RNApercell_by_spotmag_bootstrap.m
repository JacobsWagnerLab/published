function Xans = RNApercell_by_spotmag_bootstrap(dayList,daily_scale,daily_thresh,timecourseInfo)
%{
-About-
This code calculates mRNA numbers per cell based on total spot intensities. Spot magnitude is normalized
by a spot intensity of a single mRNA and spots within a cell is summed.
This code do bootstrap by sampling an equal number of cells from different
day's experiments.

-Inputs-
dayList: contains file path of each day
daily_scale: single mRNA normalization factor, different for each day.
daily_thresh: threshold for spot curation, different for each day.
timecouseInfo: actual time (minutes) of FISH sample; day and the file name
for the images from that time.

-varargin-


-Outputs-
Xans: Bootstrap result of mean mRNA number per cell in each time point.
standard error of mean, calculated from bootstrapping, and the number of
cells used for calculation

-Example-
Xans = RNApercell_by_spotmag_bootstrap(dayList,daily_scale,daily_thresh,timecourseInfo)
see smlacZFISHworkflow for the input
  
-Supplementary-

-Keywords-

-Dependencies-
smlacZFISHworkflow

-References-

-Author-
Sangjin Kim, 2013, January 10
%}

timePoints = timecourseInfo.time;
days = timecourseInfo.day;
files = timecourseInfo.file;
SMnormalization = 1; % 1= do normalization in spotThresholding function. Choose 0 to avoid normalization

% Go over each time point
for ti = 1:length(timePoints)
    cellindex = []; 
    cellNum = [];
    magnitudeArray1a = []; 
    magnitudeArray1b = []; 
    
    for di = 1:length(days{ti})
        dd = days{ti}(di);
        dayfolder = dayList(dd).folderpath;
        scaleZ5 = daily_scale(dd,1); scaleZ3 = daily_scale(dd,2);
        thresh1 = daily_thresh(dd,1); thresh2 = daily_thresh(dd,2);
        
        tifile = strcat('t',sprintf('%01.0f',files{ti}(di)),'-mesha');
        cellListfile = strcat(dayfolder,tifile,'Cy5.mat'); load(cellListfile); TEMPcellLista = cellList;
        cellListfile = strcat(dayfolder,tifile,'Cy3.mat'); load(cellListfile); TEMPcellListb = cellList;
        TEMPcellindex = CellCounter(TEMPcellLista,dd); %take [day, frame, cellNum] of non-empty cell in the cellList
        cellindex{dd} = TEMPcellindex(:,2:3);
        cellNum(dd) = length(TEMPcellindex); % total 'real' cell # in the cellList

        if SMnormalization == 0
            scaleZ5 = 1; scaleZ3 = 1;
        end;

        magnitudeArray1a{dd} = SpotThresholding(TEMPcellLista,1,thresh1,scaleZ5);
        magnitudeArray1b{dd} = SpotThresholding(TEMPcellListb,2,thresh2,scaleZ3);
        clear TEMP*;
    end;
   
    close all;
    nboot = min(cellNum(find(cellNum>0))); %equal number of cells withdrawn each day
    % can change to..nboot = max(cellNum(find(cellNum>0)));
    % the result does not change much
    
    numberOfIterations = 3000; %for bootstrapping
    bootspot1CHa = zeros(numberOfIterations,6);
    bootspot1CHb = zeros(numberOfIterations,6);
    parallelobject = parpool(12); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parfor i = 1:numberOfIterations   %change parfor to for if debugging
        % sample and generate test manitudeArray to calculate statistics
        % (bootspot1CH).
        [magnitudeArraya,magnitudeArrayb] = MakeTestArrays(cellindex,magnitudeArray1a,magnitudeArray1b);
        bootspot1CHa(i,:) = SpotStat1CH(magnitudeArraya{1});
        bootspot1CHb(i,:) = SpotStat1CH(magnitudeArrayb{1});
    end
    delete(parallelobject); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %bootspot1CHa(:,6) = mean total spotmag per cell (Z5), calculated 3000 times
    Xans(ti,1) = mean(bootspot1CHa(:,6)); Xans(ti,2) = std(bootspot1CHa(:,6));
    
    %bootspot1CHb(:,6) = mean total spotmag per cell (Z3), calculated 3000 times
    Xans(ti,3) = mean(bootspot1CHb(:,6)); Xans(ti,4) = std(bootspot1CHb(:,6));
    
    Xans(ti,5) = nboot*length(find(cellNum>0)); 
    %total number of cells analyzed from the equal number of cells
    %sampled each day (nboot)
end

end % for this function


function cellindex = CellCounter(cellList,n)
i = 0;
for frame=1:length(cellList)
    for cell=1:length(cellList{frame})
        if ~isempty(cellList{frame}{cell}) && length(cellList{frame}{cell}.mesh)>4
            i = i+1;
            cellindex(i,1) = n;
            cellindex(i,2) = frame;
            cellindex(i,3) = cell;
        end
    end
end
end %for cellcounted

function magnitudeArray = SpotThresholding(cellList,ch,thresh,scale)
sizeOfCellList  = cellfun(@length,cellList);
magnitudeArray  = cell(1,length(cellList)); 

for frame=1:length(cellList)
    magnitudeArray{frame}  = cell(1,sizeOfCellList(frame));
    for cellNum=1:length(cellList{frame})
        if ~isempty(cellList{frame}{cellNum}) && length(cellList{frame}{cellNum}.mesh)>4 && isfield(cellList{frame}{cellNum},'spots') 
            %---for given cell, spot threshold
            if ch == 1
                prf = cellList{frame}{cellNum}.signal1./cellList{frame}{cellNum}.steparea;
            elseif ch == 2
                prf = cellList{frame}{cellNum}.signal2./cellList{frame}{cellNum}.steparea;
            end;
            
            for i=1:3
                prf = prf*0.5+prf([1 1:end-1])*0.25+prf([2:end end])*0.25; %smoothing
            end
            
            pos = cellList{frame}{cellNum}.spots.positions;
            prf = [0;prf];
            ind = prf(pos+1)>thresh; %1 or 0
            ind2 = cellList{frame}{cellNum}.spots.magnitude>0 & cellList{frame}{cellNum}.spots.magnitude<20;
            ind2 = ind2';
            try
                ind = ind & ind2; %logic numbers
            catch
            end
            magnitudeArray{frame}{cellNum} = cellList{frame}{cellNum}.spots.magnitude(ind)/scale;
        end
    end
end
end %for spotThresholding

function [magnitudeArraya,magnitudeArrayb] = MakeTestArrays(cellindex,magnitudeArray1a,magnitudeArray1b)
for di = 1:length(cellindex)
    try
        cellNumperday(di) = length(cellindex{di});
    catch
    end
end
nboot = min(cellNumperday(find(cellNumperday>0)));
%or
%nboot = max(cellNumperday(find(cellNumperday>0)));

testcellindex = []; %day, frame, cell numbers
for di = 1:length(cellNumperday) %go over each day
    if cellNumperday(di)>0
        %randomly choose nboot cells from a given day
        sampleindx = unidrnd(cellNumperday(di),nboot,1);
        sampleindx2 = zeros(nboot,3);
        sampleindx2(:,1) = di; % for day
        sampleindx2(:,2:3) = cellindex{di}(sampleindx,:); %for frame and cell
        testcellindex = [testcellindex; sampleindx2];
    end
end

for j = 1:size(testcellindex,1) %generate test magnitude list
    magnitudeArraya{1}{j}  = magnitudeArray1a{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
    magnitudeArrayb{1}{j}  = magnitudeArray1b{testcellindex(j,1)}{testcellindex(j,2)}{testcellindex(j,3)};
end

end %for makeTestArrays

function Xans = SpotStat1CH(magnitudeArray)
% magnitudeArray is of the length of cells
% find cells with spot magnitude.
spotArrayNonEmpty = find(~cellfun('isempty',magnitudeArray)>0); 
tempMagnitudeArray = magnitudeArray(spotArrayNonEmpty);
tempMagnitudeArray = cell2mat(tempMagnitudeArray);

if length(tempMagnitudeArray)>0
    spotList = zeros(length(tempMagnitudeArray),3); %length of spots, or the number of spots
    spotList(:,1) = tempMagnitudeArray;

    %sum spot magnitudes in each cell
    totalmagpercell = cell2mat(cellfun(@sum,magnitudeArray,'uniformoutput',0));
    %count number of spots per cell
    spotspercell = cell2mat(cellfun(@length,magnitudeArray,'uniformoutput',0));
    %divide the number of cells with spot(s) by the total number of cells
    cellWithSpot = length(spotArrayNonEmpty)*100/length(magnitudeArray); %%%%%%%%

    Xans(1) = mean(spotList(:,1)); % mean spot intensity
    Xans(2) = std(spotList(:,1)); % std of spot intensity
    Xans(3) = length(spotList); % total spot number
    Xans(4) = cellWithSpot;  %cells with spot %
    Xans(5) = mean(spotspercell); % spot # per cell including zeros
    Xans(6) = mean(totalmagpercell);
else
    Xans(1) = NaN; Xans(2) = NaN; Xans(3) = 0; Xans(4) = 0; Xans(5) = 0; Xans(6) = 0;
end

end %for spotStat1CH

