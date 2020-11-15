function Xans = RNApercell_by_totalfluorescence_bootstrap(dayList,daily_scale,daily_avgbasePix,timecourseInfo)
%{
-About-
This code calculates mRNA numbers per cell based on total fluorescence signal within cell mesh.
This code do bootstrap by sampling an equal number of cells from different
day's experiments.

-Inputs-
dayList: contains file path of each day
daily_scale: single mRNA normalization factor, different for each day.
daily_avgbasePix: mean fluorescence intensity within a pixel at time zero.
this is used for background subtraction
timecouseInfo: actual time (minutes) of FISH sample; day and the file name
for the images from that time.

-varargin-

-Outputs-
Xans: Bootstrap result of mean mRNA number per cell in each time point.
standard error of mean, calculated from bootstrapping, and the number of
cells used for calculation

-Example-
Xans = RNApercell_by_totalfluorescence_bootstrap(dayList,daily_scale,daily_avgbasePix,timecourseInfo)
see smlacZFISHworkflow for the input
  
-Supplementary-

-Keywords-

-Dependencies-
smlacZFISHworkflow
MicrobeTracker

-References-

-Author-
Sangjin Kim, 2016, March 20
%}
timePoints = timecourseInfo.time;
days = timecourseInfo.day;
files = timecourseInfo.file;

bsSteps = 3000; %1 for no iteration, 3000 for real bootstrap
cellNumperday = zeros(length(timePoints),1);

for ti = 1:length(timePoints)
    int1 = []; 
    int2 = []; 
    
    for di = 1: length(days{ti})
        clear dayfolder;
        dd = days{ti}(di);
        dayfolder = dayList(dd).folderpath;
        scaleZ5 = daily_scale(dd,1); scaleZ3 = daily_scale(dd,2);
        avgbasePix1 = daily_avgbasePix(dd,1); avgbasePix2 = daily_avgbasePix(dd,2);

        tifile = strcat('t',sprintf('%01.0f',files{ti}(di)),'-meshaCy5.mat');     
                
        cellListfile = strcat(dayfolder,tifile);
        load(cellListfile);
        cellarea = AreaHist(cellList); %% "areahist" made based on lengthhist
        tempint = (inthist(cellList,'signal1')-avgbasePix1*cellarea)/scaleZ5; 
        int1{di} = tempint; 
        tempint = (inthist(cellList,'signal2')-avgbasePix2*cellarea)/scaleZ3; 
        int2{di} = tempint; 
        cellNumperday(ti,dd) = length(tempint);
    end;
    
    if sum(cellNumperday) == 0, 
        continue; 
    else
        nboot = min(cellNumperday(ti,(find(cellNumperday(ti,:)>0)))); %choose min cellNum per day.
        %or
        %nboot = max(cellNumperday(ti,(find(cellNumperday(ti,:)>0)))); %choose max cellNum per day.
    end;
            
    bootint = 0;
    for i = 1:bsSteps
        clear samplecellindex;
        tempint1 = []; tempint2 = [];
        for di = 1:length(days{ti})
            dd = days{ti}(di);
            if cellNumperday(ti,dd)>0
                samplecellindex = unidrnd(cellNumperday(ti,dd),nboot,1);
                tempint1 = [tempint1, int1{di}(samplecellindex)];
                tempint2 = [tempint2, int2{di}(samplecellindex)];
            end
        end
        % bootint is mean of int in each bs iteration
        bootint(i,1) = mean(tempint1); %mean of int1 of selected cells
        bootint(i,2) = mean(tempint2); %mean of int2 of selected cells
    end  
    Xans(ti,1) = mean(bootint(:,1)); 
    Xans(ti,2) = std(bootint(:,1));
    Xans(ti,3) = mean(bootint(:,2)); 
    Xans(ti,4) = std(bootint(:,2));
    Xans(ti,5) = nboot*length(find(cellNumperday(ti,:)>0)); 
end
end %end of the main function

% Function to calculate cell area, based on lengthhist function in
% MicrobeTracker
function varargout = AreaHist(varargin)
% lengthhist(cellList)
% lengthhist(cellList1,cellList2,...)
% lengthhist(cellList,xarray)
% lengthhist(cellList1,cellList2,xarray)
% lengthhist(cellList1,pix2mu)
% lengthhist(cellList1,pix2mu1,cellList2,pix2mu2)
% lengthhist(cellList1,cellList2,'overlap') 
% lengthlist = lengthhist(cellList)
% lengthlist(...'nooutput')
% lengthlist(...'nodisp')
% [lengthlist1,lengthlist2] = lengthhist(cellList1,cellList2)
% 
% This function plots a histogram of the length of every cell in a population
% 
% <cellList> is an array that contains the meshes. You can drag and drop 
%     the file with the data into MATLAB workspace or open it using MATLAB
%     Import Tool. The default name of the variable is cellList, but it can
%     be renamed.
% <cellList1>, <cellList2> - you can load two arrays, they will be plotted
%     together for comparison.
% <xarray> - array of x values for the histogram, which serve the the
%     centers of bins of the histogram (the boundaries will be in between,
%     for example [1 2 3 4 5] to display all the cells shorter than 1.5 
%     micron in the first bin, between 1.5 and 2.5 microns in the second, etc.).
% <pix2mu>, <pix2mu1>, <pix2mu2> - conversion factors from pixels to 
%     microns, the size of a pixel in microns. Typical value 0.064.
% 'overlap' - indicate this if you wish the histograms to overlap,
%     otherwise they will be displayed separately.
% <lengthlist>, <lengthlist1>, <lengthlist2> - arrays containing the length
%     of every cell to save and plot separately.
% 'nooutput' - blocks standard output (type of data processed, mean, and
%     standard deviation).
% 'nodisp' - suppresses displaying the results as a figure.

expr = 'value = cellList{frame}{cell}.area;';
xlabel1 = 'Cell area, pixels';
xlabel2 = 'Cell area, \mum';
areaarray = plothist(expr,xlabel1,xlabel2,varargin);
for i=1:nargout, varargout{i} = areaarray{i}; end

end