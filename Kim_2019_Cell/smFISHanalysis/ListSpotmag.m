function spotmagList = ListSpotmag(folderPath,channel,timePoints,thresh1,thresh2)
%{
-About-
To obtain a list of spot magnitude in all cells at certain time points.
This list can be used to calculate SM normalization factor or to draw
histogram of spot intensity at certain timepoint.

Note: my files are named organized as.. 
folderPath/t1-meshaCy5.mat
folderPath/t1-meshaCy3.mat
folderPath/t2-meshaCy5.mat
folderPath/t2-meshaCy3.mat
...
code is based on this format of file names. 

-Inputs-
folderPath: a folder containing spot files (from spotFinder, MicrobeTracker) from one time-course experiment.
channel: 1 for Z5, 2 for Z3.
timePoints: 1, 2,... Usually I name data from t = x min as t(x+1). t = 0 is
t1.

-varargin-

-Outputs-
spotmagList: list of spotmagnitude

-Example-
folderPath = '131101-lacZ 2CFISH\A.MG1655early\'; %last "\" is essential.
spotmagList = ListSpotmag(folderPath,5,[1])
   
-Supplementary-

-Keywords-
spotmag, SM FISH analysis

-Dependencies-
smlacZFISHWorkflow

-References-

-Author-
Sangjin Kim, 2013 November 10
%}

n1 = 0; 
totali = 1; totalf = 1;

for filei = 1:length(timePoints)
    ti = timePoints(filei);
    tifile = strcat('t',sprintf('%01.0f',ti),'-mesha'); 
    cellListfile = strcat(folderPath,tifile,'Cy5.mat'); load(cellListfile); cellList1a = cellList;
    cellListfile = strcat(folderPath,tifile,'Cy3.mat'); load(cellListfile); cellList1b = cellList;
    
    cellindextemp = CellCounter(cellList,ti); %% (ti,frame, cell#)
%     if length(cellindextemp) ~= length(meaninthist(cellList,'signal1'))
%         display('cell indext error!');
%         break;
%     end;
    n1 = length(cellindextemp);

    cellList1a = SpotThresholding(cellList1a,1,thresh1);
    cellList1b = SpotThresholding(cellList1b,2,thresh2);
    
    cellindextemp(:,2) = cellindextemp(:,2) + totalf - 1; %update the frame number
    % combine cellList from different time points.
    cellindex(totali:(totali+n1-1),:) = cellindextemp;
    cellListA(totalf:(totalf+size(cellList1a,2)-1)) = cellList1a;
    cellListB(totalf:(totalf+size(cellList1a,2)-1)) = cellList1b;
    totali = totali+n1;
    totalf = totalf + (size(cellList1a,2));
end
    
%size(cellindex,1) %total # of cells among the selected time points

if channel == 1, 
    spotmagList = SpotmagList(cellListA); 
elseif channel == 2,
    spotmagList = SpotmagList(cellListB);
end

end


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
end


function cellList = SpotThresholding(cellList,ch,thresh)
for frame=1:length(cellList)
    for cell=1:length(cellList{frame})
        if ~isempty(cellList{frame}{cell}) && length(cellList{frame}{cell}.mesh)>4 && isfield(cellList{frame}{cell},'spots') ...
                && ~isempty(cellList{frame}{cell}.spots.x)
            %---for given cell, spot threshold
            if ch == 1
                prf = cellList{frame}{cell}.signal1./cellList{frame}{cell}.steparea;
            elseif ch == 2
                prf = cellList{frame}{cell}.signal2./cellList{frame}{cell}.steparea;
            end;
            
            for i=1:3
                prf = prf*0.5+prf([1 1:end-1])*0.25+prf([2:end end])*0.25; %smoothing
            end
            
            pos = cellList{frame}{cell}.spots.positions;
            prf = [0;prf];
            
            ind = prf(pos+1)>thresh; %1 or 0
            
            ind2 = cellList{frame}{cell}.spots.magnitude>0 & cellList{frame}{cell}.spots.magnitude<5;
            ind2 = ind2';
            try
                ind = ind & ind2; %logic numbers
            catch
            end;
            
            cellList{frame}{cell}.spots.positions = cellList{frame}{cell}.spots.positions(ind); %takes only 1s
            cellList{frame}{cell}.spots.x = cellList{frame}{cell}.spots.x(ind);
            cellList{frame}{cell}.spots.y = cellList{frame}{cell}.spots.y(ind);
            cellList{frame}{cell}.spots.l = cellList{frame}{cell}.spots.l(ind);
            cellList{frame}{cell}.spots.d = cellList{frame}{cell}.spots.d(ind);
            cellList{frame}{cell}.spots.magnitude = cellList{frame}{cell}.spots.magnitude(ind);
            cellList{frame}{cell}.spots.h = cellList{frame}{cell}.spots.h(ind);
            cellList{frame}{cell}.spots.w = cellList{frame}{cell}.spots.w(ind);
            cellList{frame}{cell}.spots.b = cellList{frame}{cell}.spots.b(ind);
        end;
    end;
end;
end

function Xans = SpotmagList(cellList)
spotn = 0;
sc = 0;
cellwithspot =0; 
spotlist = [];
for frame=1:length(cellList)
    for cell=1:length(cellList{frame})
        if ~isempty(cellList{frame}{cell}) && length(cellList{frame}{cell}.mesh)>4 && isfield(cellList{frame}{cell},'spots') 
            flagspot = 0; 
            for spoti = 1:length(cellList{frame}{cell}.spots.magnitude)
                if cellList{frame}{cell}.spots.magnitude(spoti)>0 && cellList{frame}{cell}.spots.magnitude(spoti)<50
                    flagspot = flagspot+1;
                    spotn = spotn + 1; %there is a spot!
                    spotlist(spotn,1) = cellList{frame}{cell}.spots.magnitude(spoti); % spot intensity
                end
            end
        end
    end
end

if ~isempty(spotlist) 
    Xans = spotlist; 
else %incase there is no spot in the cellList.
    Xans = NaN;
end
end