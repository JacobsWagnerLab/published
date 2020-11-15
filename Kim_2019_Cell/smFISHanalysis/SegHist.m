function SegHist(cellList)
%{
-About-
Plot distribution of intensities of segments within cells in cellList
Do for both signal1 and signal2 at the same time
This must be similar to CL_segHist in old sharedRepo database.

-Inputs-
cellList: cellList containing signal1 and signal2

-varargin-

-Outputs-

-Example-
SegHist(cellList)

-Supplementary-

-Keywords-
cellList, segment intensity

-Dependencies-
smlacZFISHworkflow

-References-

-Author-
modified by Sangjin Kim, May 1, 2012
%}

segintarray1 = []; segintarray2 = [];
for frame=1:length(cellList)
    for cell=1:length(cellList{frame})
        if ~isempty(cellList{frame}{cell}) && length(cellList{frame}{cell}.mesh)>4
            segintarray1 = [segintarray1 reshape(cellList{frame}{cell}.signal1./cellList{frame}{cell}.steparea,1,[])];
            segintarray2 = [segintarray2 reshape(cellList{frame}{cell}.signal2./cellList{frame}{cell}.steparea,1,[])];
        end
    end
end

c = 0:(max(segintarray1)/sqrt(length(segintarray1))):max(segintarray1);
h1 = hist(segintarray1,c);
figure, bar(c,h1')
xlabel('Mean intensity in a segment (ch1)')
ylabel('%');

c = 0:(max(segintarray2)/sqrt(length(segintarray2))):max(segintarray2);
h2 = hist(segintarray2,c);
figure, bar(c,h2')
xlabel('Mean intensity in a segment (ch2)')
ylabel('%');



