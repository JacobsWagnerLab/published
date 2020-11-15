
%cellList_0=cellList; cellList=new2old_cellList(cellList); % to convert new cell List into old one
%cellList=process_cellList(cellList);
cellList2=peakfinder(cellList,.2,.2,[],[],[],[],'signal1',5);

% frame=1; Cell=8;vector=cellList2{frame}{Cell}.lengthvector;
% fSignal=cellList2{frame}{Cell}.signal1;
% 
% spotpos=cellList2{frame}{Cell}.spots.l;
% figure, plot(vector, fSignal), hold on, 
% plot(spotpos ,repmat(0.05,1, length(spotpos)),'o');

% data=zeros(219,1);
% count=0;
% for frame=1:length(cellList)
%      for Cell=1:length(cellList{frame})
%          count=count+1;
%          data(count)=cellList{frame}{Cell}.length*64;
%      end
%  end
% figure, hist(data)

adjust_peakfinder;% add 1 pixel to the MipZ spots to correct for offset;

%data=parAtotbetMipZ(cellList2);
[data,YY,YY_all]=parAtotbetMipZ_2(cellList2);
% data column designation: parA1, dist, lng, spotNum
% parA1: parA fluorescence within each interval
% dist: distance of the interval
% lng: the length of the corresponding cell
% spotNum: the number of spots within the corresponding cell

