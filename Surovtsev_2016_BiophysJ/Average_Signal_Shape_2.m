% to get average shape of the signal (ParA)
% small modiciation to get 3 pixels before the PC


% orient cellList according to signal2
cellList = CL_orientCellList(cellList,1,'signal2',0,1,0);
 %cellList = peakfinder_(cellList,'signal2');

n=50; % sampling over [0,1] interval

MeshData=cellList.meshData;
 IDData=cellList.cellId;

n_frames=length(MeshData);
S1_av=[]; S2_av=[];
n_cells=0; 
for frame=1:n_frames
  MeshData_1=MeshData{frame};  
  for cell=1:length(MeshData_1)
    CellData=MeshData_1{cell};
    if isfield(CellData,'spots') && ~isempty(CellData.spots) && length(CellData.spots.l)==1 && ~isnan(CellData.spots.positions);
      n_cells=n_cells+1;  
%       signal1=CellData.signal1;
%        signal2=CellData.signal2;
      signal1=CellData.signal1./CellData.steparea;
       signal2=CellData.signal2./CellData.steparea;
      
      ind=CellData.spots.positions;
      if ind>2
        ind=ind-2;  
        else
        break  
      end
      ind_fin=length(signal1)-3;   
      
      
      l=CellData.length-CellData.lengthvector(ind);
       lengthvector=CellData.lengthvector(ind:ind_fin)-CellData.lengthvector(ind);
            
%       y1=integinterp(lengthvector/l,signal1(ind:ind_fin),n); 
%        S1_all(n_cells,:)=y1/sum(y1);
%        y2=integinterp(lengthvector/l,signal2(ind:ind_fin),n); 
%        S2_all(n_cells,:)=y2/sum(y2);
      y1=integinterp(lengthvector/l,signal1(ind:ind_fin),n); 
       S1_all(n_cells,:)=y1*(ind_fin-ind+1)/sum(y1); % similar to y_norm=y*length(parA1_dist)/sum(y)
      y2=integinterp(lengthvector/l,signal2(ind:ind_fin),n); 
       S2_all(n_cells,:)=y2*(ind_fin-ind+1)/sum(y2);
    end     
    
  end
    
end

S1_av=mean(S1_all,1);
 S2_av=mean(S2_all,1);
 
 figure; hold on
 plot(0:1/(n-1):1,S1_av,'-r','LineWidth',2)
  plot(0:1/(n-1):1,S2_av,'-g','LineWidth',2)
 xlabel('Relative cell coordinate','FontSize',18)
  ylabel('Fluorescence','FontSize',18)
  set(gca,'FontSize',16)
  legend({'Signal1','Signal2'},'FontSize',14)