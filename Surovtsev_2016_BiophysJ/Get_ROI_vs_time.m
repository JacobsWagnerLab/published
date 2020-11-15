
% Intitially to get FRAP curve for ROI in specified cell

cell=1;
r=2.5;
x0=125; dx0=5; %x0=73; dx0=5;
dims=[x0-dx0,x0+dx0,r];

IDs=cellList.cellId;
Meshes=cellList.meshData;

id=find(IDs{1}==cell);
n_frames=length(IDs);

signal0=sum(double(Meshes{1}{id}.signal1));
for frame=1:n_frames
  
  cellMesh=Meshes{frame}{id}.mesh;
  cellImage=TL_Fluo(:,:,frame);
  
  [roi_indx, px_values,roi] = extractCellPixels_2(cellImage,cellMesh, dims);
  signal1=sum(double(Meshes{frame}{id}.signal1));
  
  ROI_vs_t(frame,:)=[length(roi_indx),mean(px_values),median(px_values),min(px_values),max(px_values),signal1];
  
    
end
