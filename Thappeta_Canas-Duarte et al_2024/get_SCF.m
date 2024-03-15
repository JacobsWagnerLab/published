%% Calculating the signal correlation factor for all cells in the cellList
% -Input-
%   cellList: Oufti cellList
%   c1_pth: file path to the 1st fluoresence channel
%   c2_pth: file path to the 2nd fluoresence channel
%   dx_from_center: This determines how many pixels away from the 
%       centerline you want to calculate the correlation. Must add up to an
%       integer. A value of one means that the correlation area will be two
%       pixels wide
%   pole_length: This is the number of pixels away from the pole that you 
%       want to calculate the correlation for. This is to avoid the 
%       artificial positive correlation at the pole that results from the 
%       decrease in volume.
%   field_name: the name of the field the SCF results should be stored to
%   in the cell structure. This is useful, because you may want to
%   calculate the SCF between various pairs of channels. Explicitly
%   indicate the field name, so that you can distinguish the different SCF
%   results based on different pairs of channels.
%
% -Output-
%   cellList: Oufti cellList with the addtion of SCF results
%
% -Author-
%   Yingjie Xiang, 12/30/2020

function cellList = get_SCF(cellList,c1_pth,c2_pth,dx_from_center,pole_length,field_name)
c1_im = loadimseries(c1_pth);
c2_im = loadimseries(c2_pth);

for frame = 1:length(cellList.meshData)
    for cc = 1:length(cellList.meshData{frame})
        try
            [corr_matrix,~] = Cell_Pixel_Correlation({double(c1_im(:,:,frame)),double(c2_im(:,:,frame))},cellList,frame,cc,dx_from_center,pole_length);
            scf = corr_matrix(2,1);
            cellList.meshData{frame}{cc}.(field_name) = scf;
        catch
            fprintf('Failed for Frame #%d, Cell #%d\n',frame,cc);
            continue;
        end
    end
end
end