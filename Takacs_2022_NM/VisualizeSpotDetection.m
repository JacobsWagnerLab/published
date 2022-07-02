%{
-About-
This function allows you to visualize the output from
Find_Irregular_Spots.m to determine the quality of spot detection.

-Inputs-

    cellList -          cellList output from Find_Irregular_Spots.m

    phase -             A cell array of phase images, such that phase{1}
                        corresponds to the frame in cellList.meshData{1}.

    fluor -             Cell array of the fluorescence images, such that
                        fluor{1} corresponds to the frame in
                        cellList.meshData{1}. These images should be those
                        that contain the fluorescence data used in spot
                        detection.

    number_frames -     User can determine the number of frames to scan
                        through. Generally, I'd check 5.

-varargin-


-Outputs-

    figure -            Figure that you can scan through with the phase and
                        fluor images for each cell. It also will show a
                        second copy of the fluor image with the assigned
                        spots from Find_Irregular_Spots.m.

-Example-

   
-Supplementary-


-Keywords-
spot detection, plotting

-Dependencies-

-References-

-Author-
Sander Govers (original), 2017 
Modified to its current form by Molly Scott (2018)
%}
    


% 
%=========================================================================
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
function VisualizeSpotDetection(cellList,phase,fluor,number_frames)

%Create a figure to display your subplots
figure(1);
set(gcf,'position',[38   197   1000   800]);
%Scan through each frame by clicking
for ii=1:number_frames
    img = fluor(:,:,ii); 
    for jj=1:length(cellList.meshData{ii})
    if ~isempty(cellList.meshData{ii}{jj})
        if cellList.meshData{ii}{jj}.mesh ~= 0;
            subplot(1,3,1);
            hold off;
            %establish the mesh and cell that you will display
            M=cellList.meshData{ii}{jj}.mesh;
            mBox=cellList.meshData{ii}{jj}.box;
            X=[M(:,1);flipud(M(:,3))]-mBox(1)+1;
            Y=[M(:,2);flipud(M(:,4))]-mBox(2)+1;
            imCell=phase(mBox(2):mBox(2)+mBox(4),mBox(1):mBox(1)+mBox(3),ii);
            imdb=double(imCell);
            %display that cell with phase, fluor images
            imshow(imdb,[],'initialMagnification','fit');hold on;
            plot(X,Y,'-g');
            title(['Frame = ',num2str(ii),', Cell = ', num2str(jj)]);
            subplot(1,3,2)
            imCellS=fluor(mBox(2):mBox(2)+mBox(4),mBox(1):mBox(1)+mBox(3),ii);
            imdSd=double(imCellS);
            imshow(imdSd,[],'initialMagnification','fit');hold off;        
            subplot(1,3,3)
            imCellS=fluor(mBox(2):mBox(2)+mBox(4),mBox(1):mBox(1)+mBox(3),ii);
            imdSd=double(imCellS);
            %overlay spots onto the image
            imshow(imdSd,[],'initialMagnification','fit');hold on;
            if isfield(cellList.meshData{ii}{jj},'spotPosition')
                for ss=1:size(cellList.meshData{ii}{jj}.spotPosition,1)
                spotPositionX=cellList.meshData{ii}{jj}.spotPosition(ss,1)-mBox(1)+1;
                spotPositionY=cellList.meshData{ii}{jj}.spotPosition(ss,2)-mBox(2)+1;
                spotPositionX = double(spotPositionX);
                spotPositionY = double(spotPositionY);
                plot(spotPositionX,spotPositionY,'bo');
                text(spotPositionX + 2, spotPositionY + 2,num2str(ss),'Color','r')
                end
            end
            hold off;
            pause;
        end
    end
    end
end
end

        