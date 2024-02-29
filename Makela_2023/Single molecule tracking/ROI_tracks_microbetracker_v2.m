function [movie_data] = ROI_tracks_microbetracker_v2(params,cellList_allfull, dataFilename, brightfield_filename)



appData2.dataFilename = dataFilename;
appData2.dataPathname = params.dataPathname;


appData2.data = importdata(str2mat(strcat(appData2.dataPathname, appData2.dataFilename)));

% remove showing figure
% figure;
dataIn = imread(brightfield_filename);
averageImage = double(mean(dataIn,3));   
%     imshow(averageImage,[min(min(averageImage)) mean(max(averageImage))], 'InitialMagnification' ,180 );
%     figurehandle = gcf;
%     colormap(gray);


movie_data.brightfield_image = averageImage;


for ii =  1:length(cellList_allfull)
    plgx = [cellList_allfull{ii}.mesh(1:end-1,1);flipud(cellList_allfull{ii}.mesh(:,3))];
    plgy = [cellList_allfull{ii}.mesh(1:end-1,2);flipud(cellList_allfull{ii}.mesh(:,4))];
    [a,b] = expandpoly(plgx, plgy, 0);% changed to zero from 1.2
    ROIs{ii} = [a,b];
    %ROIs{ii} = [plgx, plgy];
    
end


% [transform_filename, transform_pathname] = uigetfile('*.mat', 'Transform: ','MultiSelect', 'off');
% load([ transform_pathname, transform_filename]); 

%data  = tformfwd(T,appData2.data.data(:,2:3));  

data = appData2.data;
  
ROIData2 = selectCellROIs_copynumbers_useROI(data,params,ROIs);
cellROI_data  = ROIData2;


for j = 1:numel(cellROI_data)
    cellROI_data(j).length = cellList_allfull{j}.length;
    cellROI_data(j).lengthvector = cellList_allfull{j}.lengthvector;
    cellROI_data(j).mesh = cellList_allfull{j}.mesh;
    cellROI_data(j).area = cellList_allfull{j}.area;
end

movie_data.cellROI_data = cellROI_data;

end


%%

function ROIData = selectCellROIs_copynumbers_useROI(data,appData,ROIs)


hold all;




%initialise data structure
tempStruct = ...
    struct('localizationData',zeros(1,2),...
    'tracks',zeros(1,2),...
    'nMolecules',zeros(1,1),...
    'ROIVertices',zeros(1,2));

if isfield(data,'data')
remainingData = data.data;
else
remainingData = data;
end

%shift localisations to match brightfield meshes
if appData.pixelshift(1) ~= 0
    % %-2 will shift left by 2 pixels
    remainingData(:,2) = remainingData(:,2)+appData.pixelshift(1);
end

if appData.pixelshift(2) ~= 0
    % %-2 will shift up by 2 pixels
    remainingData(:,3) = remainingData(:,3)+appData.pixelshift(2);
end

ROIData = [];


ii = 1;




for i = 1:numel(ROIs)
    
    % plotting
%     v = axis;
    
%     h = plot(remainingData(:,2),remainingData(:,3),'b.','MarkerSize',2);
%     axis(v);
%     drawnow;
    
    ROIVertices = ROIs{i};
    
    inIndexes = find( inpolygon(remainingData(:,2),remainingData(:,3),ROIVertices(:,1),ROIVertices(:,2)) );
    inData = remainingData(inIndexes,:);
    
    remainingData(inIndexes,:) = [];
    
    inData(:,11) = ii * ones(numel(inIndexes),1); % ROI id
    
    if ~isempty(inData)
    %    track data for this ROI
    pos(:,1) = inData(:,2);
    pos(:,2) = inData(:,3);
    pos(:,3) = inData(:,1);
   
    appData.trackParams.quiet = 1;
    tracks = trackWithDummy(pos, appData.trackParams);
    nMolecules = max(tracks(:,4));
    clear pos;
    
    tempStruct.localizationData = inData;
    tempStruct.tracks = tracks;
    tempStruct.nMolecules = nMolecules;
    tempStruct.ROIVertices = ROIVertices;
    
    ROIData = [ROIData; tempStruct];
    else 
    tempStruct.localizationData = [];
    tempStruct.tracks = [];
    tempStruct.nMolecules = [];
    tempStruct.ROIVertices = ROIVertices;
    
    ROIData = [ROIData; tempStruct];    
    warning('empty cells');
    end
    
%     delete(h);
    
    ii = ii + 1;
    
    if ii == round((numel(ROIs)./4))
        disp('25%');
    elseif ii == round((numel(ROIs)./2))
        disp('50%')
    elseif ii == round((3.*numel(ROIs)./4))
        disp('75%')
    end
        
    
end


% v = axis;

% h = plot(remainingData(:,2),remainingData(:,3),'b.','MarkerSize',1);
% axis(v);

end