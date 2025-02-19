function SVMCuratedCellDetection(oufti_analysis,organism,pth)
% author:     Sander K Govers
% date:       02/20/2017
% copyright:  Yale University
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SVMCuratedCellDetection cleans up cell meshes detected by Oufti based on
% pre-trained SVM models. Cells can either be E. coli or C. crescentus
% (since different SVM models have been trained for both organisms).
% Function first generates a set of predictors (different for each
% organism) which are subsequently used by the SVM model to classify a
% given cell detection as 'good' or 'bad'. The latter type of cell
% detections is then removed from the Oufti analysis and the entire thing is saved as a separate analysis file (_SVMcurated) in the same folder.

%The predictors used to classify cell detections are generated based on the
%properties of the cell meshes in combination with the respective phase
%contrast images.

%Both SVM models were constructed using MATLAB's fitcsvm function in which
%a linear kernel together with a cost matrix (to reduce the number of false
%positives) were employed. Performace details (cross-validation, AUROC, confusion matrix) of both SVM models can be
%found in a separate pdf file accompanying this function.

%IMPORTANT: The SVM classifiers might be (aka are) sensitive to the parameter set used during Oufti analysis. Therefore, one
%should use the specified parameter sets (location) for optimal
%performance. They can be found in the same folder that contains the SVM
%models (\\aunt\common\sharedRepo\cellList_curation\). 
%Only the threshold factor should be adapted upon using the specified
%parameter sets.
%For E. coli, also adapt the cell width parameter (depending on the size of the cells or the
%nutrient condition they were grown in).
%Although the predictors used for SVM classification should in theory be
%size-independent, the SVM models themselves were trained on WT cells and
%as such will only perform optimal on WT-looking cells. Mutants displaying
%aberrant cell morphologies might be incorrectly identified as 'bad' cell
%detections.

%IMPORTANT2: You'll need Matlab 2014 or more recent given that this type of SVM
%formulation was only introduced in Matlab R2014a

%INPUT
%-oufti_analysis    Analysis file obtained from Oufti (do not append .mat)
%-organism         Organism used. Either 'EC' or 'CC', for E. coli and
%                  C. crescentus respectively.
%-pth              Path to folder containing oufti analysis file and phase contrast images (the latter should be placed in a separate folder named 'c1' in the path folder). 

%OUPUT
%-ouftiAnalysis_SVMcurated    Analysis file in which cell detections classified as
%                             'bad' by the SVM models are removed. File is automatically saved in folder
%                             containing the original Oufti analysis file. This file can then be
%                             loaded back into Oufti for visual inspection of the
%                             curation procedure. 

if nargin == 3
    if strcmp(organism, 'EC')
        do_Ecoli=true;
        do_Caulobacter=false;
    elseif strcmp(organism, 'CC')
        do_Ecoli=false;
        do_Caulobacter=true;
    else 
        disp('Please specify correct type of organism');
    end
else
    disp('Input parameter missing')
end

%Identify all phase images
files = dir([pth '\c1' '/*.tif*']);

%Adjust pixel size accordingly
pixel_size=6.4/100;

%load cellList from oufti
oufti_analysis_file=load([pth,'\',oufti_analysis]);
cellList=oufti_analysis_file.cellList;
tic
%% Clean out cell entries ahead of time to avoid too much 'if' testing in parallel loop
kk=0;
for ii=1:length(cellList.meshData)
    for jj=1:length(cellList.meshData{ii})
        if isempty(cellList.meshData{ii}{jj}) ||...
            ~isfield(cellList.meshData{ii}{jj},'mesh') ||...
            length(cellList.meshData{ii}{jj}.mesh)<=4
            cellList.meshData{ii}{jj}=[];
            kk=kk+1;
        end
    end
end
disp(['Number of cell entries cleaned out : ',num2str(kk)]);

%% E. coli part
if do_Ecoli
    w = waitbar(0,'Generating E. coli predictors...');
    %Generate structure containing predictor info
    tab = struct('perimeter',[],'phase_intensity',[],'contour_intensity',[],'min_contour_intensity', [], 'max_curvature', [], 'min_curvature', [],...
    'expanded_contour_intensity',[],'ratio_expanded_contour_intensity',[],'max_contour_intensity_variability',[],'max_expanded_contour_intensity',[],'frame',[],'cell',[]);
    
    for ii=1:length(cellList.meshData)
        clear meshData tmpCell image
        tmp_structure = struct('perimeter',[],'phase_intensity',[],'contour_intensity',[],'min_contour_intensity', [], 'max_curvature', [], 'min_curvature', [],...
        'expanded_contour_intensity',[],'ratio_expanded_contour_intensity',[],'max_contour_intensity_variability',[],'max_expanded_contour_intensity',[],'frame',[],'cell',[]);
        
        if ~isempty(cellList.meshData{ii})
            %Go through all frames
            meshData = cellList.meshData{ii};
            image=imread([pth '\c1' '/' files(ii).name]);
                for jj=1:length(meshData)
                if ~isempty(cellList.meshData{ii}{jj})
                    %Process all cells within frame ii
                    tmp_structure(jj).frame=ii;
                    tmp_structure(jj).cell=jj;
                    
                    mesh=double(meshData{jj}.mesh);
                    
                    %Separate mesh features for perimeter calculation
                    x_left_mesh=mesh(:,1);
                    y_left_mesh=mesh(:,2);
                    x_right_mesh=mesh(:,3);
                    y_right_mesh=mesh(:,4);

                 %Calculate perimeter
                    coordinates_mesh = [x_left_mesh y_left_mesh;flipud([x_right_mesh y_right_mesh])];
                    distance_between_seperate_coordinates=coordinates_mesh(1:end-1,:)-coordinates_mesh(2:end,:);
                    perimeter=bsxfun(@hypot,distance_between_seperate_coordinates(:,1),distance_between_seperate_coordinates(:,2));
                    tmp_structure(jj).perimeter = sum(perimeter).*pixel_size;
                    
                 %Calculate phase intensity 
                    %Cell with predictor rows containing NaN, indicating
                    %something went wrong during predictor generation
                    %will be automatically removed at the end
                    tmp_structure(jj).phase_intensity = NaN;
                    if isfield(meshData{jj},'signal0') && ~isempty(meshData{jj}.signal0)
                        %Get total phase contrast intensity  
                        tmp_structure(jj).phase_intensity = sum(meshData{jj}.signal0);
                    end
                    
                 %Calculate cell contour intensity
                    %Get cell box amd corresponding cropped image 
                    cell_box=double(meshData{jj}.box);
                    image_cell=image(cell_box(2):cell_box(2)+cell_box(4),cell_box(1):cell_box(1)+cell_box(3));
                    image_cell_double=double(image_cell);
                    %Cell mesh coordinates
                    X_cell_contour=double([mesh(:,1);flipud(mesh(:,3))])-cell_box(1)+1;
                    Y_cell_contour=double([mesh(:,2);flipud(mesh(:,4))])-cell_box(2)+1;
                    %Improfile function (built-in matlab function) to
                    %calculate intensity values of pixels along a cell
                    %contour
                    contour_intensity=improfile(image_cell,X_cell_contour,Y_cell_contour);
                    %Normalize contour intensity values
                    contour_intensity=(contour_intensity-min(image_cell_double(:)))./(max(image_cell_double(:))-min(image_cell_double(:)));
                    %Average contour intensity
                    tmp_structure(jj).contour_intensity=mean(contour_intensity);
                    
                    %Sort contour intensity values to be able te extract
                    %minimum values
                    sort_contour_intensity=sort(improfile(image_cell,X_cell_contour,Y_cell_contour), 'descend');
                    %Calculate minimum contour intensity values based on
                    %lowest 10 values
                    if length(sort_contour_intensity)>10
                        tmp_structure(jj).min_contour_intensity=mean(sort_contour_intensity(end-10:end));
                    else
                        tmp_structure(jj).min_contour_intensity=NaN;
                    end

                    %Calculate variance along contour intensities by
                    %comparing one point with another one 10 steps further
                    extended_contour_intensity=vertcat(contour_intensity(1:end),contour_intensity(1),contour_intensity(2),contour_intensity(3),contour_intensity(4),contour_intensity(5),contour_intensity(6),contour_intensity(7),contour_intensity(8),contour_intensity(9),contour_intensity(10));
                    %Generate vector containing contour intensity variances
                    contour_intensity_variability =[];
                    for ss=2:length(extended_contour_intensity)-10
                        ten_contour_points_further = std([extended_contour_intensity(ss) extended_contour_intensity(ss+10)]);
                        contour_intensity_variability = [contour_intensity_variability;ten_contour_points_further];
                    end
                    %Maximum variance of contour intensities
                    tmp_structure(jj).max_contour_intensity_variability=max(contour_intensity_variability);
                    
                 %Calculate contour curvature
                    model=double(meshData{jj}.model);
                    model_curvature = LineCurvature2D(model);
                    tmp_structure(jj).max_curvature=max(model_curvature);
                    tmp_structure(jj).min_curvature=min(model_curvature);

                 %Calculate expanded cell polygon and extract its
                 %properties
                    %Create expanded cell polygon       
                    mask = poly2mask(X_cell_contour,Y_cell_contour,cell_box(:,4),cell_box(:,3));
                    mask_dilated = imdilate(mask,strel('square',3));
                    %Get coordinates of expanded cell contour
                    model_dilated = cell2mat(bwboundaries(mask_dilated,'noholes'));
                    
                    if ~isempty(model_dilated)
                        %X and Y coordinates of expanded cell contour (X
                        %and Y are switched after using cell2mat)
                        Y_dilated_cell_contour=double(model_dilated(:,1));
                        X_dilated_cell_contour=double(model_dilated(:,2));
                        %Calculate pixel intensity along cell contour
                        contour_intensity_dilated=improfile(image_cell,X_dilated_cell_contour,Y_dilated_cell_contour);
                        %Normalize
                        contour_intensity_dilated=(contour_intensity_dilated-min(image_cell_double(:)))./(max(image_cell_double(:))-min(image_cell_double(:)));
                        %Average expanded contour intensity
                        tmp_structure(jj).expanded_contour_intensity=mean(contour_intensity_dilated);
                        %Maximum expanded contour intensity
                        tmp_structure(jj).max_expanded_contour_intensity=max(contour_intensity_dilated);
                        %Compare two halves of expaned contour
                        half_index = round(length(contour_intensity_dilated)/2);
                        %Ratio of halve's intensities
                        tmp_structure(jj).ratio_expanded_contour_intensity = sum(contour_intensity_dilated(1:half_index))/sum(contour_intensity_dilated(half_index+1:end));
                    else
                        tmp_structure(jj).expanded_contour_intensity=NaN;
                        tmp_structure(jj).max_expanded_contour_intensity=NaN;
                        tmp_structure(jj).ratio_expanded_contour_intensity=NaN;
                    end
                end
                end
        end
        tab=[tab,tmp_structure];

        waitbar(ii/length(cellList.meshData));
    end
    tab(1)=[];
    
    clear meshData img mesh ii jj kk x1 x2 y1 y2 a b perim
    clear mBox imCell imdb X Y contInt sortcontInt extcontInt contIntVar
    clear model modelCurv mask maskDil modelDil YY XX contIntDil ixExp
    clear cc dd ss tenContPoints
    
    table_array(:,1)=cat(1,tab.perimeter);
    table_array(:,2)=cat(1,tab.phase_intensity);
    table_array(:,3)=cat(1,tab.contour_intensity);
    table_array(:,4)=cat(1,tab.min_contour_intensity);
    table_array(:,5)=cat(1,tab.max_curvature);
    table_array(:,6)=cat(1,tab.min_curvature);
    table_array(:,7)=cat(1,tab.expanded_contour_intensity);
    table_array(:,8)=cat(1,tab.ratio_expanded_contour_intensity);
    table_array(:,9)=cat(1,tab.max_contour_intensity_variability);
    table_array(:,10)=cat(1,tab.max_expanded_contour_intensity);
    table_array(:,11)=cat(1,tab.frame); 
    table_array(:,12)=cat(1,tab.cell);     
    
    out=table_array;
    
    %Create array with normalized values. Normalization is performed by substracting the median and dividing by mad (median absolute deviation)
    %Non-normal distributions are first log-transformed

    normalized_array = zeros(length(out), 12);
    %perimeter
    normalized_array(:,1)=cat(1,(log(out(:,1))-nanmedian(log(out(:,1))))/mad(log(out(:,1)),1));
    %phase intensity
    normalized_array(:,2)=cat(1,(out(:,2)-nanmedian(out(:,2)))/mad(out(:,2),1));
    %contour intensity
    normalized_array(:,3)=cat(1,(log(out(:,3))-nanmedian(log(out(:,3))))/mad(log(out(:,3)),1));
    %minimum contour intensity
    normalized_array(:,4)=cat(1,(log(out(:,4))-nanmedian(log(out(:,4))))/mad(log(out(:,4)),1));
    %maximum curvature
    normalized_array(:,5)=cat(1,(log(out(:,5))-nanmedian(log(out(:,5))))/mad(log(out(:,5)),1));
    %minimum curvature
    normalized_array(:,6)=cat(1,(log(out(:,6))-nanmedian(log(out(:,6))))/mad(log(out(:,6)),1));
    %expanded contour intensity
    normalized_array(:,7)=cat(1,(log(out(:,7))-nanmedian(log(out(:,7))))/mad(log(out(:,7)),1));
    %ratio expanded contour intensity
    normalized_array(:,8)=cat(1,(log(out(:,8))-nanmedian(log(out(:,8))))/mad(log(out(:,8)),1));
    %maximum variance of contour intensity
    normalized_array(:,9)=cat(1,(log(out(:,9))-nanmedian(log(out(:,9))))/mad(log(out(:,9)),1));
    %maximum expanded contour intensity
    normalized_array(:,10)=cat(1,(log(out(:,10))-nanmedian(log(out(:,10))))/mad(log(out(:,10)),1));
    %frame
    normalized_array(:,11)=cat(1,out(:,11));
    %celll
    normalized_array(:,12)=cat(1,out(:,12));
    
    predictor_array=real(normalized_array(:,1:10)); 
    
    %Check for NaN values, to set these rows artificially to 0 after SVM prediction
    delete_row1=any(isnan(predictor_array),2);
    delete_row2=any(isinf(predictor_array),2);
    delete_row=bsxfun(@or,delete_row1,delete_row2);
close (w);
end



%% C. crescentus part
if do_Caulobacter
    wa = waitbar(0,'Generating C. crescentus predictors...');
    %Generate structure containing predictor info
    tab = struct('volume',[],'contour_intensity',[],'max_contour_intensity',[],'min_contour_intensity', [], 'max_curvature', [], 'min_curvature', [],...
    'min_expanded_contour_intensity',[],'mesh_gradient',[],'max_contour_intensity_variability',[],'max_phase_intensity',[],'midcell_phase_intensity_kurtosis',[],'frame',[],'cell',[]);

    for ii=1:length(cellList.meshData)
        clear meshData tmpCell
        tmp_structure = struct('volume',[],'contour_intensity',[],'max_contour_intensity',[],'min_contour_intensity', [], 'max_curvature', [], 'min_curvature', [],...
    'min_expanded_contour_intensity',[],'mesh_gradient',[],'max_contour_intensity_variability',[],'max_phase_intensity',[],'midcell_phase_intensity_kurtosis',[],'frame',[],'cell',[]);

        if ~isempty(cellList.meshData{ii})
            %Go through all frames
            meshData = cellList.meshData{ii};
            image=imread([pth '\c1' '/' files(ii).name]);
            %Calculation of mesh gradient requires normalized images
            image_sorted_pixels = sort(image(:));
            pixel_number = length(image_sorted_pixels);
            min_intensity = mean(image_sorted_pixels(1:round(0.02*pixel_number)));
            max_intensity = mean(image_sorted_pixels(round(0.98*pixel_number:end))); 
            image_normalized = (double(image)-min_intensity)./(max_intensity-min_intensity);

            for jj=1:length(meshData)
            if ~isempty(cellList.meshData{ii}{jj})
                %Process all cells within frame ii
                tmp_structure(jj).frame=ii;
                tmp_structure(jj).cell=jj;
                
                mesh=double(meshData{jj}.mesh);
             %Calculate cell volume
                %Get steplength based on mesh in order to be able to
                %calculate cell volume
                steplength = edist(mesh(2:end,1)+mesh(2:end,3),mesh(2:end,2)+mesh(2:end,4),...
                mesh(1:end-1,1)+mesh(1:end-1,3),mesh(1:end-1,2)+mesh(1:end-1,4))/2;
                %Calculate cell volume
                d = edist(mesh(:,1),mesh(:,2),mesh(:,3),mesh(:,4));
                stepvolume = (d(1:end-1).*d(2:end) + (d(1:end-1)-d(2:end)).^2/3).*steplength*pi/4;
                tmp_structure(jj).volume = sum(stepvolume)*pixel_size*pixel_size*pixel_size;

             %Calculate cell contour intensity
                %Get cell box amd corresponding cropped image 
                cell_box=double(meshData{jj}.box);
                image_cell=image(cell_box(2):cell_box(2)+cell_box(4),cell_box(1):cell_box(1)+cell_box(3));
                image_cell_double=double(image_cell);
                %Cell mesh coordinates
                X_cell_contour=double([mesh(:,1);flipud(mesh(:,3))])-cell_box(1)+1;
                Y_cell_contour=double([mesh(:,2);flipud(mesh(:,4))])-cell_box(2)+1;
                %Improfile function (built-in matlab function) to
                %calculate intensity values of pixels along a cell
                %contour
                contour_intensity=improfile(image_cell,X_cell_contour,Y_cell_contour);
                %Normalize contour intensity values
                contour_intensity=(contour_intensity-min(image_cell_double(:)))./(max(image_cell_double(:))-min(image_cell_double(:)));
                 %Average contour intensity
                tmp_structure(jj).contour_intensity=mean(contour_intensity);

                %Sort contour intensity values to be able te extract
                %maximum and minimum values
                sort_contour_intensity=sort(improfile(image_cell,X_cell_contour,Y_cell_contour), 'descend');
                %Calculate maximum contour intensity values
                tmp_structure(jj).max_contour_intensity=mean(sort_contour_intensity(1));
                %Calculate minimum contour intensity values based on
                %lowest 10 values
                if length(sort_contour_intensity)>10
                    tmp_structure(jj).min_contour_intensity=mean(sort_contour_intensity(end-10:end));
                else
                    tmp_structure(jj).min_contour_intensity=NaN;
                end

                %Calculate variance along contour intensities by
                %comparing one point with another one 10 steps further
                extended_contour_intensity=vertcat(contour_intensity(1:end),contour_intensity(1),contour_intensity(2),contour_intensity(3),contour_intensity(4),contour_intensity(5),contour_intensity(6),contour_intensity(7),contour_intensity(8),contour_intensity(9),contour_intensity(10));
                %Generate vector containing contour intensity variances
                contour_intensity_variability =[];
                for ss=2:length(extended_contour_intensity)-10
                    ten_contour_points_further = std([extended_contour_intensity(ss) extended_contour_intensity(ss+10)]);
                    contour_intensity_variability = [contour_intensity_variability;ten_contour_points_further];
                end
                %Maximum variance of contour intensities
                tmp_structure(jj).max_contour_intensity_variability=max(contour_intensity_variability);

             %Calculate contour curvature
                model=double(meshData{jj}.model);
                model_curvature = LineCurvature2D(model);
                tmp_structure(jj).max_curvature=max(model_curvature);
                tmp_structure(jj).min_curvature=min(model_curvature);

             %Calculate expanded cell polygon and extract its
             %minimum contour intensity
                %Create expanded cell polygon       
                mask = poly2mask(X_cell_contour,Y_cell_contour,cell_box(:,4),cell_box(:,3));
                mask_dilated = imdilate(mask,strel('square',3));
                %Get coordinates of expanded cell contour
                model_dilated = cell2mat(bwboundaries(mask_dilated,'noholes'));

                if ~isempty(model_dilated)
                    %X and Y coordinates of expanded cell contour (X
                    %and Y are switched after using cell2mat)
                    Y_dilated_cell_contour=double(model_dilated(:,1));
                    X_dilated_cell_contour=double(model_dilated(:,2));
                    %Calculate pixel intensity along cell contour
                    contour_intensity_dilated=improfile(image_cell,X_dilated_cell_contour,Y_dilated_cell_contour);
                    %Normalize
                    contour_intensity_dilated=(contour_intensity_dilated-min(image_cell_double(:)))./(max(image_cell_double(:))-min(image_cell_double(:)));
                    %Minimum expanded contour intensity
                    tmp_structure(jj).min_expanded_contour_intensity=min(contour_intensity_dilated);
                else
                    tmp_structure(jj).min_expanded_contour_intensity=NaN;
                end

             %Calculate mesh gradient
                %using Brad's cellEdgeGradient function)
                try
                    grad = cellEdgeGradient(mesh, 1, image_normalized);
                catch
                    grad= NaN;
                end
                tmp_structure(jj).mesh_gradient = grad;

             %Get maximum phase intensity
                pixel_values = double(image_cell);
                [I,J] = meshgrid(1:size(image_cell,2), 1:size(image_cell,1));
                %Determine which pixels of cropped image fall within
                %cell contour
                IN = inpolygon(I,J,X_cell_contour,Y_cell_contour);
                pixel_values_cell = pixel_values(IN);
                %Get the maximum pixel value of pixels within the cell
                %contour
                if ~isempty(pixel_values_cell)
                    tmp_structure(jj).max_phase_intensity = max(pixel_values_cell);
                else
                    tmp_structure(jj).max_phase_intensity = NaN;
                end

             %Get midcell region pixel intensities
                %Midcell is taken here as the middle 30% of the cell
                ix35=round(length(mesh)*0.35);
                ix65=round(length(mesh)*0.65);
                midcell= mesh(ix35:ix65,:);
                %X and Y coordinates of polygon for midcell region
                X_midcell_contour=double([midcell(:,1);flipud(midcell(:,3));midcell(1,1)])-cell_box(1)+1;
                Y_midcell_contour=double([midcell(:,2);flipud(midcell(:,4));midcell(1,2)])-cell_box(2)+1;
                %Identify pixels in midcell region
                pixels_in_midcell = inpolygon(I,J,X_midcell_contour,Y_midcell_contour);
                pixel_values_midcell = pixel_values(pixels_in_midcell);
                %The kurtosis (fourth moment) of midcell pixel values is
                %used here because of its potential to identify the
                %"tailedness" of a given distribution
                if ~isempty(pixel_values_midcell)
                    tmp_structure(jj).midcell_phase_intensity_kurtosis = kurtosis(pixel_values_midcell);
                else
                    tmp_structure(jj).midcell_phase_intensity_kurtosis =NaN;
                end
            end
            end
        end
    tab=[tab,tmp_structure];

    waitbar(ii/length(cellList.meshData));
    end
    tab(1)=[];
    
    clear meshData img imgSort ss maxInt minInt imgNorm mesh steplength ii jj kk x1 x2 y1 y2 a b perim
    clear d stepvolume grad ix35 ix65 midCell XmidCell YmidCell INmidCell pxMidCell
    clear mBox imCell imdb X Y contInt sortcontInt extcontInt contIntVar
    clear model modelCurv mask maskDil modelDil YY XX contIntDil ixExp
    
    table_array(:,1)=cat(1,tab.volume);
    table_array(:,2)=cat(1,tab.contour_intensity);
    table_array(:,3)=cat(1,tab.max_contour_intensity);
    table_array(:,4)=cat(1,tab.min_contour_intensity);
    table_array(:,5)=cat(1,tab.max_curvature);
    table_array(:,6)=cat(1,tab.min_curvature);
    table_array(:,7)=cat(1,tab.min_expanded_contour_intensity);
    table_array(:,8)=cat(1,tab.mesh_gradient);
    table_array(:,9)=cat(1,tab.max_contour_intensity_variability);
    table_array(:,10)=cat(1,tab.max_phase_intensity);
    table_array(:,11)=cat(1,tab.midcell_phase_intensity_kurtosis);    
    table_array(:,12)=cat(1,tab.frame); 
    table_array(:,13)=cat(1,tab.cell); 
    
    out=table_array;
    
    %Create array with normalized values. Normalization is performed by substracting the median and dividing by mad (median absolute deviation)
    %Non-normal distributions are first log-transformed

    normalized_array = zeros(length(out), 13);
    %volume
    normalized_array(:,1)=cat(1,(log(out(:,1))-nanmedian(log(out(:,1))))/mad(log(out(:,1)),1));
    %contour intensity
    normalized_array(:,2)=cat(1,(out(:,2)-nanmedian(out(:,2)))/mad(out(:,2),1));
    %maximum contour intensity
    normalized_array(:,3)=cat(1,(log(out(:,3))-nanmedian(log(out(:,3))))/mad(log(out(:,3)),1));
    %minimum contour intensity
    normalized_array(:,4)=cat(1,(log(out(:,4))-nanmedian(log(out(:,4))))/mad(log(out(:,4)),1));
    %maximum curvature
    normalized_array(:,5)=cat(1,(log(out(:,5))-nanmedian(log(out(:,5))))/mad(log(out(:,5)),1));
    %minimum curvature
    normalized_array(:,6)=cat(1,(log(out(:,6))-nanmedian(log(out(:,6))))/mad(log(out(:,6)),1));
    %minimum expanded contour intensity
    normalized_array(:,7)=cat(1,(log(out(:,7))-nanmedian(log(out(:,7))))/mad(log(out(:,7)),1));
    %mesh gradient
    normalized_array(:,8)=cat(1,(log(out(:,8))-nanmedian(log(out(:,8))))/mad(log(out(:,8)),1));
    %maximum variance of contour intensity
    normalized_array(:,9)=cat(1,(log(out(:,9))-nanmedian(log(out(:,9))))/mad(log(out(:,9)),1));
    %maximum phase intensity
    normalized_array(:,10)=cat(1,(log(out(:,10))-nanmedian(log(out(:,10))))/mad(log(out(:,10)),1));
    %kurtosis of phase intensity at midcell
    normalized_array(:,11)=cat(1,(log(out(:,11))-nanmedian(log(out(:,11))))/mad(log(out(:,11)),1));    
    %frame
    normalized_array(:,12)=cat(1,out(:,12));
    %celll
    normalized_array(:,13)=cat(1,out(:,13));
    
    predictor_array=real(normalized_array(:,1:11)); 
    
    %Check for NaN values, to set these rows artificially to 0 after SVM prediction
    delete_row1=any(isnan(predictor_array),2);
    delete_row2=any(isinf(predictor_array),2);
    delete_row=bsxfun(@or,delete_row1,delete_row2);
    close (wa);
end

toc


%% Perform SVM predictions
wb = waitbar(0,'Classifying and saving...');
if do_Ecoli
    SVMModel = load('\\aunt\common\sharedRepo\cellList_curation\SVMModel_EC.mat','SVMModel');
    label = predict(SVMModel.SVMModel,predictor_array);
    
    %Identify incorrect/bad cell detections
    missCells=[];
    for ii=1:length(out)
        if label(ii)==0 || delete_row(ii)==1
            missCells=[missCells; normalized_array(ii,11) normalized_array(ii,12)];
        end
    end
end

if do_Caulobacter
    SVMModel = load('\\aunt\common\sharedRepo\cellList_curation\SVMModel_CC.mat','SVMModel');
    label = predict(SVMModel.SVMModel,predictor_array);
    
    %Identify incorrect/bad cell detections
    missCells=[];
    for ii=1:length(out)
        if label(ii)==0 || delete_row(ii)==1
            missCells=[missCells; normalized_array(ii,12) normalized_array(ii,13)];
        end
    end
end

%%Clear incorrect cell detections in cellList
for jj=1:length(missCells)
    cellList.meshData{missCells(jj,1)}{missCells(jj,2)}=[];
end
clear ii jj

for ii=1:length(cellList.meshData)
    cellList.cellId{ii}=cellList.cellId{ii}(~cellfun('isempty',cellList.meshData{ii}));
    cellList.meshData{ii}=cellList.meshData{ii}(~cellfun('isempty',cellList.meshData{ii}));
    waitbar(ii/length(cellList.meshData));
end


save([pth,'\',oufti_analysis,'_SVMcurated'],'-struct', 'oufti_analysis_file');
save([pth,'\',oufti_analysis,'_SVMcurated'],'cellList', '-append'); 
close (wb);
end
                    
                    
%% Additional function to calculate euclidean distance between two points
function d=edist(x1,y1,x2,y2)
    % computes the length between 2 points
    d=sqrt((x2-x1).^2+(y2-y1).^2);
end

%% Additional function required to calculate curvature

function k=LineCurvature2D(Vertices,Lines)
%source: https://www.mathworks.com/matlabcentral/fileexchange/32696-2d-line-curvature-and-normals/content/LineCurvature2D.m 

%This function calculates the curvature of a 2D line. It first fits 
% polygons to the points. Then calculates the analytical curvature from
% the polygons;
%
%  k = LineCurvature2D(Vertices,Lines)
% 
% inputs,
%   Vertices : A M x 2 list of line points.
%   (optional)
%   Lines : A N x 2 list of line pieces, by indices of the vertices
%         (if not set assume Lines=[1 2; 2 3 ; ... ; M-1 M])
%
% outputs,
%   k : M x 1 Curvature values
%
%
%
% Example, Circle
%  r=sort(rand(15,1))*2*pi;
%  Vertices=[sin(r) cos(r)]*10;
%  Lines=[(1:size(Vertices,1))' (2:size(Vertices,1)+1)']; Lines(end,2)=1;
%  k=LineCurvature2D(Vertices,Lines);
%
%  figure,  hold on;
%  N=LineNormals2D(Vertices,Lines);
%  k=k*100;
%  plot([Vertices(:,1) Vertices(:,1)+k.*N(:,1)]',[Vertices(:,2) Vertices(:,2)+k.*N(:,2)]','g');
%  plot([Vertices(Lines(:,1),1) Vertices(Lines(:,2),1)]',[Vertices(Lines(:,1),2) Vertices(Lines(:,2),2)]','b');
%  plot(sin(0:0.01:2*pi)*10,cos(0:0.01:2*pi)*10,'r.');
%  axis equal;
%
% Example, Hand
%  load('testdata');
%  k=LineCurvature2D(Vertices,Lines);
%
%  figure,  hold on;
%  N=LineNormals2D(Vertices,Lines);
%  k=k*100;
%  plot([Vertices(:,1) Vertices(:,1)+k.*N(:,1)]',[Vertices(:,2) Vertices(:,2)+k.*N(:,2)]','g');
%  plot([Vertices(Lines(:,1),1) Vertices(Lines(:,2),1)]',[Vertices(Lines(:,1),2) Vertices(Lines(:,2),2)]','b');
%  plot(Vertices(:,1),Vertices(:,2),'r.');
%  axis equal;
%
% Function is written by D.Kroon University of Twente (August 2011)

% If no line-indices, assume a x(1) connected with x(2), x(3) with x(4) ...
if(nargin<2)
    Lines=[(1:(size(Vertices,1)-1))' (2:size(Vertices,1))'];
end

% Get left and right neighbor of each points
Na=zeros(size(Vertices,1),1); Nb=zeros(size(Vertices,1),1);
Na(Lines(:,1))=Lines(:,2); Nb(Lines(:,2))=Lines(:,1);

% Check for end of line points, without a left or right neighbor
checkNa=Na==0; checkNb=Nb==0;
Naa=Na; Nbb=Nb;
Naa(checkNa)=find(checkNa); Nbb(checkNb)=find(checkNb);

% If no left neighbor use two right neighbors, and the same for right... 
Na(checkNa)=Nbb(Nbb(checkNa)); Nb(checkNb)=Naa(Naa(checkNb));

% Correct for sampeling differences
Ta=-sqrt(sum((Vertices-Vertices(Na,:)).^2,2));
Tb=sqrt(sum((Vertices-Vertices(Nb,:)).^2,2)); 

% If no left neighbor use two right neighbors, and the same for right... 
Ta(checkNa)=-Ta(checkNa); Tb(checkNb)=-Tb(checkNb);

% Fit a polygons to the vertices 
% x=a(3)*t^2 + a(2)*t + a(1) 
% y=b(3)*t^2 + b(2)*t + b(1) 
% we know the x,y of every vertice and set t=0 for the vertices, and
% t=Ta for left vertices, and t=Tb for right vertices,  
x = [Vertices(Na,1) Vertices(:,1) Vertices(Nb,1)];
y = [Vertices(Na,2) Vertices(:,2) Vertices(Nb,2)];
M = [ones(size(Tb)) -Ta Ta.^2 ones(size(Tb)) zeros(size(Tb)) zeros(size(Tb)) ones(size(Tb)) -Tb Tb.^2];
invM=inverse3(M);
a(:,1)=invM(:,1,1).*x(:,1)+invM(:,2,1).*x(:,2)+invM(:,3,1).*x(:,3);
a(:,2)=invM(:,1,2).*x(:,1)+invM(:,2,2).*x(:,2)+invM(:,3,2).*x(:,3);
a(:,3)=invM(:,1,3).*x(:,1)+invM(:,2,3).*x(:,2)+invM(:,3,3).*x(:,3);
b(:,1)=invM(:,1,1).*y(:,1)+invM(:,2,1).*y(:,2)+invM(:,3,1).*y(:,3);
b(:,2)=invM(:,1,2).*y(:,1)+invM(:,2,2).*y(:,2)+invM(:,3,2).*y(:,3);
b(:,3)=invM(:,1,3).*y(:,1)+invM(:,2,3).*y(:,2)+invM(:,3,3).*y(:,3);

% Calculate the curvature from the fitted polygon
k = 2*(a(:,2).*b(:,3)-a(:,3).*b(:,2)) ./ ((a(:,2).^2+b(:,2).^2).^(3/2));

end

function  Minv  = inverse3(M)
% This function does inv(M) , but then for an array of 3x3 matrices
adjM(:,1,1)=  M(:,5).*M(:,9)-M(:,8).*M(:,6);
adjM(:,1,2)=  -(M(:,4).*M(:,9)-M(:,7).*M(:,6));
adjM(:,1,3)=  M(:,4).*M(:,8)-M(:,7).*M(:,5);
adjM(:,2,1)=  -(M(:,2).*M(:,9)-M(:,8).*M(:,3));
adjM(:,2,2)=  M(:,1).*M(:,9)-M(:,7).*M(:,3);
adjM(:,2,3)=  -(M(:,1).*M(:,8)-M(:,7).*M(:,2));
adjM(:,3,1)=  M(:,2).*M(:,6)-M(:,5).*M(:,3);
adjM(:,3,2)=  -(M(:,1).*M(:,6)-M(:,4).*M(:,3));
adjM(:,3,3)=  M(:,1).*M(:,5)-M(:,4).*M(:,2);
detM=M(:,1).*M(:,5).*M(:,9)-M(:,1).*M(:,8).*M(:,6)-M(:,4).*M(:,2).*M(:,9)+M(:,4).*M(:,8).*M(:,3)+M(:,7).*M(:,2).*M(:,6)-M(:,7).*M(:,5).*M(:,3);
Minv=bsxfun(@rdivide,adjM,detM);
end