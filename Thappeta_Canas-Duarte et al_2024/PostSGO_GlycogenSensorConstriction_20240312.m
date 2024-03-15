function PostSGO_GlycogenSensorConstriction_20240312()
% Author: Silvia, November 10 2023
% Importing segmentation results after SGO and constriction analysis
% Evaluate each section on a consecutive manner for debugging/ optimize parameters

clc
clear all
close all

Dir = uigetdir();% get the directory containing all the files
cd(Dir);
nameExp='GlycogenSensor_Transition'; %Name of the experiment, will be used to name the output files


allsubs = dir(fullfile(Dir, '**'));

isdir = [allsubs.isdir] & ~ismember({allsubs.name}, {'.', '..','...','phase','seg','fluor1','fluor2','fluor3','masks','cp_output','cell'});

allsubdirs = fullfile({allsubs(isdir).folder}, {allsubs(isdir).name});


% Enter sub_folder names
input_maskfolder = '\masks';
input_phasefolder = '\phase';
input_fluor1folder = '\fluor1';
input_fluor2folder = '\fluor2';
input_fluor3folder = '\fluor3';
input_cellfolder = '\cell';


%% Loop over positions

DATAAll=struct([]);

for h=1:length(allsubdirs)% reading files from one position at a time

    close all

    input_maskfolder = '\masks';
    input_phasefolder = '\phase';
    input_cellfolder = '\cell';
    input_fluor1folder = '\fluor1';
    input_fluor2folder = '\fluor2';
    input_fluor3folder = '\fluor3';

    position_folder=allsubdirs(h);
    mask_folder = strcat(allsubdirs(h), input_maskfolder);
    phase_folder = strcat(allsubdirs(h), input_phasefolder);
    cell_folder = strcat(allsubdirs(h), input_cellfolder);
    fluor1_folder = strcat(allsubdirs(h), input_fluor1folder);
    fluor2_folder = strcat(allsubdirs(h), input_fluor2folder);
    fluor3_folder = strcat(allsubdirs(h), input_fluor3folder);

    cd(string(position_folder));
    clist=dir('*clist.mat');
    data=load(clist.name);

    cd(string(mask_folder));
    labels=dir('*.png');
    L=imread(labels.name);


    DATAAll(h).data = data.data;
    DATAAll(h).statsMask = regionprops(L,'all');

    cd(string(cell_folder));
    cells=dir('*.mat');
    cells = natsortfiles({cells.name});

    for j=1:length(cells)
        load(cells{j});
        DATAAll(h).Cell(j).coord.orientation=CellA{1, 1}.coord.orientation;
        DATAAll(h).Cell(j).mask=CellA{1, 1}.mask;
        DATAAll(h).Cell(j).phase=CellA{1, 1}.phase;
        DATAAll(h).Cell(j).fluor1=CellA{1, 1}.fluor1;
        DATAAll(h).Cell(j).fluor2=CellA{1, 1}.fluor2;
        DATAAll(h).Cell(j).fluor3=CellA{1, 1}.fluor3;
    end

end
cd(Dir);
save([nameExp,'_AllCells.mat'],'DATAAll')

%% Filters potentially constricting cells based on shape factors
Pre_ConstrictionDATA=struct([]);
cellcount=0;

for pos=1:length(DATAAll)
    pos   %for monitoring the run, can be commented out
    for cell=1:length(DATAAll(pos).Cell)
        ShapeF=(DATAAll(pos).statsMask(cell).ConvexArea/DATAAll(pos).statsMask(cell).Area);

        if DATAAll(pos).statsMask(cell).MajorAxisLength >50 && DATAAll(pos).statsMask(cell).Circularity > 0.3 && ShapeF < 1.4 %Values used for this filtering were determined from plotting the histrograms of DATAAll
            cellcount=cellcount+1;
            Pre_ConstrictionDATA(cellcount).Position=pos;
            Pre_ConstrictionDATA(cellcount).CellID=cell;
            Pre_ConstrictionDATA(cellcount).Area=DATAAll(pos).data(cell,14);
            Pre_ConstrictionDATA(cellcount).Centroid=DATAAll(pos).statsMask(cell).Centroid;
            Pre_ConstrictionDATA(cellcount).Circularity=DATAAll(pos).statsMask(cell).Circularity;
            Pre_ConstrictionDATA(cellcount).ShapeF=ShapeF;
            Pre_ConstrictionDATA(cellcount).Orientation=DATAAll(pos).Cell(cell).coord.orientation;
            Pre_ConstrictionDATA(cellcount).Eccentricity=DATAAll(pos).statsMask(cell).Eccentricity;
            Pre_ConstrictionDATA(cellcount).MajorAxisLength=DATAAll(pos).statsMask(cell).MajorAxisLength;
            Pre_ConstrictionDATA(cellcount).MinorAxisLength=DATAAll(pos).statsMask(cell).MinorAxisLength;
            Pre_ConstrictionDATA(cellcount).Mask=DATAAll(pos).Cell(cell).mask;
            Pre_ConstrictionDATA(cellcount).Phase=DATAAll(pos).Cell(cell).phase;
            Pre_ConstrictionDATA(cellcount).Fluor1=DATAAll(pos).Cell(cell).fluor1;
            Pre_ConstrictionDATA(cellcount).Fluor2=DATAAll(pos).Cell(cell).fluor2;
            Pre_ConstrictionDATA(cellcount).Fluor3=DATAAll(pos).Cell(cell).fluor3;

        end
    end
end

save([nameExp,'_PreCellConstriction.mat'],'Pre_ConstrictionDATA','-v7.3')

%% Curation of Constricting cells
% Displays an overlay of each potentially constricting cells and asks user if it is indeed a constricting cell or not.

ConstrictionAnalysis=struct([]);
constricting_cell=0;

for i=1:length(Pre_ConstrictionDATA)

    clear A P mask
    close all

    A = imresize(Pre_ConstrictionDATA(i).Mask,25,'lanczos3');
    P = imresize(Pre_ConstrictionDATA(i).Phase,25,'lanczos3');
    mask = boundarymask(A);

    figure(1);
    imshow(labeloverlay(P,mask, 'ColorMap', 'spring', 'Transparency', 0.5));
    zoom on;
    pause;
    zoom off;

    % Ask the user if the cells are constricting
    answer = questdlg(sprintf('Is %d a constricting cell?', i), 'Constricting?', 'Yes', 'No', 'Yes');

    if strcmp(answer, 'Yes')
        constricting_cell=constricting_cell+1;

        ConstrictionAnalysis(constricting_cell).Position=Pre_ConstrictionDATA(i).Position;
        ConstrictionAnalysis(constricting_cell).CellID=Pre_ConstrictionDATA(i).CellID;
        ConstrictionAnalysis(constricting_cell).Area=Pre_ConstrictionDATA(i).Area;
        ConstrictionAnalysis(constricting_cell).Centroid=Pre_ConstrictionDATA(i).Centroid;
        ConstrictionAnalysis(constricting_cell).Circularity=Pre_ConstrictionDATA(i).Circularity;
        ConstrictionAnalysis(constricting_cell).ShapeF=Pre_ConstrictionDATA(i).ShapeF;
        ConstrictionAnalysis(constricting_cell).Orientation=Pre_ConstrictionDATA(i).Orientation;
        ConstrictionAnalysis(constricting_cell).Eccentricity=Pre_ConstrictionDATA(i).Eccentricity;
        ConstrictionAnalysis(constricting_cell).MajorAxisLength=Pre_ConstrictionDATA(i).MajorAxisLength;
        ConstrictionAnalysis(constricting_cell).MinorAxisLength=Pre_ConstrictionDATA(i).MinorAxisLength;
        ConstrictionAnalysis(constricting_cell).Mask=Pre_ConstrictionDATA(i).Mask;
        ConstrictionAnalysis(constricting_cell).Phase=Pre_ConstrictionDATA(i).Phase;
        ConstrictionAnalysis(constricting_cell).Fluor1=Pre_ConstrictionDATA(i).Fluor1;
        ConstrictionAnalysis(constricting_cell).Fluor2=Pre_ConstrictionDATA(i).Fluor2;
        ConstrictionAnalysis(constricting_cell).Fluor3=Pre_ConstrictionDATA(i).Fluor3;
    end

end

%%
% This section computes the midcell axis and identifies the constriction plane for each constricting cell.
% Segmentation of the glycogen sensor signal, as well as an alternative lenght estimation are also performed here. 
% Defective masks can trigger the code to error out, which can be fixed by
% removing that cell from the analysis using ConstrictionAnalysis(k) = [];

nameExp='GlycogenSensor_Transition'; %Redundant, used in case only this part of the code is being runned

ConstrictingCells_Profiles=struct([]);
ConstrictingCells_Profiles2=struct([]);
count_last=0;
count_last2=0;

for k=1:length(ConstrictionAnalysis)

    clear A A1 AF1 AF2 AP mask rows cols y x1 x2 LengthUniqueY1 UniqueY1 distance max_distance
    clear equidistant_y equidistant_x p fit_y fit_x rows2 row_idx LeftHalf RightHalf DistanceLeft DistanceRight
    clear uniqueRightHalf uniqueLeftHalf invLeft_idx2 invLeft_idx
    close all

    A1 = imresize(ConstrictionAnalysis(k).Mask,25,'lanczos3');
    A = imrotate(A1,ConstrictionAnalysis(k).Orientation+90);
    AF1=imresize(ConstrictionAnalysis(k).Fluor1,25,'lanczos3');
    AF1 = imrotate(AF1,ConstrictionAnalysis(k).Orientation+90);
    AF2=imresize(ConstrictionAnalysis(k).Fluor2,25,'lanczos3');
    AF2 = imrotate(AF2,ConstrictionAnalysis(k).Orientation+90);
    AP=imresize(ConstrictionAnalysis(k).Phase,25,'lanczos3');
    AP = imrotate(AP,ConstrictionAnalysis(k).Orientation+90);

    mask = boundarymask(A);

    % Compute the horizontal line that goes through the middle of the label
    [rows, cols] = find(mask);
    y = round(mean(rows));
    x1 = min(cols);
    x2 = max(cols);
    LengthUniqueY1=length(unique(rows));
    UniqueY1=unique(rows);

    % Compute the distance matrix from each boundary pixel to the horizontal line
    distances = abs(rows - y);

    % Compute the maximum distance
    max_distance = max(distances);

    % Compute the equidistant line
    equidistant_y = zeros(size(rows));
    for i = 1:length(rows)
        row = rows(i);
        col = cols(i);
        if mask(row, col)
            distance = distances(i);
            equidistant_y(i) = y + max_distance * sign(row - y) * distance / max_distance;
        else
            equidistant_y(i) = NaN;
        end
    end
    equidistant_x = cols;

    % Compute the bivariate fit to the center line based on the distance
    p = polyfit(equidistant_y(~isnan(equidistant_y)), equidistant_x(~isnan(equidistant_y)), 2);
    fit_y = linspace(min(equidistant_y), max(equidistant_y), LengthUniqueY1);
    fit_x = polyval(p, fit_y);

    %Finds constriction point
    [rows2,row_idx]=sort(rows);

    count_left=0;
    count_right=0;
    currentY=0;

    for n=1:length(rows2)
        for j=1:LengthUniqueY1
            if currentY< rows2(n) && rows2(n)<=fit_y(j)
                if cols(row_idx(n))<fit_x(j)
                    count_left=count_left+1;
                    LeftHalf(1,count_left)=cols(row_idx(n));
                    LeftHalf(2,count_left)=rows2(n);
                    DistanceLeft(1,count_left)=sqrt((LeftHalf(1,count_left)-fit_x(j))^2);
                    DistanceLeft(2,count_left)=LeftHalf(2,count_left);
                else
                    count_right=count_right+1;
                    RightHalf(1,count_right)=cols(row_idx(n));
                    RightHalf(2,count_right)=rows2(n);
                    DistanceRight(1,count_right)=sqrt((RightHalf(1,count_right)-fit_x(j))^2);
                    DistanceRight(2,count_right)=RightHalf(2,count_right);

                end
            end
            currentY=fit_y(j);
        end
    end

    DistanceLeft(1,:)=smoothdata(DistanceLeft(1,:),'gaussian',250);
    DistanceRight(1,:)=smoothdata(DistanceRight(1,:),'gaussian',250);

    [uniqueDistanceLeft,invLeft_idx,invLeft_idx2]=unique(DistanceLeft(2,:));
    uniqueDistanceLeft(1,:)=DistanceLeft(1,invLeft_idx);
    uniqueDistanceLeft(2,:)=DistanceLeft(2,invLeft_idx);
    uniqueLeftHalf(1,:)=LeftHalf(1,invLeft_idx);
    uniqueLeftHalf(2,:)=LeftHalf(2,invLeft_idx);

    [uniqueDistanceRight,invRight_idx,invRight_idx2]=unique(DistanceRight(2,:));
    uniqueDistanceRight(1,:)=DistanceRight(1,invRight_idx);
    uniqueDistanceRight(2,:)=DistanceRight(2,invRight_idx);
    uniqueRightHalf(1,:)=RightHalf(1,invRight_idx);
    uniqueRightHalf(2,:)=RightHalf(2,invRight_idx);

    invDistanceLeft= max(uniqueDistanceLeft(1,:)) - uniqueDistanceLeft(1,:);
    [pks_left,locs_left,w_left,p_left]=findpeaks(invDistanceLeft(25:end-25),uniqueDistanceLeft(2,(25:end-25)),'MinPeakProminence',20,'MinPeakDistance',15,'Annotate','extents','WidthReference','halfheight');

    for n=1:length(invDistanceLeft)
        if invDistanceLeft(n)==pks_left(1)
            coord_pk_left(1,1)=uniqueLeftHalf(1,n);
            coord_pk_left(2,1)=uniqueLeftHalf(2,n);
        end
    end

    invDistanceRight= max(uniqueDistanceRight(1,:)) - uniqueDistanceRight(1,:);
    [pks_right,locs_right,w_right,p_right]=findpeaks(invDistanceRight(25:end-25),uniqueDistanceRight(2,(25:end-25)),'MinPeakProminence',20,'MinPeakDistance',15,'Annotate','extents','WidthReference','halfheight');

    for n=1:length(invDistanceRight)
        if invDistanceRight(n)==pks_right(1)
            coord_pk_right(1,1)=uniqueRightHalf(1,n);
            coord_pk_right(2,1)=uniqueRightHalf(2,n);
        end
    end

    x_divplane=[coord_pk_left(2,1) coord_pk_right(2,1)];
    y_divplane=[coord_pk_left(1,1) coord_pk_right(1,1)];

    [x_int, y_int] = polyxpoly(x_divplane, y_divplane, fit_y, fit_x); %x_int

    x_divpoint=0;

    for n=1:length(fit_y)
        if fit_y(n)==round(x_int)
            x_divpoint=n;
        end
    end

%     figure(1)
%     hold on
%     %scatter(fit_x,fit_y)
%     %scatter(cols,rows)
%     scatter(RightHalf(1,:),RightHalf(2,:))
%     scatter(LeftHalf(1,:),LeftHalf(2,:))
%     % scatter(uniqueDistanceLeft(2,:),uniqueDistanceLeft(1,:))
%     % scatter(uniqueDistanceRight(2,:),uniqueDistanceRight(1,:))
%     hold off

    %Display the label and the midcell line
    %Also saves figure2 for checking accurate cell mapping

    figure(2)
    imshow(AP);
    hold on;
    plot(fit_x, fit_y, 'b', 'LineWidth', 2);
    plot(y_divplane,x_divplane, 'r', 'LineWidth', 2);
    scatter(y_int,x_int,'y','LineWidth',3)
    hold off;
    fig2_name = sprintf('_Figure2_%d.png', k);
    saveas(gcf,[nameExp,fig2_name]);

    %Segmentation of Glycogen sensor signal

    AF2_1=AF2*2;
    v0 = max(AF2_1(:));
    if v0 > 0.975*65535
        AF2_1=AF2_1*0.75;
    elseif v0 < 0.475*65535
        AF2_1=AF2_1*2;
    elseif v0 < 0.275*65535
        AF2_1=AF2_1*3;
    end

    %Apply background subtraction
    background = imopen(AF2_1, strel('disk', 300));
    AF2_a = AF2_1 - background;

    % Signal amplification
    AF2_b=AF2_a+AF2_a;
    v1 = max(AF2_b(:));
    if v1 > 0.975*65535
        AF2_b=AF2_b*0.75;
    elseif v1 < 0.475*65535
        AF2_b=AF2_b*2;
    elseif v1 < 0.275*65535
        AF2_b=AF2_b*3;
    end

    % Contrast adjustment

    AF2_c = locallapfilt(AF2_b, 0.4, 0.8);

    v2 = max(AF2_c(:));
    if v2 > 0.975*65535
        AF2_c=AF2_c*0.75;
    elseif v2 < 0.475*65535
        AF2_c=AF2_c*2;
    elseif v2 < 0.275*65535
        AF2_c=AF2_c*3;
    end

    AF2_d = medfilt2(AF2_c, [25 25]);
    AF2_e = (0.85*AF2_c) + (AF2_d*0.15);

    v3 = max(AF2_e(:));
    if v3 > 0.975*65535
        AF2_e=AF2_e*0.75;
    elseif v3 < 0.475*65535
        AF2_e=AF2_e*2;
    elseif v3 < 0.275*65535
        AF2_e=AF2_e*3;
    end
    background2 = imopen(AF2_e, strel('disk', 300));
    AF2_f = AF2_e - background2;

    % Apply Otsu's method to find the global threshold
    global_thresh = graythresh(AF2_f);

    % Apply adaptive thresholding to segment the image

    binary_AF2_1 = imbinarize(AF2_f,1.2*global_thresh);
    binary_AF2_2 = imbinarize(AF2_f,"adaptive", 'Sensitivity', 0.5); 

    label_matrix = logical(binary_AF2_1);
    regions=regionprops(label_matrix,'all');

    masked_AF2_1 = immultiply(A, binary_AF2_1);
    masked_AF2_2 = immultiply(A, binary_AF2_2);

    %Remove small objects from the image
    clean_AF2_1 = bwareaopen(masked_AF2_1, 15000);
    clean_AF2_2 = bwareaopen(masked_AF2_2, 15000);

    se=strel('disk',25);
    closed_AF2 = imclose(clean_AF2_2,se);
    filled_AF2 = imfill(closed_AF2, 'holes');
    opened_AF2 = imopen(filled_AF2,se);

    se=strel('disk',15);
    eroded_AF2_2 = imerode(opened_AF2, se);
    se=strel('disk',25);
    eroded_AF2_1 = imerode(clean_AF2_1, se);

    final_AF2=eroded_AF2_1;

    if length(regions) <2
        final_AF2=eroded_AF2_2;
    elseif length(regions) >=2
        regionsArea=0;
        for q=1:length(regions)
            regionsArea=regionsArea+regions(q).Area;
        end
        if regionsArea/(25*25) > 0.5*ConstrictionAnalysis(k).Area
            final_AF2=eroded_AF2_2;
        end
    else
        final_AF2=eroded_AF2_1;
    end

    final_AF2 = immultiply(A, final_AF2);

    %Display the original image and the segmented image side by side

    fig3=figure(3);
    hold on;
    subplot(1, 3, 1);
    hold on;
    imshow(AP);
    plot(fit_x, fit_y, 'b', 'LineWidth', 1);
    plot(y_divplane,x_divplane, 'r', 'LineWidth', 1);
    scatter(y_int,x_int,'y','LineWidth',1)
    title('Original Image');
    hold off;
    subplot(1, 3, 2);
    imshow(AF2_f);
    title('LaPlacian filter');
    subplot(1, 3, 3);
    hold on;
    imshow(final_AF2);
    plot(fit_x, fit_y, 'b', 'LineWidth', 1);
    plot(y_divplane,x_divplane, 'r', 'LineWidth', 1);
    scatter(y_int,x_int,'y','LineWidth',1)
    title('Segmented Image');
    hold off;
    hold off;
    fig3_name = sprintf('_Figure3_%d.png', k);
    saveas(gcf,[nameExp,fig3_name]);
    close(fig3);

    valuesF1 = interp2(double(AF1), fit_x, fit_y);
    valuesF2 = interp2(double(AF2), fit_x, fit_y);
    valuesP = interp2(double(AP), fit_x, fit_y);

    d1_valuesF1 = interp2(double(AF1), fit_x(1:x_divpoint), fit_y(1:x_divpoint));
    d1_valuesF2 = interp2(double(AF2), fit_x(1:x_divpoint), fit_y(1:x_divpoint));
    d1_valuesP = interp2(double(AP), fit_x(1:x_divpoint), fit_y(1:x_divpoint));

    d2_valuesF1 = interp2(double(AF1), fit_x(x_divpoint:end), fit_y(x_divpoint:end));
    d2_valuesF2 = interp2(double(AF2), fit_x(x_divpoint:end), fit_y(x_divpoint:end));
    d2_valuesP = interp2(double(AP), fit_x(x_divpoint:end), fit_y(x_divpoint:end));

    norm_valuesF1=valuesF1./max(valuesF1);
    norm_valuesF2=valuesF2./max(valuesF2);
    norm_valuesP=valuesF1./max(valuesP);

    norm_d1_valuesF1=d1_valuesF1./max(valuesF1);
    norm_d1_valuesF2=d1_valuesF2./max(valuesF2);
    norm_d1_valuesP=d1_valuesF1./max(valuesP);

    norm_d2_valuesF1=d2_valuesF1./max(valuesF1);
    norm_d2_valuesF2=d2_valuesF2./max(valuesF2);
    norm_d2_valuesP=d2_valuesF1./max(valuesP);

    %Estimate glycogen caps lenght as an alternative to using the segmented areas

    rev_norm_d2_valuesF1=fliplr(norm_d2_valuesF1);
    rev_norm_d2_valuesF2=fliplr(norm_d2_valuesF2);
    rev_norm_d2_valuesP=fliplr(norm_d2_valuesP);

    d1_valleyF2=mean(norm_d1_valuesF2(end-round(0.5*length(norm_d1_valuesF2)):end-round(0.25*length(norm_d1_valuesF2))));
    d2_valleyF2=mean(rev_norm_d2_valuesF2(end-round(0.5*length(rev_norm_d2_valuesF2)):end-round(0.25*length(rev_norm_d2_valuesF2))));

    gFactor=1.1;

    valleyF2=gFactor*((d1_valleyF2+d2_valleyF2)/2);

    [pks_cF2,locs_cF2,w_cF2,p_cF2]=findpeaks(norm_valuesF2,'MinPeakProminence',0.15,'MinPeakDistance',500,'Annotate','extents','WidthReference','halfheight');

    if ~isempty(w_cF2)
        c_widthF2=w_cF2(1);
    else
        c_widthF2=0;
    end

    [pks_d1F2,locs_d1F2,w_d1F2,p_d1F2]=findpeaks(norm_d1_valuesF2,'MinPeakProminence',0.15,'MinPeakDistance',500,'Annotate','extents','WidthReference','halfheight');

    if ~isempty(locs_d1F2)
        d1_locF2=locs_d1F2(1);
    else
        d1_locF2=1;
    end

    [pks_d2F2,locs_d2F2,w_d2F2,p_d2F2]=findpeaks(rev_norm_d2_valuesF2,'MinPeakProminence',0.2,'MinPeakDistance',250,'Annotate','extents','WidthReference','halfheight');
    if ~isempty(locs_d2F2)
        d2_locF2=locs_d2F2(1);
    else
        d2_locF2=1;
    end

    d1_widthF2=0;

    if d1_locF2 > 1
        n=d1_locF2;
        d1_widthF2=d1_locF2;
        limit_d1=length(norm_d1_valuesF2)-2;
        while (norm_d1_valuesF2(n)+norm_d1_valuesF2(n+1)+norm_d1_valuesF2(n+2))/3 > gFactor*d1_valleyF2
            if (norm_d1_valuesF2(n)+norm_d1_valuesF2(n+1)+norm_d1_valuesF2(n+2))/3 <= gFactor*d1_valleyF2
                break
            end
            if n>= limit_d1
                break
            end
            d1_widthF2=d1_widthF2+1;
            n=n+1;
        end
    elseif d1_locF2 <=1
        for n=1:round(0.5*length(norm_d1_valuesF2))-1
            if (norm_d1_valuesF2(n)+norm_d1_valuesF2(n+1))/2 >= gFactor*d1_valleyF2
                d1_widthF2=d1_widthF2+1;
            end
        end
    end
    d1_widthF2=(d1_widthF2/1);

    d2_widthF2=0;

    if d2_locF2 > 1
        n=d2_locF2;
        d2_widthF2=d2_locF2;
        limit_d2=length(rev_norm_d2_valuesF2)-2;
        while (rev_norm_d2_valuesF2(n)+rev_norm_d2_valuesF2(n+1)+rev_norm_d2_valuesF2(n+2))/3 > gFactor*d2_valleyF2
            if (rev_norm_d2_valuesF2(n)+rev_norm_d2_valuesF2(n+1)+rev_norm_d2_valuesF2(n+2))/3 <= gFactor*d2_valleyF2
                break
            end
            if n>=limit_d2
                break
            end
            d2_widthF2=d2_widthF2+1;
            n=n+1;
        end
    elseif d2_locF2 <=1
        for n=1:round(0.5*length(rev_norm_d2_valuesF2))-1
            if (rev_norm_d2_valuesF2(n)+rev_norm_d2_valuesF2(n+1))/2 >= gFactor*d2_valleyF2
                d2_widthF2=d2_widthF2+1;
            end
        end
    end
    d2_widthF2=(d2_widthF2/1);

    diff_widthF2=sqrt((d2_widthF2-d1_widthF2)^2);

    d1_length=(length(fit_y(1:x_divpoint))/1);
    d2_length=(length(fit_y(x_divpoint:end))/1);

    diff_lenght=sqrt((d1_length-d2_length)^2);

    if d1_length>d2_length
        abs_diff_widthF2=d1_widthF2-d2_widthF2;
        abs_diff_lenght=d1_length-d2_length;
    else
        abs_diff_widthF2=d2_widthF2-d1_widthF2;
        abs_diff_lenght=d2_length-d1_length;
    end

    norm_diff_widthF2=sqrt((d2_widthF2-d1_widthF2)^2)/(d1_length+d2_length);
    norm_diff_lenght=sqrt((d1_length-d2_length)^2)/(d1_length+d2_length);

    %Measure areas of Glycogen caps and daugther cells

    d1_area=0;
    d2_area=0;
    d1_GS_area=0;
    d2_GS_area=0;

    se=strel('disk',1);
    eroded_A = imerode(A, se);


    for m=1:size(A,2)
        for n=1:round(x_int)
            if eroded_A(n,m)==1
                d1_area=d1_area+1;
            end
        end
        for n=round(x_int):size(A,1)
            if eroded_A(n,m)==1
                d2_area=d2_area+1;
            end
        end
    end

    % Convert the binary image to a label matrix
    glylabel_matrix = bwlabel(final_AF2);

    GlyRegions=regionprops(glylabel_matrix,'all');

    GS_Regions=zeros(length(GlyRegions),2);
    for i=1:length(GlyRegions)
        GS_Regions(i,1)=GlyRegions(i).Area;
        GS_Regions(i,2)=GlyRegions(i).Centroid(2);
    end

    [~, idx] = sort(GS_Regions(:,end),1);
    sorted_GS_Regions = GS_Regions(idx,:);

    if length(sorted_GS_Regions) >1
        d1_GS_area=sorted_GS_Regions(1,1);
        d2_GS_area=sorted_GS_Regions(end,1);
    end

    d1_area=d1_area*(0.07*0.07)/(25*25);
    d2_area=d2_area*(0.07*0.07)/(25*25);
    d1_GS_area=d1_GS_area*(0.07*0.07)/(25*25);
    d2_GS_area=d2_GS_area*(0.07*0.07)/(25*25);

    if d1_GS_area == d2_GS_area
        diff_area=abs(d1_area-d2_area);
        diff_GS_area=d1_GS_area-10*(0.07*0.07);
    elseif d1_area>d2_area && d1_GS_area >0 && d2_GS_area >0
        diff_area=d1_area-d2_area;
        diff_GS_area=d1_GS_area-d2_GS_area;
    elseif d2_area>d1_area && d1_GS_area >0 && d2_GS_area >0
        diff_area=d2_area-d1_area;
        diff_GS_area=d2_GS_area-d1_GS_area;
    else
        diff_area=NaN;
        diff_GS_area=NaN;
    end

    ConstrictionAnalysis(k).LenghtDiff=diff_lenght;
    ConstrictionAnalysis(k).F2WidthDiff=diff_widthF2;
    ConstrictionAnalysis(k).d1_widthF2=d1_widthF2;
    ConstrictionAnalysis(k).d1_lengthF2=d1_length;
    ConstrictionAnalysis(k).d2_widthF2=d2_widthF2;
    ConstrictionAnalysis(k).d2_lengthF2=d2_length;
    ConstrictionAnalysis(k).Abs_LenghtDiff=abs_diff_lenght;
    ConstrictionAnalysis(k).Abs_F2WidthDiff=abs_diff_widthF2;
    ConstrictionAnalysis(k).normLenghtDiff=norm_diff_lenght;
    ConstrictionAnalysis(k).normF2WidthDiff=norm_diff_widthF2;

    ConstrictionAnalysis(k).D1valuesF1=d1_valuesF1;
    ConstrictionAnalysis(k).D2valuesF1=d2_valuesF1;
    ConstrictionAnalysis(k).D1valuesF1=d1_valuesF2;
    ConstrictionAnalysis(k).D2valuesF1=d2_valuesF2;
    ConstrictionAnalysis(k).D1valuesP=d1_valuesP;
    ConstrictionAnalysis(k).D2valuesP=d2_valuesP;

    ConstrictionAnalysis(k).normD1valuesF1=norm_d1_valuesF1;
    ConstrictionAnalysis(k).normD2valuesF1=norm_d2_valuesF1;
    ConstrictionAnalysis(k).normD1valuesF1=norm_d1_valuesF2;
    ConstrictionAnalysis(k).normD2valuesF1=norm_d2_valuesF2;
    ConstrictionAnalysis(k).normD1valuesP=norm_d1_valuesP;
    ConstrictionAnalysis(k).normD2valuesP=norm_d2_valuesP;

    ConstrictionAnalysis(k).normF1=norm_valuesF1;
    ConstrictionAnalysis(k).normF2=norm_valuesF2;
    ConstrictionAnalysis(k).meanF2=mean(valuesF2);
    ConstrictionAnalysis(k).medianF2=median(valuesF2);
    ConstrictionAnalysis(k).maxF2=max(valuesF2);
    ConstrictionAnalysis(k).normP=norm_valuesP;

    ConstrictionAnalysis(k).d1Valley=d1_valleyF2;
    ConstrictionAnalysis(k).d2Valley=d2_valleyF2;

    ConstrictionAnalysis(k).cPksF2=pks_cF2;
    ConstrictionAnalysis(k).Num_cPksF2=length(pks_cF2);
    ConstrictionAnalysis(k).cWidthsF2=w_cF2;

    ConstrictionAnalysis(k).CalcCellLength=length(fit_y(1:end));

    ConstrictionAnalysis(k).AP=AP;
    ConstrictionAnalysis(k).AP_maks=mask;
    ConstrictionAnalysis(k).AF2_b=AF2_b;
    ConstrictionAnalysis(k).AF2_f=AF2_f;
    ConstrictionAnalysis(k).AF2_maks=final_AF2;
    ConstrictionAnalysis(k).AreaMicrons=ConstrictionAnalysis(k).Area*(0.07*0.07);

    ConstrictionAnalysis(k).d1_area=d1_area;
    ConstrictionAnalysis(k).d2_area=d2_area;
    ConstrictionAnalysis(k).diff_area=diff_area;
    ConstrictionAnalysis(k).d1_GS_area=d1_GS_area;
    ConstrictionAnalysis(k).d2_GS_area=d2_GS_area;
    ConstrictionAnalysis(k).diff_GS_area=diff_GS_area;

    for h=1:length(fit_y(1:end))
        ConstrictionAnalysis(k).NormCellLength(h)=h/length(fit_y(1:end));
    end

    if ConstrictionAnalysis(k).Num_cPksF2 >0 && ConstrictionAnalysis(k).Num_cPksF2<3 && ConstrictionAnalysis(k).medianF2 < 10000 && ConstrictionAnalysis(k).medianF2 > 3250

        count_last=count_last+1;
        ConstrictingCells_Profiles(count_last).Position=ConstrictionAnalysis(k).Position;
        ConstrictingCells_Profiles(count_last).CellID=ConstrictionAnalysis(k).CellID;
        ConstrictingCells_Profiles(count_last).normLenghtDiff=ConstrictionAnalysis(k).normLenghtDiff;
        ConstrictingCells_Profiles(count_last).normF2WidthDiff=ConstrictionAnalysis(k).normF2WidthDiff;
        ConstrictingCells_Profiles(count_last).LenghtDiff=ConstrictionAnalysis(k).LenghtDiff;
        ConstrictingCells_Profiles(count_last).F2WidthDiff=ConstrictionAnalysis(k).F2WidthDiff;
        ConstrictingCells_Profiles(count_last).Abs_LenghtDiff=ConstrictionAnalysis(k).Abs_LenghtDiff;
        ConstrictingCells_Profiles(count_last).Abs_F2WidthDiff=ConstrictionAnalysis(k).Abs_F2WidthDiff;
        ConstrictingCells_Profiles(count_last).d1_area=ConstrictionAnalysis(k).d1_area;
        ConstrictingCells_Profiles(count_last).d2_area=ConstrictionAnalysis(k).d2_area;
        ConstrictingCells_Profiles(count_last).diff_area=ConstrictionAnalysis(k).diff_area;
        ConstrictingCells_Profiles(count_last).d1_GS_area=ConstrictionAnalysis(k).d1_GS_area;
        ConstrictingCells_Profiles(count_last).d2_GS_area=ConstrictionAnalysis(k).d2_GS_area;
        ConstrictingCells_Profiles(count_last).diff_GS_area=ConstrictionAnalysis(k).diff_GS_area;

    end

    if ConstrictionAnalysis(k).Num_cPksF2 >1 && ConstrictionAnalysis(k).Num_cPksF2<3 && ConstrictionAnalysis(k).medianF2 < 10000 && ConstrictionAnalysis(k).medianF2 > 3250

        count_last2=count_last2+1;
        ConstrictingCells_Profiles2(count_last2).Position=ConstrictionAnalysis(k).Position;
        ConstrictingCells_Profiles2(count_last2).CellID=ConstrictionAnalysis(k).CellID;
        ConstrictingCells_Profiles2(count_last2).normLenghtDiff=ConstrictionAnalysis(k).normLenghtDiff;
        ConstrictingCells_Profiles2(count_last2).normF2WidthDiff=ConstrictionAnalysis(k).normF2WidthDiff;
        ConstrictingCells_Profiles2(count_last2).LenghtDiff=ConstrictionAnalysis(k).LenghtDiff;
        ConstrictingCells_Profiles2(count_last2).F2WidthDiff=ConstrictionAnalysis(k).F2WidthDiff;
        ConstrictingCells_Profiles2(count_last2).Abs_LenghtDiff=ConstrictionAnalysis(k).Abs_LenghtDiff;
        ConstrictingCells_Profiles2(count_last2).Abs_F2WidthDiff=ConstrictionAnalysis(k).Abs_F2WidthDiff;
        ConstrictingCells_Profiles2(count_last2).d1_area=ConstrictionAnalysis(k).d1_area;
        ConstrictingCells_Profiles2(count_last2).d2_area=ConstrictionAnalysis(k).d2_area;
        ConstrictingCells_Profiles2(count_last2).diff_area=ConstrictionAnalysis(k).diff_area;
        ConstrictingCells_Profiles2(count_last2).d1_GS_area=ConstrictionAnalysis(k).d1_GS_area;
        ConstrictingCells_Profiles2(count_last2).d2_GS_area=ConstrictionAnalysis(k).d2_GS_area;
        ConstrictingCells_Profiles2(count_last2).diff_GS_area=ConstrictionAnalysis(k).diff_GS_area;

    end

end

for i=1:length(ConstrictionAnalysis)
    ConstrictionAnalysis(i).RatioAreas=ConstrictionAnalysis(i).diff_area/ConstrictionAnalysis(i).diff_GS_area;
    ConstrictionAnalysis(i).DiffAreasPercentage=100*ConstrictionAnalysis(i).diff_area/ConstrictionAnalysis(i).AreaMicrons;
end

%%
AsymmetricCells=struct([]);
count_asymmetric=0;

for i=1:length(ConstrictionAnalysis)
    if ConstrictionAnalysis(i).DiffAreasPercentage>3.75 % Filters cells in which the difference in area of future daugthers size is under a threshold (i.e, 3.75)
        count_asymmetric=count_asymmetric+1;
        AsymmetricCells(count_asymmetric).Position=ConstrictionAnalysis(i).Position;
        AsymmetricCells(count_asymmetric).CellID=ConstrictionAnalysis(i).CellID;
        AsymmetricCells(count_asymmetric).Area=ConstrictionAnalysis(i).Area;
        AsymmetricCells(count_asymmetric).AreaMicrons=ConstrictionAnalysis(i).AreaMicrons;
        AsymmetricCells(count_asymmetric).d1_area=ConstrictionAnalysis(i).d1_area;
        AsymmetricCells(count_asymmetric).d2_area=ConstrictionAnalysis(i).d2_area;
        AsymmetricCells(count_asymmetric).diff_area=ConstrictionAnalysis(i).diff_area;
        AsymmetricCells(count_asymmetric).d1_GS_area=ConstrictionAnalysis(i).d1_GS_area;
        AsymmetricCells(count_asymmetric).d2_GS_area=ConstrictionAnalysis(i).d2_GS_area;
        AsymmetricCells(count_asymmetric).diff_GS_area=ConstrictionAnalysis(i).diff_GS_area;
        AsymmetricCells(count_asymmetric).DiffAreaPercentage=ConstrictionAnalysis(i).DiffAreasPercentage;
        AsymmetricCells(count_asymmetric).Number1=i;
        AsymmetricCells(count_asymmetric).Number2=count_asymmetric;

        fig4=figure(4);
        hold on;
        subplot(1, 3, 1);
        hold on;
        imshow(ConstrictionAnalysis(i).AP);
        title('Phase Image');
        hold off;
        subplot(1, 3, 2);
        imshow(ConstrictionAnalysis(i).AF2_maks);
        title('GS mask');
        subplot(1, 3, 3);
        hold on;
        imshow(labeloverlay(ConstrictionAnalysis(i).AF2_b,ConstrictionAnalysis(i).AP_maks, 'ColorMap', 'spring', 'Transparency', 0));
        title('Segmented Image');
        hold off;
        hold off;
        fig4_name = sprintf('_Figure4_%d.png', count_asymmetric);
        saveas(gcf,[nameExp,fig4_name]);
        close(fig4);

    end
end

%% Plotting 

x_1=[AsymmetricCells.diff_GS_area]';

x1=linspace(0,1,1000);
y1=x1;

yfit1=0.99754119*x_1 + 0.12870224; %Values from principal component regression analysis (PCR), other linear fits can be used instead

figure('Name','5','WindowState','maximized');
hold on
q5=tiledlayout(1,1); % Requires R2019b or later

nexttile
hold on;

scatter([AsymmetricCells.diff_GS_area],[AsymmetricCells.diff_area]);
plot(x_1,yfit1,'k-.');
set(gca,'FontSize',12)
set(gca,'FontName','Arial')
ylabel('\DeltaArea','FontSize',12)
xlabel('\Delta Glycogen sensor Area','FontSize',12)
xlim([-0.05 1]);
ylim([0 1]);
hold off;


[CCf CPf]=corrcoef([AsymmetricCells.diff_GS_area],[AsymmetricCells.diff_area]);

[rho_f, pval_f] = corr([AsymmetricCells.diff_GS_area]',[AsymmetricCells.diff_area]', 'Type', 'Spearman');






