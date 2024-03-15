function PostSGO_GlycogenSensorAllCells_20240312()
% Author: Silvia, February 15 2024
% Importing segmentation results after SGO and constriction analysis

%%% PARAMETERS %%%
clc
clear all
close all

Dir = uigetdir();% get the directory containing all the files
cd(Dir);
nameExp='GlycogenSensor_Transition';


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

%%
Pre_AllCellsDATA=struct([]);
cellcount=0;

for pos=1:length(DATAAll)
    for cell=1:length(DATAAll(pos).Cell)
        ShapeF=(DATAAll(pos).statsMask(cell).ConvexArea/DATAAll(pos).statsMask(cell).Area);

        if DATAAll(pos).statsMask(cell).MajorAxisLength >25 && DATAAll(pos).statsMask(cell).MajorAxisLength <65 && DATAAll(pos).statsMask(cell).MinorAxisLength >13.5 && DATAAll(pos).statsMask(cell).MinorAxisLength <17.5 && DATAAll(pos).statsMask(cell).Circularity > 0.3 && ShapeF < 1.4
            cellcount=cellcount+1;
            Pre_AllCellsDATA(cellcount).Position=pos;
            Pre_AllCellsDATA(cellcount).CellID=cell;
            Pre_AllCellsDATA(cellcount).Area=DATAAll(pos).data(cell,14);
            Pre_AllCellsDATA(cellcount).Centroid=DATAAll(pos).statsMask(cell).Centroid;
            Pre_AllCellsDATA(cellcount).Circularity=DATAAll(pos).statsMask(cell).Circularity;
            Pre_AllCellsDATA(cellcount).ShapeF=ShapeF;
            Pre_AllCellsDATA(cellcount).Orientation=DATAAll(pos).Cell(cell).coord.orientation;
            Pre_AllCellsDATA(cellcount).Eccentricity=DATAAll(pos).statsMask(cell).Eccentricity;
            Pre_AllCellsDATA(cellcount).MajorAxisLength=DATAAll(pos).statsMask(cell).MajorAxisLength;
            Pre_AllCellsDATA(cellcount).MinorAxisLength=DATAAll(pos).statsMask(cell).MinorAxisLength;
            Pre_AllCellsDATA(cellcount).Mask=DATAAll(pos).Cell(cell).mask;
            Pre_AllCellsDATA(cellcount).Phase=DATAAll(pos).Cell(cell).phase;
            Pre_AllCellsDATA(cellcount).Fluor1=DATAAll(pos).Cell(cell).fluor1;
            Pre_AllCellsDATA(cellcount).Fluor2=DATAAll(pos).Cell(cell).fluor2;
            Pre_AllCellsDATA(cellcount).Fluor3=DATAAll(pos).Cell(cell).fluor3;

        end
    end
end

save([nameExp,'_PreAllCells.mat'],'Pre_AllCellsDATA','-v7.3')

%%
AllCellsDATA=Pre_AllCellsDATA; %Renames the structure for next analyses

%%
% This section computes the midcell axis and identifies the constriction centroid for each constricting cell.
% Segmentation of the glycogen sensor signal, as well the nucleoid (from HU-mCherry signal) and the nucleoid centroid position.
% Defective masks can trigger the code to error out, which can be fixed by removing that cell from the analysis using AllCellsDATA(k) = [];

for k=1:length(AllCellsDATA)

    clear A A1 AF1 AF2 AF3 AP mask rows cols y x1 x2 LengthUniqueY1 UniqueY1 distance max_distance
    clear equidistant_y equidistant_x p fit_y fit_x rows2 row_idx LeftHalf RightHalf DistanceLeft DistanceRight
    clear uniqueRightHalf uniqueLeftHalf invLeft_idx2 invLeft_idx
    close all

    A1 = imresize(AllCellsDATA(k).Mask,25,'lanczos3');
    A = imrotate(A1,AllCellsDATA(k).Orientation+90);
    AF1=imresize(AllCellsDATA(k).Fluor1,25,'lanczos3');
    AF1 = imrotate(AF1,AllCellsDATA(k).Orientation+90);
    AF2=imresize(AllCellsDATA(k).Fluor2,25,'lanczos3');
    AF2 = imrotate(AF2,AllCellsDATA(k).Orientation+90);
    AF3=imresize(AllCellsDATA(k).Fluor3,25,'lanczos3');
    AF3 = imrotate(AF3,AllCellsDATA(k).Orientation+90);
    AP=imresize(AllCellsDATA(k).Phase,25,'lanczos3');
    AP = imrotate(AP,AllCellsDATA(k).Orientation+90);

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

    MidPoint_MidCell=round(length(fit_y)/2);

    for n=1:length(fit_y)
        if fit_y(n)==MidPoint_MidCell
            x_int=n;
        end
    end

    % Calculation of cell's centroid
    MidCellX=fit_y(MidPoint_MidCell);
    MidCellY=fit_x(MidPoint_MidCell);

    fig1=figure(1);
    hold on
    imshow(AP);
    hold on;
    plot(fit_x, fit_y, 'b', 'LineWidth', 2);
    scatter(MidCellY,MidCellX,'y','LineWidth',3)
    hold off;
    fig1_name = sprintf('_Figure1_%d.png', k);
    saveas(gcf,[nameExp,fig1_name]);
    close(fig1);

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
    background = imopen(AF2_1, strel('disk', 100));
    AF2_a = imsubtract(AF2_1,background);

    % Signal amplification
    AF2_b=imadd(AF2_a,AF2_a);
    AF2_b=imadd(AF2_b,AF2_a);

    v1 = max(AF2_b(:));
    if v1 > 0.975*65535
        AF2_b=immultiply(AF2_b,0.75);
    elseif v1 < 0.475*65535
        AF2_b=immultiply(AF2_b,2);
    elseif v1 < 0.275*65535
        AF2_b=immultiply(AF2_b,3);
    end

    AF2_c = locallapfilt(AF2_b, 0.4, 0.8);

    v2 = max(AF2_c(:));
    if v2 > 0.975*65535
        AF2_c=immultiply(AF2_c,0.75);
    elseif v2 < 0.475*65535
        AF2_c=immultiply(AF2_c,2);
    elseif v2 < 0.275*65535
        AF2_c=immultiply(AF2_c,3);
    end

    AF2_d = medfilt2(AF2_c, [25 25]);
    AF2_e = imadd(immultiply(AF2_c,0.85),immultiply(AF2_d,0.15));

    v3 = max(AF2_e(:));
    if v3 > 0.975*65535
        AF2_e=immultiply(AF2_e,0.75);
    elseif v3 < 0.475*65535
        AF2_e=immultiply(AF2_e,2);
    elseif v3 < 0.275*65535
        AF2_e=immultiply(AF2_e,3);
    end
    background2 = imopen(AF2_e, strel('disk', 150));
    AF2_f = imsubtract(AF2_e,background2);
    AF2_g = AF2_f;
    AF2_h = locallapfilt(AF2_g, 0.4, 0.85);

    % Apply Otsu's method to find the global threshold
    global_thresh = graythresh(AF2_h);

    % Apply adaptive thresholding to segment the image

    binary_AF2_1 = imbinarize(AF2_h,1.25*global_thresh);
    binary_AF2_2 = imbinarize(AF2_h,"adaptive", 'Sensitivity', 0.5); %Previously 0.5

    label_matrix = logical(binary_AF2_1);
    regions=regionprops(label_matrix,'all');

    masked_AF2_1 = immultiply(A, binary_AF2_1);
    masked_AF2_2 = immultiply(A, binary_AF2_2);

    %Remove small objects from the image
    clean_AF2_1 = bwareaopen(masked_AF2_1, 150); %Before 15k
    clean_AF2_2 = bwareaopen(masked_AF2_2, 150); %Before 15k

    se=strel('disk',25);
    closed_AF2 = imclose(clean_AF2_2,se);
    filled_AF2 = imfill(closed_AF2, 'holes');
    opened_AF2 = imopen(filled_AF2,se);
    closed_AF2_1 = imclose(clean_AF2_1,se);
    filled_AF2_1 = imfill(closed_AF2_1, 'holes');
    opened_AF2_1 = imopen(filled_AF2_1,se);


    final_AF2=opened_AF2_1;

    final_AF2 = immultiply(A, final_AF2);

    %Segmentation of nucleoid

    AF3_1=immultiply(AF3,2);
    v0 = max(AF3_1(:));
    if v0 > 0.975*65535
        AF3_1=immultiply(AF3_1,0.75);
    elseif v0 < 0.475*65535
        AF3_1=immultiply(AF3_1,2);
    elseif v0 < 0.275*65535
        AF3_1=immultiply(AF3_1,3);
    end

    %Apply background subtraction
    background3 = imopen(AF3_1, strel('disk', 300));
    AF3_a = imsubtract(AF3_1,background3);

    AF3_c = locallapfilt(AF3_a, 0.4, 0.5);

    v2 = max(AF3_c(:));
    if v2 > 0.975*65535
        AF3_c=immultiply(AF3_c,0.75);
    elseif v2 < 0.475*65535
        AF3_c=immultiply(AF3_c,2);
    elseif v2 < 0.275*65535
        AF3_c=immultiply(AF3_c,3);
    end

    AF3_d = medfilt2(AF3_c, [25 25]);
    AF3_e = imadd(immultiply(AF3_c,0.85),immultiply(AF3_d,0.15));

    v3 = max(AF3_e(:));
    if v3 > 0.975*65535
        AF3_e=immultiply(AF3_e,0.75);
    elseif v3 < 0.475*65535
        AF3_e=immultiply(AF3_e,2);
    elseif v3 < 0.275*65535
        AF3_e=immultiply(AF3_e,3);
    end
    background4 = imopen(AF3_e, strel('disk', 300));
    AF3_f = imsubtract(AF3_e,background4);

    % Apply Otsu's method to find the global threshold
    global_thresh = graythresh(AF3_f);
    binary_AF3_1 = imbinarize(AF3_f,1.25*global_thresh);

    A_split=A;
    A_split(round(x_int),:) = 0;
    masked_AF3_1 = immultiply(A_split, binary_AF3_1);

    %Remove small objects from the image
    clean_AF3_1 = bwareaopen(masked_AF3_1, 150);

    se=strel('disk',25);
    opened_AF3_1 = imopen(clean_AF3_1, se);

    final_AF3 =opened_AF3_1;
    final_AF3 = immultiply(A, final_AF3);

    valuesF1 = interp2(double(AF1), fit_x, fit_y);
    valuesF2 = interp2(double(AF2), fit_x, fit_y);
    valuesF3 = interp2(double(AF3), fit_x, fit_y);
    valuesP = interp2(double(AP), fit_x, fit_y);

    norm_valuesF1=valuesF1./max(valuesF1);
    norm_valuesF2=valuesF2./max(valuesF2);
    norm_valuesF3=valuesF3./max(valuesF3);
    norm_valuesP=valuesF1./max(valuesP);

    %Measure areas of Glycogen caps and daugther cells

    d1_area=0;
    d2_area=0;
    d1_GS_area=0;
    d2_GS_area=0;

    for m=1:size(A,2)
        for n=1:MidPoint_MidCell-1
            if A(n,m)==1
                d1_area=d1_area+1;
            end
        end
        for n=MidPoint_MidCell+1:size(A,1)
            if A(n,m)==1
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

    d1_area=d1_area/(25*25);
    d2_area=d2_area/(25*25);
    d1_GS_area=d1_GS_area/(25*25);
    d2_GS_area=d2_GS_area/(25*25);

    if d1_area>d2_area && d1_GS_area >0 && d2_GS_area >0
        diff_area=d1_area-d2_area;
        diff_GS_area=d1_GS_area-d2_GS_area;
    elseif d2_area>d1_area && d1_GS_area >0 && d2_GS_area >0
        diff_area=d2_area-d1_area;
        diff_GS_area=d2_GS_area-d1_GS_area;
    else
        diff_area=NaN;
        diff_GS_area=NaN;
    end

    if d1_GS_area >0 && d2_GS_area >0
        Abs_diff_GS_area=sqrt((d1_GS_area-d2_GS_area)^2);
    else
        Abs_diff_GS_area=NaN;
    end


    % Calculate the nucleoid centroid position
    nucleoidlabel_matrix = bwlabel(final_AF3);
    NucleoidRegions=regionprops(nucleoidlabel_matrix,'all');

    if length(NucleoidRegions)==1
        NucleoidMidPointX=NucleoidRegions(1).Centroid(1,1);
        NucleoidMidPointY=NucleoidRegions(1).Centroid(1,2);

    elseif length(NucleoidRegions)>1
        NucleoidMidPointX=(NucleoidRegions(1).Centroid(1,1)+NucleoidRegions(2).Centroid(1,1))/2;
        NucleoidMidPointY=(NucleoidRegions(1).Centroid(1,2)+NucleoidRegions(2).Centroid(1,2))/2;
    end

    CellRegions=regionprops(A,'all');
    NormEuclidian_NucleoidCentroidDistance = sqrt((MidCellX - NucleoidMidPointY)^2+((MidCellY - NucleoidMidPointX)^2))/length(fit_y(1:end));
    Norm_NucleoidCentroidDistance = (MidCellX - NucleoidMidPointY)/length(fit_y(1:end));

    AllCellsDATA(k).CalcCellLength=length(fit_y(1:end));
    AllCellsDATA(k).AP=AP;
    AllCellsDATA(k).AP_mask=mask;
    AllCellsDATA(k).AF2_f=AF2_f;
    AllCellsDATA(k).AF2_mask=final_AF2;
    AllCellsDATA(k).AF3_a=AF3_a;
    AllCellsDATA(k).AF3_f=AF3_f;
    AllCellsDATA(k).AF3_mask=final_AF3;
    AllCellsDATA(k).AreaMicrons=AllCellsDATA(k).Area*(0.07*0.07);
    AllCellsDATA(k).NormEuclidian_NucleoidCentroidDistance=NormEuclidian_NucleoidCentroidDistance;
    AllCellsDATA(k).Norm_NucleoidCentroidDistance=Norm_NucleoidCentroidDistance;
    AllCellsDATA(k).d1_area=d1_area;
    AllCellsDATA(k).d2_area=d2_area;
    AllCellsDATA(k).diff_area=diff_area;
    AllCellsDATA(k).d1_GS_area=d1_GS_area;
    AllCellsDATA(k).d2_GS_area=d2_GS_area;
    AllCellsDATA(k).diff_GS_area=diff_GS_area;
    AllCellsDATA(k).Abs_diff_GS_area=Abs_diff_GS_area;


    for h=1:length(fit_y(1:end))
        AllCellsDATA(k).NormCellLength(h)=h/length(fit_y(1:end));
    end

    fig3=figure(3);
    hold on
    f3=tiledlayout(2,3); % Requires R2019b or later
    nexttile
    hold on;
    imshow(AP*2);
    title('Phase');
    nexttile
    hold on;
    imshow(A);
    scatter(MidCellY,MidCellX,'r','LineWidth',1)
    title('Cell Mask');
    nexttile
    hold on;
    imshow(final_AF3+mask);
    scatter(MidCellY,MidCellX,'r','LineWidth',1)
    scatter(NucleoidMidPointX,NucleoidMidPointY,'b','LineWidth',1)
    title('Segmented Nucleoids');
    nexttile
    hold on;
    imshow(AF3_c);
    title('Raw HU-mCherry');
    nexttile
    hold on;
    imshow(AF2_h);
    title('Raw Glycogen Sensor');
    nexttile
    hold on;
    imshow(final_AF2);
    title('Segmented Glycogen Sensor');
    hold off;
    f3.TileSpacing = 'compact';
    f3.Padding = 'compact';
    fig3_name = sprintf('_Figure3_%d.png', k);
    saveas(gcf,[nameExp,fig3_name]);
    close(fig3);

end

%% Curation of segmented cells

AllCells_Curated=struct([]);
count_curated=0;

for u=1:length(AllCellsDATA)
    if AllCellsDATA(u).diff_GS_area ~= 0 %removes cells in which the calculation of glycogen sensor segmentation failed
        if ~isnan(AllCellsDATA(u).diff_GS_area)
            Ratio=AllCellsDATA(u).Norm_NucleoidCentroidDistance/(AllCellsDATA(u).diff_GS_area*(0.07*0.07)/AllCellsDATA(u).AreaMicrons);
            if Ratio > -50 && Ratio < 50 %filters out extreme outliers due to incorrect calculation of centroid/glycogen area difference
                count_curated=count_curated+1;
                AllCells_Curated(count_curated).CalcCellLength=AllCellsDATA(u).CalcCellLength;
                AllCells_Curated(count_curated).AP=AllCellsDATA(u).AP;
                AllCells_Curated(count_curated).AP_mask=AllCellsDATA(u).AP_mask;
                AllCells_Curated(count_curated).AF2_f=AllCellsDATA(u).AF2_f;
                AllCells_Curated(count_curated).AF2_mask=AllCellsDATA(u).AF2_mask;
                AllCells_Curated(count_curated).AF3_a=AllCellsDATA(u).AF3_a;
                AllCells_Curated(count_curated).AF3_f=AllCellsDATA(u).AF3_f;
                AllCells_Curated(count_curated).AF3_mask=AllCellsDATA(u).AF3_mask;
                AllCells_Curated(count_curated).AreaMicrons=AllCellsDATA(u).AreaMicrons;
                AllCells_Curated(count_curated).Norm_NucleoidCentroidDistance=AllCellsDATA(u).Norm_NucleoidCentroidDistance;
                AllCells_Curated(count_curated).NormAbs_NucleoidCentroidDistance=sqrt((AllCellsDATA(u).Norm_NucleoidCentroidDistance)^2);
                AllCells_Curated(count_curated).NormEuclidian_NucleoidCentroidDistance=AllCellsDATA(u).NormEuclidian_NucleoidCentroidDistance;
                AllCells_Curated(count_curated).d1_area=AllCellsDATA(u).d1_area;
                AllCells_Curated(count_curated).d2_area=AllCellsDATA(u).d2_area;
                AllCells_Curated(count_curated).diff_area=AllCellsDATA(u).diff_area;
                AllCells_Curated(count_curated).d1_GS_area=AllCellsDATA(u).d1_GS_area;
                AllCells_Curated(count_curated).d2_GS_area=AllCellsDATA(u).d2_GS_area;
                AllCells_Curated(count_curated).diff_GS_area=AllCellsDATA(u).diff_GS_area;
                AllCells_Curated(count_curated).Abs_diff_GS_area=AllCellsDATA(u).Abs_diff_GS_area;
                AllCells_Curated(count_curated).Norm_diff_GS_area=AllCellsDATA(u).diff_GS_area*(0.07*0.07)/AllCellsDATA(u).AreaMicrons;
                AllCells_Curated(count_curated).NormAbs_diff_GS_area=AllCellsDATA(u).Abs_diff_GS_area*(0.07*0.07)/AllCellsDATA(u).AreaMicrons;
                AllCells_Curated(count_curated).Ratio=AllCellsDATA(u).Norm_NucleoidCentroidDistance/(AllCellsDATA(u).diff_GS_area*(0.07*0.07)/AllCellsDATA(u).AreaMicrons);

                fig4=figure(4);
                hold on
                f4=tiledlayout(1,4); % Requires R2019b or later
                nexttile
                hold on;
                imshow(AllCellsDATA(u).AP*2);
                title('Phase');
                nexttile
                hold on;
                imshow(AllCellsDATA(u).AP_mask);
                title('Cell Mask');
                nexttile
                hold on;
                imshow(AllCellsDATA(u).AF2_mask);
                title('Segmented GS');
                nexttile
                hold on;
                imshow(AllCellsDATA(u).AF3_mask);
                title('Segmented Nucleoid');
                hold off;
                f4.TileSpacing = 'compact';
                f4.Padding = 'compact';
                fig4_name = sprintf('_Figure4_%d.png', u);
                saveas(gcf,[nameExp,fig4_name]);
                close(fig4);
            end

        end
    end
end



%%
x_gs=[AllCells_Curated.Norm_diff_GS_area]';

x_gs=linspace(-0.15,0.15,1000);
y_gs=x1;

yfit1=1.85822286*x_gs+ 0.01655165; %Fit from PCR analysis, can be replaced by other linear fits

figure(5)
hold on
scatter([AllCells_Curated.Norm_diff_GS_area],[AllCells_Curated.Norm_NucleoidCentroidDistance])
hold off

[CC_wt CP_wt]=corrcoef([AllCells_Curated.Norm_NucleoidCentroidDistance],[AllCells_Curated.Norm_diff_GS_area]);
[rho_wt, pval_wt] = corr([AllCells_Curated.Norm_NucleoidCentroidDistance]',[AllCells_Curated.Norm_diff_GS_area]', 'Type', 'Spearman');
%%
figure(6)
hold on
scatter([AllCells_Curated.NormAbs_diff_GS_area],[AllCells_Curated.NormAbs_NucleoidCentroidDistance])
hold off

[CC_Abs_wt CP_Abs_wt]=corrcoef([AllCells_Curated.NormAbs_NucleoidCentroidDistance],[AllCells_Curated.NormAbs_diff_GS_area]);
[rho_Abs_wt, pval_Abs_wt] = corr([AllCells_Curated.NormAbs_NucleoidCentroidDistance]',[AllCells_Curated.NormAbs_diff_GS_area]', 'Type', 'Spearman');




