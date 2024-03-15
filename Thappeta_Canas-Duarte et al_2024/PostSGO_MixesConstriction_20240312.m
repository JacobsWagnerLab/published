function PostSGO_MixesConstriction_20240312()
% Author: Silvia, February 05 2024
% Importing segmentation results after SGO and constriction analysis

%%% PARAMETERS %%%
clc
clear all
close all

Dir = uigetdir();% get the directory containing all the files
cd(Dir);
nameExp='Mix2R3_Transition';

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
    %for h=1:2
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
Pre_ConstrictionDATA=struct([]);
cellcount=0;

for pos=1:length(DATAAll)

    for cell=1:length(DATAAll(pos).Cell)
        ShapeF=(DATAAll(pos).statsMask(cell).ConvexArea/DATAAll(pos).statsMask(cell).Area);

        if DATAAll(pos).statsMask(cell).MajorAxisLength >50 && DATAAll(pos).statsMask(cell).Circularity > 0.5 && ShapeF < 1.2 && DATAAll(pos).statsMask(cell).MinorAxisLength >12 && DATAAll(pos).statsMask(cell).MinorAxisLength <19
            %Values used for this filtering were determined from plotting the histrograms of DATAAll

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

%ConstrictionAnalysis=struct([]); For rerunning from when stopped
constricting_cell=0;
%constricting_cell=890; If stopped during curation, value of last curated constricting cell can be placed here for resuming from stop

for i=1:length(Pre_ConstrictionDATA)

    clear A P mask
    close all

    A = imresize(Pre_ConstrictionDATA(i).Mask,25,'lanczos3');
    P = imresize(Pre_ConstrictionDATA(i).Phase,25,'lanczos3');
    mask = boundarymask(A);

    figure(1);
    imshow(labeloverlay(P*10,mask, 'ColorMap', 'spring', 'Transparency', 0.45));
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
nameExp='Mix2R4Transition';

for k=1:length(ConstrictionAnalysis)

    clear A A1 AF1 AF2 AF3 AP mask rows cols y x1 x2 LengthUniqueY1 UniqueY1 distance max_distance
    clear equidistant_y equidistant_x p fit_y fit_x rows2 row_idx LeftHalf RightHalf DistanceLeft DistanceRight
    clear uniqueRightHalf uniqueLeftHalf invLeft_idx2 invLeft_idx
    close all

    A1 = imresize(ConstrictionAnalysis(k).Mask,25,'lanczos3');
    A = imrotate(A1,ConstrictionAnalysis(k).Orientation+90);
    AF1=imresize(ConstrictionAnalysis(k).Fluor1,25,'lanczos3');
    AF1 = imrotate(AF1,ConstrictionAnalysis(k).Orientation+90);
    AF2=imresize(ConstrictionAnalysis(k).Fluor2,25,'lanczos3');
    AF2 = imrotate(AF2,ConstrictionAnalysis(k).Orientation+90);
    AF3=imresize(ConstrictionAnalysis(k).Fluor3,25,'lanczos3');
    AF3 = imrotate(AF3,ConstrictionAnalysis(k).Orientation+90);
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

    figure(1)
    hold on
    %scatter(fit_x,fit_y)
    %scatter(cols,rows)
    scatter(RightHalf(1,:),RightHalf(2,:))
    scatter(LeftHalf(1,:),LeftHalf(2,:))
    % scatter(uniqueDistanceLeft(2,:),uniqueDistanceLeft(1,:))
    % scatter(uniqueDistanceRight(2,:),uniqueDistanceRight(1,:))
    hold off

    %Display the label and the line
    %Also saves figure2 for checking accurate cell mapping

    figure(2)
    fig2=figure(2);
    imshow(AP*10);
    hold on;
    plot(fit_x, fit_y, 'b', 'LineWidth', 2);
    plot(y_divplane,x_divplane, 'r', 'LineWidth', 2);
    scatter(y_int,x_int,'y','LineWidth',3)
    hold off;
    fig2_name = sprintf('_Figure2_%d.png', k);
    saveas(gcf,[nameExp,fig2_name]);
    close(fig2);

    %Segmentation of nucleoid

    AF3_1=AF3;
    v0 = max(AF3_1(:));
    if v0 > 0.975*65535
        AF3_1=AF3_1*0.75;
    elseif v0 < 0.475*65535
        AF3_1=AF3_1*2;
    elseif v0 < 0.275*65535
        AF3_1=AF3_1*3;
    end

    %Apply background subtraction
    background = imopen(AF3_1, strel('disk', 300));
    AF3_a = AF3_1 - background;

    AF3_c = locallapfilt(AF3_a, 0.4, 0.5);

    v2 = max(AF3_c(:));
    if v2 > 0.975*65535
        AF3_c=AF3_c*0.75;
    elseif v2 < 0.475*65535
        AF3_c=AF3_c*2;
    elseif v2 < 0.275*65535
        AF3_c=AF3_c*3;
    end

    AF3_d = medfilt2(AF3_c, [25 25]);
    AF3_e = (0.85*AF3_c) + (AF3_d*0.15);

    v3 = max(AF3_e(:));
    if v3 > 0.975*65535
        AF3_e=AF3_e*0.75;
    elseif v3 < 0.475*65535
        AF3_e=AF3_e*2;
    elseif v3 < 0.275*65535
        AF3_e=AF3_e*3;
    end
    background2 = imopen(AF3_e, strel('disk', 300));
    AF3_f = AF3_e - background2;

    % Apply Otsu's method to find the global threshold
    global_thresh = graythresh(AF3_f);
    binary_AF3_1 = imbinarize(AF3_f,1.25*global_thresh);

    A_split=A;
    A_split(round(x_int),:) = 0;
    masked_AF3_1 = immultiply(A_split, binary_AF3_1);

    %Remove small objects from the image
    clean_AF3_1 = bwareaopen(masked_AF3_1, 15000); 

    se=strel('disk',25);
    opened_AF3_1 = imopen(clean_AF3_1, se);

    final_AF3=opened_AF3_1;
    final_AF3 = immultiply(A, final_AF3);


    %Display the original image and the segmented image side by side

    fig3=figure(3);
    hold on
    f3=tiledlayout(1,3); % Requires R2019b or later
    nexttile
    hold on;
    imshow(AP*2.5);
    plot(fit_x, fit_y, 'b', 'LineWidth', 1);
    plot(y_divplane,x_divplane, 'r', 'LineWidth', 1);
    scatter(y_int,x_int,'y','LineWidth',1)
    title('Original Image');
    nexttile
    hold on;
    imshow(AF3_1);
    title('Raw image');
    nexttile
    hold on;
    imshow(final_AF3);
    plot(fit_x, fit_y, 'b', 'LineWidth', 1);
    plot(y_divplane,x_divplane, 'r', 'LineWidth', 1);
    scatter(y_int,x_int,'y','LineWidth',1)
    title('Segmented Image');
    hold off;
    f3.TileSpacing = 'compact';
    f3.Padding = 'compact';
    fig3_name = sprintf('_Figure3_%d.png', k);
    saveas(gcf,[nameExp,fig3_name]);
    close(fig3);

    % Convert the binary image to a label matrix
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
    MidPoint_MidCell=round(length(fit_y)/2);
    MidPoint_MidCellX=fit_y(MidPoint_MidCell);
    MidPoint_MidCellY=fit_x(MidPoint_MidCell);
    CFP_CellRegions=regionprops(A,AF1,'all');
    YFP_CellRegions=regionprops(A,AF2,'all');

    NormAbs_NucleoidCentroidDistance = sqrt((MidPoint_MidCellX - NucleoidMidPointY)^2)/length(fit_y(1:end));
    NormAbs_ConstrictionSiteDistance = sqrt((MidPoint_MidCellX - x_int)^2)/length(fit_y(1:end));

    Norm_NucleoidCentroidDistance = (MidPoint_MidCellX - NucleoidMidPointY)/length(fit_y(1:end));
    Norm_ConstrictionSiteDistance = (MidPoint_MidCellX - x_int)/length(fit_y(1:end));

    fig4=figure(4);
    hold on
    f4=tiledlayout(2,3); % Requires R2019b or later
    nexttile
    hold on;
    imshow(AP*2);
    title('Phase');
    nexttile
    hold on;
    imshow(A);
    scatter(y_int,x_int,'b','LineWidth',1)
    scatter(MidPoint_MidCellY,MidPoint_MidCellX,'r','LineWidth',1)
    title('Cell Mask');
    nexttile
    hold on;
    imshow(final_AF3+mask);
    scatter(MidPoint_MidCellY,MidPoint_MidCellX,'r','LineWidth',1)
    scatter(NucleoidMidPointX,NucleoidMidPointY,'b','LineWidth',1)
    title('Segmented Nucleoids');
    nexttile
    hold on;
    imshow(AF3_a);
    title('Raw HU-mCherry');
    nexttile
    hold on;
    imshow(AF1*5);
    title('Raw mSCFP3');
    nexttile
    hold on;
    imshow(AF2);
    title('Raw mVenus');
    hold off;
    f4.TileSpacing = 'compact';
    f4.Padding = 'compact';
    fig4_name = sprintf('_Figure4_%d.png', k);
    saveas(gcf,[nameExp,fig4_name]);
    close(fig4);

    ConstrictionAnalysis(k).CalcCellLength=length(fit_y(1:end));

    ConstrictionAnalysis(k).AP=AP;
    ConstrictionAnalysis(k).AP_maks=mask;
    ConstrictionAnalysis(k).AF3_a=AF3_a;
    ConstrictionAnalysis(k).AF3_f=AF3_f;
    ConstrictionAnalysis(k).AF3_maks=final_AF3;
    ConstrictionAnalysis(k).AreaMicrons=ConstrictionAnalysis(k).Area*(0.07*0.07);
    ConstrictionAnalysis(k).CFP_meanIntensity=CFP_CellRegions(1).MeanIntensity;
    ConstrictionAnalysis(k).YFP_meanIntensity=YFP_CellRegions(1).MeanIntensity;
    ConstrictionAnalysis(k).NormAbs_NucleoidCentroidDistance=NormAbs_NucleoidCentroidDistance;
    ConstrictionAnalysis(k).NormAbs_ConstrictionSiteDistance=NormAbs_ConstrictionSiteDistance;
    ConstrictionAnalysis(k).Norm_NucleoidCentroidDistance=Norm_NucleoidCentroidDistance;
    ConstrictionAnalysis(k).Norm_ConstrictionSiteDistance=Norm_ConstrictionSiteDistance;


    for h=1:length(fit_y(1:end))
        ConstrictionAnalysis(k).NormCellLength(h)=h/length(fit_y(1:end));
    end

end
%%
save([nameExp,'_ConstrictionAnalysis_Curatedv3.mat'],'ConstrictionAnalysis','-v7.3')

%% Classifies cells into either WT or glg deletion cells using the appropiate fluorescent marker (YFP and CFP)
% Values used for classification were determined by using the histograms of each fluorescent channel's intensity

M2R3_glg_Cells=struct([]); %naming here indicates the mix (either 1 or 2) and the replica 
M2R3_WT_Cells=struct([]);
count_glg=0;
count_WT=0;

for i=1:length(ConstrictionAnalysis)
    if ConstrictionAnalysis(i).AreaMicrons < 4.25 && ConstrictionAnalysis(i).CFP_meanIntensity < 3000 && ConstrictionAnalysis(i).YFP_meanIntensity > 9000
        count_glg=count_glg+1;
        M2R3_glg_Cells(count_glg).Position=ConstrictionAnalysis(i).Position;
        M2R3_glg_Cells(count_glg).CellID=ConstrictionAnalysis(i).CellID;
        M2R3_glg_Cells(count_glg).Area=ConstrictionAnalysis(i).Area;
        M2R3_glg_Cells(count_glg).AreaMicrons=ConstrictionAnalysis(i).AreaMicrons;
        M2R3_glg_Cells(count_glg).Norm_ConstrictionSiteDistance=ConstrictionAnalysis(i).Norm_ConstrictionSiteDistance;
        M2R3_glg_Cells(count_glg).Norm_NucleoidCentroidDistance=ConstrictionAnalysis(i).Norm_NucleoidCentroidDistance;
        M2R3_glg_Cells(count_glg).NormAbs_ConstrictionSiteDistance=ConstrictionAnalysis(i).NormAbs_ConstrictionSiteDistance;
        M2R3_glg_Cells(count_glg).NormAbs_NucleoidCentroidDistance=ConstrictionAnalysis(i).NormAbs_NucleoidCentroidDistance;
        M2R3_glg_Cells(count_glg).Number1=i;
        M2R3_glg_Cells(count_glg).Number2=count_glg;

    elseif ConstrictionAnalysis(i).AreaMicrons < 4.25 && ConstrictionAnalysis(i).CFP_meanIntensity > 4000 && ConstrictionAnalysis(i).CFP_meanIntensity < 8500 && ConstrictionAnalysis(i).YFP_meanIntensity < 9000
        count_WT=count_WT+1;
        M2R3_WT_Cells(count_WT).Position=ConstrictionAnalysis(i).Position;
        M2R3_WT_Cells(count_WT).CellID=ConstrictionAnalysis(i).CellID;
        M2R3_WT_Cells(count_WT).Area=ConstrictionAnalysis(i).Area;
        M2R3_WT_Cells(count_WT).AreaMicrons=ConstrictionAnalysis(i).AreaMicrons;
        M2R3_WT_Cells(count_WT).Norm_ConstrictionSiteDistance=ConstrictionAnalysis(i).Norm_ConstrictionSiteDistance;
        M2R3_WT_Cells(count_WT).Norm_NucleoidCentroidDistance=ConstrictionAnalysis(i).Norm_NucleoidCentroidDistance;
        M2R3_WT_Cells(count_WT).NormAbs_ConstrictionSiteDistance=ConstrictionAnalysis(i).NormAbs_ConstrictionSiteDistance;
        M2R3_WT_Cells(count_WT).NormAbs_NucleoidCentroidDistance=ConstrictionAnalysis(i).NormAbs_NucleoidCentroidDistance;
        M2R3_WT_Cells(count_WT).Number1=i;
        M2R3_WT_Cells(count_WT).Number2=count_WT;
        M2R3_WT_Cells(count_WT).Ratio=ConstrictionAnalysis(i).Norm_ConstrictionSiteDistance/ConstrictionAnalysis(i).Norm_NucleoidCentroidDistance;

    end
end

figure(5)
hold on
histogram([M2R3_WT_Cells.Norm_NucleoidCentroidDistance],'Normalization','probability')
histogram([M2R3_glg_Cells.Norm_NucleoidCentroidDistance],'Normalization','probability')

hold off

figure(6)
hold on

ksdensity([M2R3_glg_Cells.Norm_ConstrictionSiteDistance])
ksdensity([M2R3_WT_Cells.Norm_ConstrictionSiteDistance])

hold off

figure(7)
hold on
scatter([M2R3_WT_Cells.Norm_ConstrictionSiteDistance],[M2R3_WT_Cells.Norm_NucleoidCentroidDistance])
scatter([M2R3_glg_Cells.Norm_ConstrictionSiteDistance],[M2R3_glg_Cells.Norm_NucleoidCentroidDistance])
hold off

[CC_wt_m2r3 CP_wt_m2r3]=corrcoef([M2R3_WT_Cells.Norm_ConstrictionSiteDistance],[M2R3_WT_Cells.Norm_NucleoidCentroidDistance]);
[rho_wt_m2r3, pval_wt_m2r3] = corr([M2R3_WT_Cells.Norm_ConstrictionSiteDistance]',[M2R3_WT_Cells.Norm_NucleoidCentroidDistance]', 'Type', 'Spearman');

[CC_glg_m2r3 CP_glg_m2r3]=corrcoef([M2R3_glg_Cells.Norm_ConstrictionSiteDistance],[M2R3_glg_Cells.Norm_NucleoidCentroidDistance]);
[rho_glg_m2r3, pval_glg_m2r3] = corr([M2R3_glg_Cells.Norm_ConstrictionSiteDistance]',[M2R3_glg_Cells.Norm_NucleoidCentroidDistance]', 'Type', 'Spearman');


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
fig4_name = sprintf('_Figure4_%d.png', count_cfp);
saveas(gcf,[nameExp,fig4_name]);
close(fig4);

%%
x_1=[M2R3_glg_Cells.diff_GS_area]';

x1=linspace(0,1,1000);
y1=x1;

yfit1=0.99754119*x_1 + 0.12870224; %From PCR analysis, can be replaced by values from any other lineat fit used

figure('Name','5','WindowState','maximized');
hold on
q5=tiledlayout(1,1); % Requires R2019b or later

nexttile
hold on;

scatter([M2R3_glg_Cells.diff_GS_area],[M2R3_glg_Cells.diff_area]);

plot(x_1,yfit1,'k-.');
set(gca,'FontSize',12)
set(gca,'FontName','Arial')
ylabel('\DeltaArea','FontSize',12)
xlabel('\Delta Glycogen sensor Area','FontSize',12)
xlim([-0.05 1]);
ylim([0 1]);
hold off;






