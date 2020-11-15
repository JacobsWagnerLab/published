function [output, T] = autoSelectv4(input, Phase)
%--------------------------------------------------------------------------
%function [output, T] = autoSelectv4(input, Phase)
%@author:  Manuel Campos
%@date:    June, 2014
%@copyright 2013-2014 Yale University
%==========================================================================
%**********output********:
%output:    Structure with cell informations organized into fields
%           (see line 45)
%T:         Structure array containing extra information about mother cell
%           and daughter cells for consistency check of the dataset (T for
%           Test)
%           
%**********input********:
%input:     Structure array as created by dotOutScan.m, which reads data
%           from *.out files (Highthroughput format from Oufti)
%Phase:     Phase images as a 3D uint16 array (nb col x nb line x nb frame)
%           
%==========================================================================
% This function weeds out from the detected cells in the microfluidic
% device that do not meet a number of criteria:
% 1- remain away from the microfluidic chamber exits by the end of the life
% span (at division).
% 2- positive growth rate
% 3- have one ancestor which divided at a "reasonable" site
% 4- have 2 descendants exactly resulting from a division at a "reasonable"
% location along the cell body (away from the pole)
%-------------------------------------------------------------------------- 

% In order to avoid lengthy lines of code, split input fields into several
% variables:
frameArray = input.frameArray;
ancArray = input.ancArray;
progArray = input.progArray;
meshArray = input.meshArray;
cellIdArray = input.cellIdArray;
polarityArray = input.polarityArray;

imgSize=size(Phase,1);

% Get the number of cells (one cell has one ID)
uniqId=unique(cellIdArray);

counter=1;

output(1)=struct('id',[],'ancestor',[],'progeny',[],'polarity',[],'frames',[],'meshes',[],'Lb',[],'Ld',[],...
    'rPm',[],'rpProf',[],'cProf',[]);
% Attribute space for output structure;
testV = nan(length(uniqId),5);
T(1) = struct('Lb',[],'Ld',[],'L2',[],'Lm',[]);

w = waitbar(0, 'Extrapolating dimensions from fits, please wait...');
hw = findobj(w, 'Type', 'Patch');
set(hw, 'EdgeColor', [0.8 .2 0], 'FaceColor', [.8 .2 0])


for cc = 1:length(uniqId) 
    disp(num2str(cc));
%     wantedCell=IDs==uniqId(cc);
%     if sum(wantedCell)>0
    % Initialize & re-initialize tmpSt at each loop
    tmpSt = struct('id', [], 'ancestor', [], 'progeny', [], 'frames', [],...
        'meshes', [], 'lengthes', [], 'extPos', [], 'polarity', []);
    testV(cc,1) = -1;
    ix = find(cellIdArray==uniqId(cc));
    % A minimum of 2 frames in which a given cell is detected is required
    % for further processing of this cell
    if (length(ix)>2) %&& str2double(frameArray{ix(1)})<2000
%         disp(num2str(cc));
        TT = length(ix);
        testV(cc,1) = 0;
        if (~isempty(ancArray{ix(1)})) && (~isempty(progArray{ix(end)}))
            testV(cc,1) = 1;%disp(num2str(cc));
            % Fill in temporary structure
            tmpSt.id = uniqId(cc);
            tmpSt.ancestor = str2num(ancArray{ix(1)});
            tmpSt.progeny = str2num(progArray{ix(end)});
%             tmpSt.lengthes=zeros(1,length(ix));
%             tmpSt.extPos=zeros(length(ix),2);
            for span = 1:length(ix)
                % Store the frame id at which the cell was detected as a
                % double array in the temporary structure
                tmpSt.frames(span) = str2double(frameArray{ix(span)});
                % Store the cell meshes information for the frames at which
                % the cell was detected as a cell array (nb frames x 1)
                M = str2num(meshArray{ix(span), 1})';
                tmpSt.meshes{span} = M;
                % If the size of the mesh is more than 4 points, store the
                % outmost vertices positions and the calculte cell length
                if size(tmpSt.meshes{span}, 1)>4
                    tmpSt.extPos(span,:)=[min([M(:,2);(M(:,4))]),max([M(:,2);(M(:,4))])];
                    steplength = edist(M(2:end,1)+M(2:end,3),M(2:end,2)+M(2:end,4),...
                               M(1:end-1,1)+M(1:end-1,3),M(1:end-1,2)+M(1:end-1,4))/2;
                    tmpSt.lengthes(span)=sum(steplength);
                end
            end
            % Store the information about the position of the cell poles
            % for further testing (below)
            poles=[tmpSt.meshes{2}(1,1), tmpSt.meshes{2}(1,2) ;...
                tmpSt.meshes{2}(end,1), tmpSt.meshes{2}(end,2)];
            
   % The rest of the code is only there to test the criteria 2-4
   % and store the cell's information if all tests are passed
            
            % Test for good cell growth (average growth over the last
            % 99 frames or less must be positive)
            rates=tmpSt.lengthes(2:end)-tmpSt.lengthes(1:end-1);
            if TT>99
                testGrowth=mean(rates(end-98:end))>0;
            else
                testGrowth=mean(rates)>0;
            end
            
            % Test for aberrant size shifts & cell detection at every frame
            % Initialize the test result array
            testSiz = ones(1, 5);
            % Test for aberrant growth by testing whether a cell grew by
            % more than 8 pixels between 2 frames (typically a few seconds)
            % a 8 pixels length shift in a few seconds is a lot!
            testSiz(1)=sum(rates>8); 
            normSiz=(tmpSt.lengthes-tmpSt.lengthes(1))./(tmpSt.lengthes(end)...
                -tmpSt.lengthes(1));
            % Test for aberrant cell length peaks at more than 130% of
            % the final cell length
            testSiz(2) = sum(normSiz > 1.3);
            % Test for aberrant cell length lows at less than -30% of
            % the final cell length
            testSiz(3) = sum(normSiz<-0.3);
            % Check for missing frames in the cell lifetime (non-continuous
            % array of frame ids)
            shiftT = tmpSt.frames(2:end)-tmpSt.frames(1:end-1);
            testSiz(4)=sum(shiftT>1);
            % Test for an aberrant initial cell length (16 pixels id really
            % small, with a pixel size tuypically around 64nm)
            testSiz(5) = tmpSt.lengthes(1)<16;
            % Each cell must pass all 5 tests (sum of all tests equal to 0)
            testS = sum(testSiz)==0;
            
            % Test for closeness to microfluidic chamber edges (bottom & up)
            if TT>46 
                TTedg = 45;
            else
                TTedg = TT-1;
            end
            testBottom = sum(tmpSt.extPos(end-TTedg:end,1)<8)==0;
            testTop = sum(tmpSt.extPos(end-TTedg:end,2)>imgSize-8)==0;
            % Each cell must pass both tests (sum of test results = 2)
            testEdge = (testBottom+testTop)==2;
            
            % Test for a good ancestor
            % That is, the last cell mesh for the ancestor is good. This
            % test also removes all the cells detected at frame 1, for
            % which we do not know the real time of birth.
            % Get the cell id of the ancestor
            idAnc = str2num(ancArray{ix(1)});
            indexAnc = find(cellIdArray==idAnc(end));
            % Initialize the "ancestor test" value to 0
            testAnc=0;
            if ~isempty(indexAnc) 
                % Get the last mesh of the ancestor cell
                mesh=str2num(meshArray{indexAnc(end),1})';
                if length(mesh)>4
                    % Calculate the last steplength array for this ancestor cell
                    steplength = edist(mesh(2:end,1)+mesh(2:end,3),mesh(2:end,2)+mesh(2:end,4),...
                                   mesh(1:end-1,1)+mesh(1:end-1,3),mesh(1:end-1,2)+mesh(1:end-1,4))/2;
                    % Calculate the last width array for this ancestor cell
                    width = sqrt((mesh(:,1)-mesh(:,3)).^2+(mesh(:,2)-mesh(:,4)).^2);
                    % Check that division occured at a reasonnable location
                    % along the cell body. The initial idea was to discard
                    % cells resulting from a division event at less than
                    % one cell radius away from the pole. However, to
                    % reduce the amount to "crappy" cells, I relaxed this
                    % constrain to a minimal distance between the division
                    % site and the ancestor cell's pole corresponding.
                    % The length ratio between the cell and its ancestor
                    % must be larger than the ratio between the ancestor
                    % cell radius and length, and smaller than 1 minus this
                    % r/L ratio.
                    lrR = max(width)/(2*sum(steplength));
                    testAnc=(tmpSt.lengthes(1)/sum(steplength)>1*lrR) &&...
                        (tmpSt.lengthes(1)/sum(steplength)<(1-1*lrR));
                    meshM = mesh;
                    % Store the pole location of the ancestor cell
                    polesM = [mesh(1,1), mesh(1,2) ; mesh(end,1), mesh(end,2)];
                    % Calculate the constriction degree of the ancestor
                   [relPos, ~] = microConstriction(mesh,...
                       Phase(:, :, str2double(frameArray{indexAnc(end)})),...
                       0, 0, lrR);
                end
            end
            
            % Test for a good descendants
            % That is a number of sescendants equal to 2, a good first cell
            % mesh for the descendants and a reasonnable division site location
            % Get the cell id of the descendants
            idProgeny = str2num(progArray{ix(end)});
            % Initialize the "descendant test" values to 0
            testPro = [0 0];
            % Loop through the 2 descendants, if there are exactly 2 of
            % them.
            if length(idProgeny) == 2
                for pro = 1:2 % "pro" for "progeny"
                    indexProgeny = find(cellIdArray==idProgeny(pro));
                    if ~isempty(indexProgeny)
                        % Get the mesh for descendant number 'pro'
                        mesh=str2num(meshArray{indexProgeny(1),1})';
                        % Test for the location of teh division site within
                        % the interval of relative cell length defined by:
                        % [lrR, 1-lrR]
                        if length(mesh)>4;
                            steplength = edist(mesh(2:end,1)+mesh(2:end,3),mesh(2:end,2)+mesh(2:end,4),...
                                       mesh(1:end-1,1)+mesh(1:end-1,3),mesh(1:end-1,2)+mesh(1:end-1,4))/2;
                            width = sqrt((mesh(:,1)-mesh(:,3)).^2+(mesh(:,2)-mesh(:,4)).^2);
                            lrR=max(width)/(2*sum(steplength));
                            testPro(pro)=(sum(steplength)/tmpSt.lengthes(end)>1*lrR) && (sum(steplength)/tmpSt.lengthes(end)<(1-1*lrR));
                            % Store the mesh information for each
                            % decsendant in a 2x1 cell array
                            meshD{pro}=mesh;
                        end
                    end
                end
            end
            % Both descendants must pass the test
            testDesc=sum(testPro)==2;
            
            % If all tests are okay for this cell, store informations about
            % this cell in the output structure array
            if testGrowth && testAnc && testDesc && testEdge && testS
                disp(num2str(cc));
                output(counter).id = uniqId(cc);
                output(counter).ancestor = idAnc;
                output(counter).progeny = idProgeny;
                output(counter).polarity = polarityArray(ix(end));
                % Get the initial cell mesh and calculate the length at
                % birth
                mesh=str2num(meshArray{ix(1),1})';
                steplength = edist(mesh(2:end,1)+mesh(2:end,3),mesh(2:end,2)+mesh(2:end,4),...
                           mesh(1:end-1,1)+mesh(1:end-1,3),mesh(1:end-1,2)+mesh(1:end-1,4))/2;
                T(counter).Lb=sum(steplength)*0.064;
                % Get the last cell mesh and calculate the length at
                % division
                mesh=str2num(meshArray{ix(end),1})';
                steplength = edist(mesh(2:end,1)+mesh(2:end,3),mesh(2:end,2)+mesh(2:end,4),...
                           mesh(1:end-1,1)+mesh(1:end-1,3),mesh(1:end-1,2)+mesh(1:end-1,4))/2;
                T(counter).Ld=sum(steplength)*0.064;output(counter).Ld=sum(steplength)*0.064;
                % Calculate the constriction degree and relatove position
                % of the division site
                width = sqrt((mesh(:,1)-mesh(:,3)).^2+(mesh(:,2)-mesh(:,4)).^2);
                lrR=max(width)/(2*sum(steplength));
                [rPcell,~]=microConstriction(mesh,Phase(:,:,str2double(frameArray{ix(end)})),0,0,lrR);
                output(counter).rpProf=rPcell;
                % Store the information of the mother cell's length at
                % division
                mesh=meshM;
                steplength = edist(mesh(2:end,1)+mesh(2:end,3),mesh(2:end,2)+mesh(2:end,4),...
                           mesh(1:end-1,1)+mesh(1:end-1,3),mesh(1:end-1,2)+mesh(1:end-1,4))/2;
                T(counter).Lm=sum(steplength)*0.064;
                % Mother cell length @ division * septum positiom
                % Decides which side of the mother cell is to be considered
                dP1=pdist([polesM(1,:);poles],'euclidean');
                dP2=pdist([polesM(2,:);poles],'euclidean');
                minP1=min(dP1(1:2));minP2=min(dP2(1:2));
                [~,polId]=min([minP1,minP2]);
                if polId==2
                    relPos=1-relPos;
                end
                % Store the length at birth as the half mother cell
                % corresponding to this cell. This 'trick' prevents from
                % having small but significant bias in cell length at
                % birth due to overlapping daughter cells (~1/2 of a pixel)
                T(counter).LbDC = sum(steplength)*0.064*relPos;
                output(counter).Lb = sum(steplength)*0.064*relPos;
                % Store the relative position of the division site of the
                % mother cell
                T(counter).rP = relPos;
                output(counter).rPm = relPos;
                % Coming back to the descendants:
                mesh = meshD{1};
                steplength1 = edist(mesh(2:end,1)+mesh(2:end,3),mesh(2:end,2)+mesh(2:end,4),...
                           mesh(1:end-1,1)+mesh(1:end-1,3),mesh(1:end-1,2)+mesh(1:end-1,4))/2;
                mesh=meshD{2};
                steplength2 = edist(mesh(2:end,1)+mesh(2:end,3),mesh(2:end,2)+mesh(2:end,4),...
                           mesh(1:end-1,1)+mesh(1:end-1,3),mesh(1:end-1,2)+mesh(1:end-1,4))/2;
                % Save the sum of the lengthes of both descendants for
                % comparison with the length at birth of this cell
                T(counter).L2 = sum(steplength1)*0.064+sum(steplength2)*0.064;
                
% This part can be used to run the function in a parallelized manner
%                 % Prepare parfor for actual cell constriction profile
%                 frAr=str2double(frameArray(ix));meAr=meshArray(ix);
%                 imgPar=Phase(:,:,frAr);%huPar=Signal1(:,:,frAr);
%                 Mpar=cell(1,length(ix));
%                 rpProf=zeros(1,length(ix));cProf=zeros(1,length(ix));
% %                 rpProf2=zeros(1,length(ix));cProf2=zeros(1,length(ix));
%                 parfor span=1:length(ix)
%                     Mpar{span}=str2num(meAr{span})';
%                     [rpProf(span),cProf(span)]=microConstriction(str2num(meAr{span})',imgPar(:,:,span),1,rPcell,lrR);
% %                     [rpProf2(span),cProf2(span),~]=sigConstriction(str2num(meAr{span})',huPar(:,:,span),1,rPcell,lrR);
%                 end
%                 output(counter).frames=frAr;
%                 output(counter).meshes=Mpar;
%                 output(counter).rpProf=rpProf;
%                 output(counter).cProf=cProf;
% %                 output(counter).rpProfAnc=rpProf2;
% % %                 output(counter).cProfAnc=cProf2;
                
                for span=1:length(ix)
                    output(counter).frames(span)=str2double(frameArray{ix(span)});
                    output(counter).meshes{span}=str2num(meshArray{ix(span),1})';
                end
                
                counter = counter+1;
                testV(cc,:) = [uniqId(cc),2,2,2,testEdge];
            else
                testV(cc,:) = [uniqId(cc),testGrowth,testAnc,testDesc,testEdge];
            end
        end
    end
    waitbar(cc/length(uniqId));
%     end
end
close(w);

end

function d=edist(x1,y1,x2,y2)
    % complementary to "getextradata", computes the length between 2 points
    d = sqrt((x2-x1).^2+(y2-y1).^2);
end

function [relPos,DC]=microConstriction(mesh,img,flagpos,divPos,lrRatio)
%--------------------------------------------------------------------------
%function [relPos, DC] = microConstriction(mesh, img, flagpos, divPos, lrRatio)
%@author:  Manuel Campos
%@origin:  Inspired by the corresponding Oufti in-built function
%@date:    June, 2014
%@copyright 2013-2014 Yale University
%==========================================================================
%**********output********:
%mesh:      Structure with cell informations organized into fields that are
%           sorted similarly
%img:       Array
%flagpos:   Array
%divPos:    Array
%lrRatio:   Array
%           
%**********input********:
%relPos:    Structure array as created by dotOutScan.m, which reads data
%           from *.out files (Highthroughput format from Oufti)
%DC:        Phase images as a 3D uint16 array (nb col x nb line x nb frame)
%           
%==========================================================================
% This function calculates the constriction degree and the relative
% position of the division site from the cell mesh and the the associated
% phase constrast image. This is the same approach as in Oufti except that
% the input data structure is not constrained by the cellList format.
% The additional input parameters allows the user to direct the search of
% the division site around a specifi position (+/-  10% cell length) by
% setting flagPos to true and specifying the relative position to target
% (divPos). Otherwise, the division site is searched within a cell length
% fraction [lrRatio, 1-lrRatio] which is typically set to the ratio of the
% cell radius over cell length. This is a way to exclude the cell poles
% from the constriction site analysis.
%--------------------------------------------------------------------------

% Calculate a bounding box for the cell from its cell mesh
mBox=round([min([mesh(:,1);(mesh(:,3))])-25,max([mesh(:,1);(mesh(:,3))])+25;...
    min([mesh(:,2);(mesh(:,4))])-25,max([mesh(:,2);(mesh(:,4))])+25]);
% Check that the bounding box is contained in the image
% That is, values above 1
checkZ=mBox<1;
mBox(checkZ)=1;
% and below the image length and width
checkZ=mBox(1,:)>size(img,2);
mBox(1,checkZ)=size(img,2);
checkZ=mBox(2,:)>size(img,1);
mBox(2,checkZ)=size(img,1);
% Extract the signal profile fro the cell in the pahse constrast image
% using the getOneSignalC.m function
img=im2double(img);
img=max(max(img))-img;
signalProfile = getOneSignalC(mesh,mBox,img,1);
%Smooth twice (as in the original function)
signalProfile = 0.5*signalProfile + 0.25*(signalProfile([1 1:end-1])+signalProfile([2:end end]));
signalProfile = 0.5*signalProfile + 0.25*(signalProfile([1 1:end-1])+signalProfile([2:end end]));

sortProfile=sort(signalProfile,'descend');
maxProfile=mean(sortProfile(1:ceil(length(sortProfile)/3)));
if flagpos
    lowBb=round(length(sortProfile)*(divPos-0.1));
    upBb=round(length(sortProfile)*(divPos+0.1));
else
    lowBb=ceil(length(sortProfile)*1*lrRatio);
    upBb=floor(length(sortProfile)*(1-1*lrRatio));
end
% Get the location and the value of the lowest intensity point in the
% profile
[minProfile,pos]=min(signalProfile(lowBb:upBb));
pos = pos+lowBb-1;
relPos = pos/length(signalProfile);
DC = (maxProfile-minProfile)/maxProfile;
end
