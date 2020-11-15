function [data,S,adjProg,adjAnc] = sizeExtract(input,dt,dx)
%--------------------------------------------------------------------------
%function [data,S,adjProg,adjAnc] = sizeExtract(input,dt,dx)
%@author:  Manuel Campos
%@date:    May, 2014
%@copyright 2013-2014 Yale University
%==========================================================================
%**********intput********:
% From the result of autoSelect(v4) (input) and the time interval dt
% between frames in sec sec and the pixel size dx in um,
% sizeExtract generates:
%
%**********output********:
% A structure data with 12 fields:
% 1. cell id [double]
% 2. ancestor id [double]
% 3. progenies ids [double]
% 4. frames at which the cell is detected [cell array]
% 5. polarity [double]    6.cell length in um at each frame [cell array]
% 7. cell width at each frame [cell array]
% 8. Instantaneous growth rates between frames um/min [cell array]
% 9. Growth rate from instantaneous growth rates in um/min[double]
% 10. Length at birth in um from initial mesh [double]
% 11. length at division in um from last mesh [double]
% 12. Cell cycle time in min [double]
%-
% An array S for data extrapolated from fitting growth curves with 8 columns:
% 1. Lb in um
% 2. growth rate in um/min
% 3. Cell cycle time in min
% 4. Ld in um
% 5. deltaL in um
% 6. confidence interval for growth rate
% 7. confidence interval for Lb
% 8. mean instantaneous growth rate in um/min
%-
% Two adjacency matrices, one for descendants adjProg, one for ancestors
% adjAnc so as to reconstruct genealogic trees.

% idT=cat(1,input.id);
% maxId=max(idT);
progT=cat(1,input.progeny);
maxProg=max(progT);
adjProg=zeros(maxProg,maxProg);
adjAnc=zeros(maxProg,maxProg);
IDs=cat(1,input.id);

w1 = waitbar(0, 'Computing dimensions from meshes, please wait...');
for cc=1:length(input)
%    disp(num2str(cc));
    data(cc).id=input(cc).id;
    data(cc).ancestor=input(cc).ancestor;
    data(cc).progeny=input(cc).progeny;
    data(cc).frames=input(cc).frames;
    data(cc).polarity=input(cc).polarity;
    for ii=1:length(input(cc).meshes)
        % Length
%         disp(num2str(ii));
        mesh=input(cc).meshes{ii};
        lng = size(mesh,1)-1;
        steplength = edist(mesh(2:end,1)+mesh(2:end,3),mesh(2:end,2)+mesh(2:end,4),...
                           mesh(1:end-1,1)+mesh(1:end-1,3),mesh(1:end-1,2)+mesh(1:end-1,4))/2;
        data(cc).length(ii)=sum(steplength)*dx;
%         % area
%         steparea = zeros(lng,1);
%         for i=1:lng, steparea(i)=polyarea([mesh(i:i+1,1);mesh(i+1:-1:i,3)],[mesh(i:i+1,2);mesh(i+1:-1:i,4)]); end
%         data(cc).area(ii) = sum(steparea)*dx*dx;
        % Width
        width = sort(sqrt((mesh(:,1)-mesh(:,3)).^2+(mesh(:,2)-mesh(:,4)).^2),'descend').*dx;
        data(cc).width(ii)=mean(width(1:floor(length(width)/3)));
        % Surface area
%         data(cc).surfArea(ii)=sum(pi*steparea*dx*dx);
        % Volume
        d = edist(mesh(:,1),mesh(:,2),mesh(:,3),mesh(:,4));
        stepvolume = (d(1:end-1).*d(2:end) + (d(1:end-1)-d(2:end)).^2/3).*steplength*pi/4;
        data(cc).volume(ii) = sum(stepvolume);
    end
   % Rates of growth
   data(cc).rates=(data(cc).length(2:end)-data(cc).length(1:end-1))./(data(cc).length(1:end-1)*dt/60);
   % Growth rate
   data(cc).gr=nanmean(data(cc).rates);
   % Length at birth & at division
   data(cc).Lb=data(cc).length(1);data(cc).Ld=data(cc).length(end);
   % Cell cycle time
   data(cc).cct=length(data(cc).frames)*dt/60;
   % Division ratio;
   if isfield(input(cc),'rpProf')
       if length(input(cc).rpProf)==1
           data(cc).rP=input(cc).rpProf;
           data(cc).rPm=input(cc).rPm;
       end
   end
   
   % Construct ancestor adjacency matrix
   if ~isempty(input(cc).ancestor)
       if any(IDs==input(cc).ancestor(end))
           adjAnc(input(cc).id,input(cc).ancestor(end))=1;
       end
   end
   % Construct progeny adjacency matrix
   if ~isempty(input(cc).progeny)
       for pro=1:length(input(cc).progeny)
           if any(IDs==input(cc).progeny(pro))
               adjProg(input(cc).id,input(cc).progeny(pro))=1;
           end
           if any(IDs==input(cc).progeny(pro))
               adjProg(input(cc).id,input(cc).progeny(pro))=1;
           end
       end
   end
   waitbar(cc/length(input));
end
close(w1);

w2 = waitbar(0, 'Extrapolating dimensions from fits, please wait...');
hw2=findobj(w2,'Type','Patch');set(hw2,'EdgeColor',[0 1 0],'FaceColor',[0 1 0])
% adjAnc=1;adjProg=1;
S=zeros(length(data),12);%1.Lb 2.gr 3.cct 4.Ld 5.DL 6.ci gr 7.ci log(Lb) 8.rmse
% 9. Xpos at birth 10. Ypos at birth 11.Xpos at division 12, Ypos at
% division
for ii=1:length(data)
    lb=floor(0.5*length(data(ii).frames));%disp(num2str(ii));
    ub=ceil(1*length(data(ii).frames));
    Tim=(data(ii).frames(end)-data(ii).frames(1)+1)'.*dt/60+dt/60;
    Abs=(data(ii).frames(lb:ub)-data(ii).frames(1)+1).*dt/60;
    if size(Abs,1)<size(Abs,2) Abs=Abs'; end
    Ord=(data(ii).length(lb:ub))';
    % Get rid of consecutive identical meshes
    DOrd=Ord(2:end)-Ord(1:end-1);ix0=DOrd~=0;ix1=DOrd~=0 & abs(DOrd)<0.32;
    if sum(ix1)>1
        XX=Abs(ix1);YY=log(Ord(ix1));
    rates=data(ii).rates;ixR=data(ii).length(2:end)>3;
    rates2=rates(ixR);rtinliers=(abs(rates(ixR))<2*std(rates(ixR)));
    rates2=rates2(rtinliers);
    
    if length(XX)>2
    ft = fittype( 'poly1' );    %opts = fitoptions( 'Method', 'LinearLeastSquares' );
%     ft = fittype('a*exp(b*x)+c','dependent',{'y'},'independent',{'x'},'coefficients',{'a', 'b', 'c'});
    [~,~,R]=fit(XX,YY,ft);%,'startpoint',[3 0.025 0.1]);%
    inliers=(abs(R.residuals)<2*std(R.residuals));
    [Lf,gof]=fit(XX(inliers),YY(inliers),ft);%,'startpoint',[3 0.025 0.1]);
    p1=Lf.p1;p2=Lf.p2;%p1=Lf.b;p2=Lf.a;%Lrmse=gof.rmse;
    
    S(ii,1)=exp(p2);%Lb
    S(ii,2)=p1;%gr
    S(ii,3)=Tim+dt/60;%cct
    S(ii,4)=S(ii,1).*exp(p1.*Tim);%Ld
%     ci=confint(Lf);
    S(ii,6)=mean(rates(inliers(1:end-1)));%ci(2,1)-ci(1,1);
    S(ii,7)=gof.rmse;%Lf.c;%
    S(ii,8)=mean(data(ii).rates);%S(ii,6)=sum(DOrd==0);S(ii,7)=sum(DOrd>0.32);
    S(ii,9:12)=[mean([input(ii).meshes{1}(:,1);input(ii).meshes{1}(:,3)]),...
        mean([input(ii).meshes{1}(:,2);input(ii).meshes{1}(:,4)]),...
        mean([input(ii).meshes{end}(:,1);input(ii).meshes{end}(:,3)]),...
        mean([input(ii).meshes{end}(:,2);input(ii).meshes{end}(:,4)]),].*dx;
    else
        S(ii,:)=NaN(1,12);
    end
    end
    waitbar(ii/length(data));
end
S(:,5)=S(:,4)-S(:,1);%DL
close(w2);
end

function d=edist(x1,y1,x2,y2)
    % complementary for "getextradata", computes the length between 2 points
    d=sqrt((x2-x1).^2+(y2-y1).^2);
end

% -helpers
%
% Lb=cat(1,D.Lb);
% Ld=cat(1,D.Ld);
% DL=Ld-Lb;
% gr=cat(1,D.gr);
% cct=cat(1,D.cct);
% aT=gr.*cct;
% pol=cat(1,D.polarity);
% Tb=[];for ii=1:length(D) Tb(ii)=D(ii).frames(1);end
% Tb=Tb';
% Td=[];for ii=1:length(D) Td(ii)=D(ii).frames(end);end
% Td=Td';

% Lb=Lb(ix);Ld=Ld(ix);gr=gr(ix);cct=cct(ix);aT=aT(ix);pol=pol(ix);Tb=Tb(ix);DL=DL(ix);