function [data,YY,YY_all,n_cells,n_int]=parAtotbetMipZ_2(cellList2)
%% same as original parAtotbetMipZ but some additional outputs added
%% this function records the amount of ParA signals in between MipZ spots
%% requires cellList that has been process such that MipZ peaks have been detected.
%% parA is signal2
%% MipZ is signal1
pixelsize= 0.06421; % micron per pixel
data=zeros(0,7);

%Quantify ParA fluorescence between MipZ focus in a population

YY_all=[];
YY=zeros(1,100);
n_int=0;
n_cells=0;

for frame=1:length(cellList2)
    for cell=1:length(cellList2{frame})
        lst=cellList2{frame}{cell};
        if ~isempty(lst)&& ~isempty(lst.signal1) && ~isempty(lst.steparea)&& ~isempty(lst.signal2)&& ~isempty(lst.spots) %&& length(lst.spots)>2
            
            n_cells=n_cells+1;
            parA=lst.signal1./(lst.steparea*pixelsize^2);
            lngv= pixelsize*cellList2{frame}{cell}.lengthvector;
            lng = lst.length*pixelsize;
            zspot=lst.spots.positions;
            spotNum=length(zspot);
            
            d_pxl=1; % shift in pixels from MipZ peak positions

            for u=2:spotNum
                %               parA1=0;
                dist= lngv(zspot(u))-lngv(zspot(u-1));
                if dist>0
                % quantify ParA signal between MipZ peaks
                    parA1=0.5*(parA(zspot(u-1))+parA(zspot(u)))...
                        +sum(parA(zspot(u-1)+d_pxl:zspot(u)-d_pxl));
                    parA2=sum(parA(zspot(u-1)+2:zspot(u)-2));
                    %parAtot1=sum(parA.*lngv);
                    parAtot1=sum(parA);
                    parAtot2=sum(lst.signal2);
                    parA1_dist=parA(zspot(u-1)+d_pxl:zspot(u)-d_pxl);
                    %                 for i=zspot(u-1):(zspot(u))
                    %                 parA1=parA1+parA(i);
                    %                 end;


                    % addition to calculate average scale-normalized ParA profile between PCs 
                    li_vect=lngv(zspot(u-1)+d_pxl:zspot(u)-d_pxl)-lngv(zspot(u-1));
                    y = integinterp(li_vect/dist,parA1_dist,100); 
                     %y_norm=y*length(parA1_dist)/sum(y);
                     % changed on 08.24.2015, I think it is more correct
                     % way to average profiles
                     y_norm=y/(sum(y)*0.01);

                    %              1    2    3     4      5      6        7
                    data=[data; [parA1 dist lng spotNum parA2 parAtot1 parAtot2]];
                    YY=YY+y_norm;
                    YY_all=[YY_all;y_norm];
                     n_int=n_int+1;
                    end
                
            end
        end
    end
end

figure
hold on;
plot(data(:,2), data(:,1),'og')
xlabel('Distance Between two parS Foci (\mum)','fontsize',14)
title(['Sample size: ' num2str(length(data), '%0.3d')], 'fontsize',14)
ylabel('Integrated ParA Fluorescence (AU)', 'fontsize', 14)
set(gca, 'fontsize', 12)
end