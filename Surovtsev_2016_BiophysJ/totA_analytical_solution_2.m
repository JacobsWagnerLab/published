% to compare exp data to ananlytical solution to diffusion-redistribution
% model of the ParA grdaient formation
% assumes data1 matrix from Hoong's analysis
%

% assumes data organized as follows:
%         1    2    3     4      5      6        7
% data=[parA1 dist lng spotNum parA2 parAtot1 parAtot2];

n_per_bin=25;

data1=data;

ii=1; YXX=[];
for ii=1:size(data1,1)
  
  li=data1(ii,2);      % inter-PC distance  
  l_cell=data1(ii,3); % cell length
  all_li=data1(data1(:,3)==l_cell,2); % all distances between PCs from that cell
   l2=sum(all_li);
   l3_sum=sum(all_li.^3);
  
 aa= data1(ii,1);   
 a_i= data1(ii,5);    % fluo between PCs  
 a_tot= data1(ii,6); % total fluo
  a_star=a_i/a_tot;        % 
   x=l3_sum/l_cell;
  %x2=l3_sum/li^3/l2/6;
  
 spotNum= data(ii,4);
  y=li^3/a_star/l_cell;
   %x3=l3_sum/6/l_cell;
            %1 2  3  4   5      6     7  8   9      10      11        12
  YXX(ii,:)=[y,x,aa,a_i,a_tot,a_star,li,l2,l3_sum,l_cell,li^3/a_star,spotNum]; 
    
    
end

% Plotting parameters
ms1=12; ms2=9; ms3=6;
fs1=18; fs2=16;fs3=12;
lw1=3; lw2=2; lw3=1;

% figure; hold on
%  tit='l_i^3/a* vs \Sigmal_i^3/l';
%  lgnd={'raw data','fit','mean-binned','fit','median-binned','fit'};
%  ind=(YXX(:,1)<Inf & YXX(:,2)<60); %(YXX(:,5)./YXX(:,10)>15 & YXX(:,5)./YXX(:,10)<50 & YXX(:,2)<17 & YXX(:,1)<40);
% 
% [Xbin1,Ybin1]=bin2_fixN(YXX(ind,2),YXX(ind,1),n_per_bin,1);
%   [Xbin2,Ybin2]=bin2_fixN(YXX(ind,2),YXX(ind,1),n_per_bin,2);
%  [coeff_r,fit_r]=oneFit(YXX(ind,2),YXX(ind,1),'lin_growth');
%    disp('fit to raw data:'); disp(fit_r); a_r= coeff_r(1); b_r=coeff_r(2);
%   [coeff_bin1,fit_bin1]=oneFit(Xbin1,Ybin1,'lin_growth');
%     disp('fit to mean-binned:'); disp(fit_bin1); a_bin1= coeff_bin1(1); b_bin1=coeff_bin1(2);
%   [coeff_bin2,fit_bin2]=oneFit(Xbin2,Ybin2,'lin_growth');
%     disp('fit to median-binned:'); disp(fit_bin2); a_bin2= coeff_bin2(1); b_bin2=coeff_bin2(2);
%  xmax=max(YXX(ind,2));   
% 
% plot (YXX(ind,2),YXX(ind,1),'xb','MarkerSize',ms3,'LineWidth',lw2)
%  plot ([0,xmax],a_r*[0,xmax]+b_r,'-b','LineWidth',lw2)
% plot (Xbin1,Ybin1,'og','MarkerSize',ms2,'LineWidth',lw2)
%  plot ([0,xmax],a_bin1*[0,xmax]+b_bin1,'-g','LineWidth',lw2)
% plot (Xbin2,Ybin2,'or','MarkerSize',ms2,'LineWidth',lw2)
%  plot ([0,xmax],a_bin2*[0,xmax]+b_bin2,'-r','LineWidth',lw2)
%  set(gca, 'FontSize',fs2)
%  %legend(lgnd,'FontSize',14)
%  xlabel('X','FontSize',fs1); ylabel('Y','FontSize',fs1); 
%  title (tit,'FontSize',fs3,'FontName','Cambria','FontAngle','Italic')
%  legend (lgnd, 'fontsize', fs3);
%  %ylim([0,1])


% figure; hold on
figure; hold on
dd=7; %number of points from end to omit for fitting
 tit='scaled ParA profile';
 x_lab='Inter-PC distance'; y_lab='Normalized fluorescence'; 
 lgnd={'mean','fit'};
dx=1/length(YY); XX=dx/2:dx:1;
YY=YY/n_int;
[coeff,fit1,YYfit]=oneFit(XX(dd+1:end-dd)',YY(dd+1:end-dd)','ax(1-x)+b');
 disp('fit to scaled profile:'); disp(fit1); a= coeff(1); b=coeff(2);
plot (XX,YY,'or','MarkerSize',ms2,'LineWidth',lw2);
 plot(XX,a*XX.*(1-XX)+b,'-k','LineWidth',lw2);

title(tit, 'FontSize',fs2)
xlabel (x_lab, 'fontsize', fs1); ylabel (y_lab, 'fontsize', fs1);
set(gca, 'fontsize', fs2)
legend (lgnd, 'fontsize', fs3); 




% to plot exp vs theor. PArA signal between PCs
% corrected for initially lost factor of 1/6
a0=YXX(:,5)./(1+kD*YXX(:,9)/12./YXX(:,10));
figure; plot(YXX(:,4),0.5*kD*a0.*(YXX(:,7).^3)./YXX(:,10)/6,'xk')