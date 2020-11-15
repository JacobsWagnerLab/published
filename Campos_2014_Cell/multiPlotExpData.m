function Out=multiPlotExpData(ExpData,dt,plots)

%---------------------------------------
% to plot multiple distribution and correlation plots from experimental data
% dt = time interval
% ExpData has 2 fields: .AA and .data
%   assumes that AA are Experimental Data are in "Manu's" format: a structure array 
%   with fields
%     id       cellID
%     ancestor 
%     progeny
%     polarity
%     frames
%     meshes
%     Lb
%     Ld
%     rPm
%     rpProf
%     cProf

%   Since AA-data lacks growth rate sinfromation, one need to supply data,
%   which is assumed to be a table with 2nd column containing growth rates organized in the same order as AA array
% plots = type of desired plots, an array with following values: 
%     % Distribution plots
% 1 = 'Growth rate distribution';
% 2 = 'deltaL  distribution';     
% 3 = 'Lb  distribution';     
% 4 = 'Ld  distribution';     
% 5 =  'Division Ratio distribution';
% 6 =  'InterDivision Times distribution';
% 7 =  'InterDivision Times distribution,LogNormal';
%      
%      % Correlation plots
% 12 = 'Scatter: Growth rate vs DeltaL';
% 21 =  'Scatter: DeltaL vs Growth rate';
% 13 = 'Scatter: Growth rate vs Lb';
% 31 = 'Scatter: Lb vs Growth rate';
% 14 =  'Scatter: Growth rate vs Ld';
% 41 =  'Scatter: Ld vs Growth rate';
% 15 =  'Scatter: Growth rate vs DR';
% 51 =  'Scatter: DR vs Growth rate';
% 16 =  'Scatter: Growth rate vs InterDivision Time';
% 61 =  'Scatter: DR vs Growth rate';
% 23 =  'Scatter: DeltaL vs Lb';
% 32 =  'Scatter: Lb vs DeltaL';
% 24 =  'Scatter: DeltaL vs Ld';
% 42 =  'Scatter: Ld vs DeltaL';
% 25 =  'Scatter: DeltaL vs DR';
% 52 =  'Scatter: DR vs DeltaL';
% 26 =  'Scatter: DeltaL vs InterDivision Time';
% 62 =  'Scatter: InterDivisionTime vs DeltaL';
% 34 =  'Scatter: Lb vs Ld';
% 43 =  'Scatter: Ld vs Lb';
% 35 =  'Scatter: Lb vs DR';
% 53 =  'Scatter: DR vs Lb';
% 36 =  'Scatter: Lb vs InterDivision Time';
% 63 =  'Scatter: InterDivision Time vs Lb';
% 45 =  'Scatter: Ld vs DR';
% 54 =  'Scatter: DR vs Ld';
% 46 =  'Scatter: Ld vs InterDivision Time';
% 64 =  'Scatter: InterDivision Time vs Ld';
% 56 =  'Scatter: DR vs InterDivision Time';
% 65 =  'Scatter: InterDivision Time vs DR';
% 
%      % Correlation between generations 
% 111 =  'Scatter: Growth rate Daughter vs Mother';
% 211 =  'Scatter: Growth rate Sister vs Sister';
% 122 =  'Scatter: DeltaL Daughter vs Mother';
% 222 =  'Scatter: DeltaL Sister vs Sister';
% 133 =  'Scatter: Lb Daughter vs Mother';
% 233 =  'Scatter: Lb Sister vs Sister';
% 144 =  'Scatter: Ld Daughter vs Mother';
% 244 =  'Scatter: Ld Sister vs Sister';
% 155 =  'Scatter: DR Daughter vs Mother';
% 255 =  'Scatter: DR Sister vs Sister';
% 166 =  'Scatter: CCT Daughter vs Mother';
% 266 = 'Scatter: CCT Sister vs Sister';
%
%
%
% version: 05.28.2015


c
 data=ExpData.data;

Cells=[AA.id]';
 MotherID=arrayfun(@(x) x.ancestor(end),AA)'; %[AA.ancestor];
 SisterID=[AA.progeny]';
 Polarity=[AA.polarity]';
Lb=[AA.Lb]';
 Ld=[AA.Ld]';
 DeltaL=Ld-Lb;
Alpha=data(:,2);
if isfield(AA,'rP')
 DR=[AA.rP]'; 
else
 DR=[AA.rPm]';   
end
Tb=((arrayfun(@(x) x.frames(1),AA)-1)*dt)';
 Td=((arrayfun(@(x) x.frames(end),AA)-1)*dt)';
 TT=Td-Tb+dt; 




  % Appearance
  n_bins=50; n_bins2=100;
   fs1=18; fs2=16; fs3=14; % font sizes for Labels, tick marks, text, respectively
   lw1=3; lw2=2; lw3=1;       % line width
    ms1=12; ms2=8; ms3=3;
   mrk1='x';mrk2='o';mrk3='.';
    col_m=[0,0,1]; col_m2=[1,0,0]; col_m3=[0,1,0];
   mrk_fit='-';
   col_f='k';

 
 
%-----------------------------------
 % GETTING GENERATIONS CORRELATIONS
 gen_1=1; gen_2=2;   % initial and final generations for correlation plot

% prepare collectors for "per gerenration" data
for kk=gen_1:gen_2 
  Lb_a{kk}=[]; Lb_d{kk}=[]; % T_m{1}=[]; alpha_m{1}=[];
   Lb_s1{kk}=[]; Lb_s2{kk}=[];
  Ld_a{kk}=[]; Ld_d{kk}=[]; % T_m{1}=[]; alpha_m{1}=[];
   Ld_s1{kk}=[]; Ld_s2{kk}=[]; 
  DeltaL_a{kk}=[]; DeltaL_d{kk}=[];
   DeltaL_s1{kk}=[]; DeltaL_s2{kk}=[];
  Alpha_a{kk}=[]; Alpha_d{kk}=[];
   Alpha_s1{kk}=[]; Alpha_s2{kk}=[];
  DR_a{kk}=[]; DR_d{kk}=[];
   DR_s1{kk}=[]; DR_s2{kk}=[];
  TT_a{kk}=[]; TT_d{kk}=[];
   TT_s1{kk}=[]; TT_s2{kk}=[]; 
end

for cell=1:length(Cells)
  mother=Cells(cell);
  gen_t=1;
  while ~isempty(mother)&& gen_t<gen_2
    mother_0=mother;
    for mm=1:length(mother_0)
      mother=[];  
      mother_1=mother_0(mm);  
      ii00=find(MotherID==mother_1,2,'first'); 
      if ~isempty(ii00)
        Lb_d{gen_t}=[Lb_d{gen_t}; Lb(ii00)];
         Lb_a{gen_t}=[Lb_a{gen_t}; Lb(cell)*ones(size(Lb(ii00)))];
        Ld_d{gen_t}=[Ld_d{gen_t}; Ld(ii00)];
         Ld_a{gen_t}=[Ld_a{gen_t}; Ld(cell)*ones(size(Lb(ii00)))]; 
        DeltaL_d{gen_t}=[DeltaL_d{gen_t}; DeltaL(ii00)];
         DeltaL_a{gen_t}=[DeltaL_a{gen_t}; DeltaL(cell)*ones(size(Lb(ii00)))];
        Alpha_d{gen_t}=[Alpha_d{gen_t}; Alpha(ii00)];
         Alpha_a{gen_t}=[Alpha_a{gen_t}; Alpha(cell)*ones(size(Lb(ii00)))];
        DR_d{gen_t}=[DR_d{gen_t}; DR(ii00)];
         DR_a{gen_t}=[DR_a{gen_t}; DR(cell)*ones(size(Lb(ii00)))]; 
        TT_d{gen_t}=[TT_d{gen_t}; TT(ii00)];
         TT_a{gen_t}=[TT_a{gen_t}; TT(cell)*ones(size(Lb(ii00)))];
                 
        mother=[mother;Cells(ii00)];
        if length(ii00)==2
          Lb_s1{gen_t}=[Lb_s1{gen_t}; Lb(ii00(1))];
           Lb_s2{gen_t}=[Lb_s2{gen_t}; Lb(ii00(2))];
          Ld_s1{gen_t}=[Ld_s1{gen_t}; Ld(ii00(1))];
           Ld_s2{gen_t}=[Ld_s2{gen_t}; Ld(ii00(2))];           
          DeltaL_s1{gen_t}=[DeltaL_s1{gen_t}; DeltaL(ii00(1))];
           DeltaL_s2{gen_t}=[DeltaL_s2{gen_t}; DeltaL(ii00(2))];
          Alpha_s1{gen_t}=[Alpha_s1{gen_t}; Alpha(ii00(1))];
           Alpha_s2{gen_t}=[Alpha_s2{gen_t}; Alpha(ii00(2))];
          DR_s1{gen_t}=[DR_s1{gen_t}; DR(ii00(1))];
           DR_s2{gen_t}=[DR_s2{gen_t}; DR(ii00(2))];
           TT_s1{gen_t}=[TT_s1{gen_t}; TT(ii00(1))];
           TT_s2{gen_t}=[TT_s2{gen_t}; TT(ii00(2))];
          %sister1=[sister1;Cells(ii00(1))]; 
          % sister2=[sister2;Cells(ii00(2))]; 
         end
      end
    end
    gen_t=gen_t+1;
  end
    
end


for ii=1:length(plots)
 pp=plots(ii);
 corr_coeff=[];     model='no fit';  
 
 switch pp
     
   case 1 % 'Growth rate distribution';
       pplot='Growth rate distribution';  disp(pplot)
         model='Gauss_normalized';
       [fit_out, fit_results, txt_x, txt_y]=plot_hist(Alpha,n_bins,model);
         change_appearance('\alpha, min^{-1}','Frequency',{'Exp data',[model, ' fit']},num2str(fit_results,3));
    
   case 2 %'deltaL  distribution';     
       pplot='Growth rate distribution';   disp(pplot)
         model='LogNormal';
        [fit_out, fit_results, txt_x, txt_y]=plot_hist(DeltaL,n_bins,model);
         change_appearance('\DeltaL, \mum','Frequency',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         
    case 3 %'Lb  distribution';     
       pplot='Lb distribution';  disp(pplot)
        model='LogNormal';
        [fit_out, fit_results, txt_x, txt_y]=plot_hist(Lb,n_bins,model);
         change_appearance('Lb, \mum','Frequency',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         
   case 4 %'Ld  distribution';     
       pplot='Ld distribution';  disp(pplot)
        model='LogNormal';
       [fit_out, fit_results, txt_x, txt_y]=plot_hist(Ld,n_bins,model);
         change_appearance('Ld, \mum','Frequency',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         
   case 5 % 'Division Ratio distribution';
       pplot='Division Ratio distribution';  disp(pplot)
         model='Gauss_normalized';
        [fit_out, fit_results, txt_x, txt_y]=plot_hist(DR,n_bins,model);
         change_appearance('DR', 'Frequency', {'Exp data',[model, ' fit']}, num2str(fit_results,3));
    
   case 6 % 'InterDivision Times distribution';
       pplot='CCTime distribution';  disp(pplot)
         model='Gauss_normalized';
        [fit_out, fit_results, txt_x, txt_y]=plot_hist(TT,n_bins,model);
         change_appearance('Interdivisiont time, min','Frequency',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         
    case 7 % 'InterDivision Times distribution,LogNormal';
       pplot='CCTime distribution';  disp(pplot)
         model='LogNormal';
        [fit_out, fit_results, txt_x, txt_y]=plot_hist(TT, n_bins, model);
         change_appearance('Interdivisiont time, min','Frequency',{'Exp data', [model, ' fit']}, num2str(fit_results,3)); 
         
     % Correlation plots
    case 12 % 'Scatter: Growth rate vs DeltaL';
       pplot='Growth rate vs DeltaL';  disp(pplot)
         %model='Gauss_normalized';
       [fit_out, fit_results, txt_x, txt_y]=plot_corr(DeltaL, Alpha, n_bins2, 1); 
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('\DeltaL, \mum','\alpha, min^{-1}',{'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
    
    case 21 % 'Scatter: DeltaL vs Growth rate';
       pplot=' DeltaL vs Growth rate';  disp(pplot)
         %model='Gauss_normalized';
        [fit_out, fit_results, txt_x, txt_y]=plot_corr(Alpha, DeltaL, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('\alpha, min^{-1}','\DeltaL, \mum',{'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
         
    case 13 % 'Scatter: Growth rate vs Lb';
       pplot='Growth rate vs Lb';  disp(pplot)
         %model='Gauss_normalized';
      [fit_out, fit_results, txt_x, txt_y]=plot_corr(Lb, Alpha, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Lb, \mum','\alpha, min^{-1}',{'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});    
         
     case 31 % 'Scatter: Lb vs Growth rate';
       pplot='Lb Growth rate ';  disp(pplot)
         %model='Gauss_normalized';
        [fit_out, fit_results, txt_x, txt_y]=plot_corr(Alpha, Lb, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('\alpha, min^{-1}', 'Lb, \mum', {'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});  
       
     case 14 % 'Scatter: Growth rate vs Ld';
       pplot='Growth rate vs Ld';  disp(pplot)
         %model='Gauss_normalized';
      [fit_out, fit_results, txt_x, txt_y]=plot_corr(Ld, Alpha, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Ld, \mum','\alpha, min^{-1}',{'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});    
         
     case 41 % 'Scatter: Ld vs Growth rate';
       pplot='Ld vs Growth rate ';  disp(pplot)
         %model='Gauss_normalized';
        [fit_out, fit_results, txt_x, txt_y]=plot_corr(Alpha, Ld, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('\alpha, min^{-1}', 'Ld, \mum', {'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});     
         
     case 15 % 'Scatter: Growth rate vs DR';
       pplot='Growth rate vs DR';  disp(pplot)
         %model='Gauss_normalized';
      [fit_out, fit_results, txt_x, txt_y]=plot_corr(DR, Alpha, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('DR','\alpha, min^{-1}',{'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});    
         
     case 51 % 'Scatter: DR vs Growth rate';
       pplot='DR vs Growth rate ';  disp(pplot)
         %model='Gauss_normalized';
        [fit_out, fit_results, txt_x, txt_y]=plot_corr(Alpha, DR, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('\alpha, min^{-1}', 'DR', {'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});       
         
     case 16 % 'Scatter: Growth rate vs InterDivision Time';
       pplot='Growth rate vs CCT';  disp(pplot)
         %model='Gauss_normalized';
      [fit_out, fit_results, txt_x, txt_y]=plot_corr(TT, Alpha, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Interdivision time, min','\alpha, min^{-1}',{'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});    
         
     case 61 % 'Scatter: DR vs Growth rate';
       pplot='CCT vs Growth rate ';  disp(pplot)
         %model='Gauss_normalized';
        [fit_out, fit_results, txt_x, txt_y]=plot_corr(Alpha, TT, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('\alpha, min^{-1}', 'Interdivision time, min', {'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});   
         
     case 23 % 'Scatter: DeltaL vs Lb';
       pplot=' DeltaL vs Lb';  disp(pplot)
         %model='Gauss_normalized';
        [fit_out, fit_results, txt_x, txt_y]=plot_corr(Lb,DeltaL, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Lb, \mum','\DeltaL, \mum',{'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
         
    case 32 % 'Scatter: Lb vs DeltaL';
       pplot='Lb vs DeltaL';  disp(pplot)
         %model='Gauss_normalized';
      [fit_out, fit_results, txt_x, txt_y]=plot_corr(DeltaL, Lb, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('\DeltaL, \mum', 'Lb, \mum', {'Single cell Data', ['binned, 1 point=',num2str(n_bins2),' cells']});
         
    case 24 % 'Scatter: DeltaL vs Ld';
       pplot=' DeltaL vs Ld';  disp(pplot)
         %model='Gauss_normalized';
        [fit_out, fit_results, txt_x, txt_y]=plot_corr(Ld, DeltaL, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Ld, \mum','\DeltaL, \mum',{'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
         
    case 42 % 'Scatter: Ld vs DeltaL';
       pplot='Ld vs DeltaL';  disp(pplot)
         %model='Gauss_normalized';
      [fit_out, fit_results, txt_x, txt_y]=plot_corr(DeltaL, Ld, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('\DeltaL, \mum', 'Ld, \mum', {'Single cell Data', ['binned, 1 point=',num2str(n_bins2),' cells']});    
         
    case 25 % 'Scatter: DeltaL vs DR';
       pplot=' DeltaL vs DR';  disp(pplot)
         %model='Gauss_normalized';
        [fit_out, fit_results, txt_x, txt_y]=plot_corr(DR, DeltaL, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('DR','\DeltaL, \mum',{'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
         
    case 52 % 'Scatter: DR vs DeltaL';
       pplot='DR vs DeltaL';  disp(pplot)
         %model='Gauss_normalized';
      [fit_out, fit_results, txt_x, txt_y]=plot_corr(DeltaL, DR, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('\DeltaL, \mum', 'DR, \mum', {'Single cell Data', ['binned, 1 point=',num2str(n_bins2),' cells']}); 
         
    case 26 % 'Scatter: DeltaL vs InterDivision Time';
       pplot=' DeltaL vs CCT';  disp(pplot)
         %model='Gauss_normalized';
        [fit_out, fit_results, txt_x, txt_y]=plot_corr(TT, DeltaL, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Interdivision time, min', '\DeltaL, \mum',{'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
         
    case 62 % 'Scatter: InterDivisionTime vs DeltaL';
       pplot='CCT vs DeltaL';  disp(pplot)
         %model='Gauss_normalized';
      [fit_out, fit_results, txt_x, txt_y]=plot_corr(DeltaL, TT, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('\DeltaL, \mum', 'Intedivision time, min', {'Single cell Data', ['binned, 1 point=',num2str(n_bins2),' cells']}); 
     
    case 34 % 'Scatter: Lb vs Ld';
       pplot=' Lb vs Ld';  disp(pplot)
         %model='Gauss_normalized';
        [fit_out, fit_results, txt_x, txt_y]=plot_corr(Ld, Lb, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Ld, \mum','Lb, \mum',{'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
         
    case 43 % 'Scatter: Ld vs Lb';
       pplot='Ld vs Lb';  disp(pplot)
         %model='Gauss_normalized';
      [fit_out, fit_results, txt_x, txt_y]=plot_corr(Lb, Ld, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Lb, \mum', 'Ld, \mum', {'Single cell Data', ['binned, 1 point=',num2str(n_bins2),' cells']});     
         
    case 35 % 'Scatter: Lb vs DR';
       pplot=' Lb vs DR';  disp(pplot)
         %model='Gauss_normalized';
        [fit_out, fit_results, txt_x, txt_y]=plot_corr(DR, Lb, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('DR','Lb, \mum',{'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
         
    case 53 % 'Scatter: DR vs Lb';
       pplot='DR vs Lb';  disp(pplot)
         %model='Gauss_normalized';
      [fit_out, fit_results, txt_x, txt_y]=plot_corr(Lb, DR, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Lb, \mum', 'DR', {'Single cell Data', ['binned, 1 point=',num2str(n_bins2),' cells']});    

    case 36 % 'Scatter: Lb vs InterDivision Time';
       pplot=' Lb vs CCT';  disp(pplot)
         %model='Gauss_normalized';
        [fit_out, fit_results, txt_x, txt_y]=plot_corr(TT, Lb, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Interdivision time, min','Lb, \mum',{'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
         
    case 63 % 'Scatter: InterDivision Time vs Lb';
       pplot='CCT vs Lb';  disp(pplot)
         %model='Gauss_normalized';
      [fit_out, fit_results, txt_x, txt_y]=plot_corr(Lb, TT, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Lb, \mum', 'Interdivision time, min', {'Single cell Data', ['binned, 1 point=',num2str(n_bins2),' cells']}); 
         
    case 45 % 'Scatter: Ld vs DR';
       pplot=' Ld vs DR';  disp(pplot)
         %model='Gauss_normalized';
        [fit_out, fit_results, txt_x, txt_y]=plot_corr(DR, Ld, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('DR','Ld, \mum',{'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
         
    case 54 % 'Scatter: DR vs Ld';
       pplot='DR vs Lb';  disp(pplot)
         %model='Gauss_normalized';
      [fit_out, fit_results, txt_x, txt_y]=plot_corr(Ld, DR, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Ld, \mum', 'DR', {'Single cell Data', ['binned, 1 point=',num2str(n_bins2),' cells']});   
         
    case 46 % 'Scatter: Ld vs InterDivision Time';
       pplot=' Ld vs CCT';  disp(pplot)
         %model='Gauss_normalized';
        [fit_out, fit_results, txt_x, txt_y]=plot_corr(TT, Ld, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Interdivision time, min','Ld, \mum',{'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
         
    case 64 % 'Scatter: InterDivision Time vs Ld';
       pplot='CCT vs Ld';  disp(pplot)
         %model='Gauss_normalized';
      [fit_out, fit_results, txt_x, txt_y]=plot_corr(Ld, TT, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Ld, \mum', 'Interdivision time, min', {'Single cell Data', ['binned, 1 point=',num2str(n_bins2),' cells']});   
   
    case 56 % 'Scatter: DR vs InterDivision Time';
       pplot=' DR vs CCT';  disp(pplot)
         %model='Gauss_normalized';
        [fit_out, fit_results, txt_x, txt_y]=plot_corr(TT, DR, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Interdivision time, min', 'DR', {'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
         
    case 65 % 'Scatter: InterDivision Time vs DR';
       pplot='CCT vs DR';  disp(pplot)
         %model='Gauss_normalized';
      [fit_out, fit_results, txt_x, txt_y]=plot_corr(DR, TT, n_bins2, 1); model='no fit';  
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('DR', 'Interdivision time, min', {'Single cell Data', ['binned, 1 point=',num2str(n_bins2),' cells']});  
         
     % Correlation between generations 
    case 111 % 'Scatter: Growth rate Daughter vs Mother';
       pplot='Growth rate: Daughter vs Mother';  disp(pplot)
         %model='Gauss_normalized';
       [fit_out, fit_results, txt_x, txt_y]=plot_corr(Alpha_a{1}, Alpha_d{1}, n_bins2, 1); 
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('\alpha_i, min^{-1}', '\alpha_{i+1}, min^{-1}', {'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
    
    case 211 % 'Scatter: Growth rate Sister vs Sister';
       pplot='Growth rate: Sister vs Sister';  disp(pplot)
         %model='Gauss_normalized';
       [fit_out, fit_results, txt_x, txt_y]=plot_corr(Alpha_s1{1}, Alpha_s2{1}, n_bins2, 1); 
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('\alpha_{s1}, min^{-1}', '\alpha_{s2}, min^{-1}', {'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});

    case 122 % 'Scatter: DeltaL Daughter vs Mother';
       pplot='DeltaL: Daughter vs Mother';  disp(pplot)
         %model='Gauss_normalized';
       [fit_out, fit_results, txt_x, txt_y]=plot_corr(DeltaL_a{1}, DeltaL_d{1}, n_bins2, 1); 
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('\DeltaL_i, \mum', '\DeltaL_{i+1}, \mum', {'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
    
    case 222 % 'Scatter: DeltaL Sister vs Sister';
       pplot='DeltaL: Sister vs Sister';  disp(pplot)
         %model='Gauss_normalized';
       [fit_out, fit_results, txt_x, txt_y]=plot_corr(DeltaL_s1{1}, DeltaL_s2{1}, n_bins2, 1); 
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('\DeltaL_{s1}, \mum', '\DeltaL_{s2}, \mum', {'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});     
    
    case 133 % 'Scatter: Lb Daughter vs Mother';
       pplot='Lb: Daughter vs Mother';  disp(pplot)
         %model='Gauss_normalized';
       [fit_out, fit_results, txt_x, txt_y]=plot_corr(Lb_a{1}, Lb_d{1}, n_bins2, 1); 
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Lb_i, \mum', 'Lb_{i+1}, \mum', {'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
    
    case 233 % 'Scatter: Lb Sister vs Sister';
       pplot='Lb: Sister vs Sister';  disp(pplot)
         %model='Gauss_normalized';
       [fit_out, fit_results, txt_x, txt_y]=plot_corr(Lb_s1{1}, Lb_s2{1}, n_bins2, 1); 
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Lb_{s1}, \mum', 'Lb_{s2}, \mum', {'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});  
         
    case 144 % 'Scatter: Ld Daughter vs Mother';
       pplot='Ld: Daughter vs Mother';  disp(pplot)
         %model='Gauss_normalized';
       [fit_out, fit_results, txt_x, txt_y]=plot_corr(Ld_a{1}, Ld_d{1}, n_bins2, 1); 
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Ld_i, \mum', 'Ld_{i+1}, \mum', {'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
    
    case 244 % 'Scatter: Ld Sister vs Sister';
       pplot='Ld: Sister vs Sister';  disp(pplot)
         %model='Gauss_normalized';
       [fit_out, fit_results, txt_x, txt_y]=plot_corr(Ld_s1{1}, Ld_s2{1}, n_bins2, 1); 
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('Ld_{s1}, \mum', 'Ld_{s2}, \mum', {'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});    
         
     case 155 % 'Scatter: DR Daughter vs Mother';
       pplot='DR: Daughter vs Mother';  disp(pplot)
         %model='Gauss_normalized';
       [fit_out, fit_results, txt_x, txt_y]=plot_corr(DR_a{1}, DR_d{1}, n_bins2, 1); 
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('DR_i', 'DR_{i+1}', {'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
    
    case 255 % 'Scatter: DR Sister vs Sister';
       pplot='DR: Sister vs Sister';  disp(pplot)
         %model='Gauss_normalized';
       [fit_out, fit_results, txt_x, txt_y]=plot_corr(DR_s1{1}, DR_s2{1}, n_bins2, 1); 
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('DR_{s1}', 'DR_{s2}', {'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});           
         
     case 166 % 'Scatter: CCT Daughter vs Mother';
       pplot='CCT: Daughter vs Mother';  disp(pplot)
         %model='Gauss_normalized';
       [fit_out, fit_results, txt_x, txt_y]=plot_corr(TT_a{1}, TT_d{1}, n_bins2, 1); 
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('T_i, min', 'T_{i+1}, min', {'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
    
    case 266 % 'Scatter: CCT Sister vs Sister';
       pplot='CCT: Sister vs Sister';  disp(pplot)
         %model='Gauss_normalized';
       [fit_out, fit_results, txt_x, txt_y]=plot_corr(TT_s1{1}, TT_s2{1}, n_bins2, 1); 
         %change_appearance('\DeltaL, \mum','\alpha, min^{-1})',{'Exp data',[model, ' fit']},num2str(fit_results,3));
         change_appearance('T_{s1}, min', 'T_{s2}, min', {'Single cell Data',['binned, 1 point=',num2str(n_bins2),' cells']});
     otherwise
       disp('I apologize for my limited abilities, but such plot is not included')
       continiue
         
 end
 
 % collect fit results
 Out(ii).plot=pplot;
 Out(ii).fit=fit_out;
 Out(ii).model=model;
 Out(ii).fit_results=fit_results;
 Out(ii).corr=corr_coeff;
 
    
end

% %get_plotting_options();
% plot_hist(Alpha,n_bins);
%  change_appearance('\alpha','Frequency','Gauss_normalized')
% plot_hist(DeltaL,n_bins);
%  change_appearance('\DeltaL','Frequency','LogNormal')
% plot_corr(DeltaL,Alpha,n_bins2,1);
%  change_appearance('\DeltaL','alpha')

 
 function  [fit_out, fit_results, txt_x, txt_y]=plot_hist(XX,n_bins,model)
    fit_out=[];  fit_results=[];
    [YY0,XX0]=hist(XX,n_bins);
    YY1=YY0'/sum(YY0)/(XX0(2)-XX0(1)); XX1=XX0';
    figure; hold on
     %plot(XX1,YY1,'Marker',mrk2,'MarkerFaceColor',col_m,'MarkerEdgeColor',col_m,'MarkerSize',ms2)
     plot(XX1,YY1,mrk2,'MarkerFaceColor',col_m,'MarkerEdgeColor',col_m,'MarkerSize',ms2)
   if nargin==3
     disp(['fiiting by ',model])  
     [coeff,fit_out,YYfit]=oneFit(XX1,YY1,model);
      ci = confint(fit_out); 
      fit_results=[coeff;ci];
       fit_results=(fit_results)';
     plot(XX1,YYfit,mrk_fit,'Color',col_f,'LineWidth',lw2)
     disp(fit_out);
     txt_x=min(XX1)+0.6*(max(XX1)-min(XX1)); txt_y=min(YY1)+0.3*(max(YY1)-min(YY1)); 
   end   
 end
 
  function [fit_out, fit_results, txt_x, txt_y]=plot_corr(XX,YY,n_bins2,choice,model)
    fit_out=[];  fit_results=[];  
    figure; hold on
    %plot(XX,YY,'Marker',mrk2,'MarkerFaceColor',col_m,'MarkerEdgeColor',col_m,'MarkerSize',ms3)
    plot(XX,YY,mrk2,'MarkerFaceColor',col_m,'MarkerEdgeColor',col_m,'MarkerSize',ms3)
    [XX1,YY1]=bin2_fixN(XX,YY,n_bins2,choice); 
    %h = denScatter2(X,Y,radius)
    %plot(XX1,YY1,'Marker',mrk2,'MarkerFaceColor',col_m,'MarkerEdgeColor',col_m,'MarkerSize',ms2)
    %plot(XX1,YY1,[mrk2,'-'],'MarkerFaceColor','w','Color',[1,1,1]-0.4*([1,1,1]-col_m),'MarkerEdgeColor',[1,1,1]-0.4*([1,1,1]-col_m),'MarkerSize',ms2,'LineWidth',lw2)
     plot(XX1,YY1,[mrk2,],'MarkerFaceColor','w','Color',0.5*col_m,'MarkerEdgeColor',0.5*col_m,'MarkerSize',ms2,'LineWidth',lw1)
    
    corr_coeff=corr(XX,YY,'type','Kendall');
     txt_x=min(XX)+0.8*(max(XX)-min(XX)); txt_y=min(YY)+0.3*(max(YY)-min(YY)); 
     text(txt_x,txt_y,['corr = ',num2str(corr_coeff,2)],'FontSize',fs3);
   if nargin==5
     disp(['fiiting by ',model])  
     [coeff,fit_out,YYfit]=oneFit(XX1,YY1,model);
       ci = confint(fit_out); 
      fit_results=[coeff;ci];
       fit_results=(fit_results)';
     plot(XX1,YYfit,mrk_fit,'Color',col_f,'LineWidth',lw2)
     disp(fit_out);
     txt_x=min(XX1)+0.6*(max(XX1)-min(XX1)); txt_y=min(YY1)+0.3*(max(YY1)-min(YY1));
   end   
  end
% 
%   function get_plotting_options()
%   % Appearance
%    fs1=18; fs2=16; fs3=14; % font sizes for Labels, tick marks, text, respectively
%    lw1=3; lw2=2; lw3=1;       % line width
%     ms1=12; ms2=8; ms3=3;
%    mrk1='x';mrk2='o';mrk3='.';
%     col_m='b';
%    mrk_fit='-';
%    col_f='k';
%   end
  
 function change_appearance(x_lab, y_lab, lgnd, fit_results, title_txt)
   set (gca, 'FontSize', fs2);
   xlabel(x_lab, 'FontSize',fs1); ylabel(y_lab, 'FontSize',fs1);
   if nargin>2
      legend(lgnd,'FontSize',fs3);
      if nargin>3
         text(txt_x, txt_y, fit_results, 'FontSize',fs3); 
         if nargin>4
           title(title_txt,'FontSize',fs3);
         end
     end  
   end
 end
      

end
  
     
 