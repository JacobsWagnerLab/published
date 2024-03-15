clc
clear all;

load('E:\Shares\Data_02\Yashna Thappeta\Figures_Data_Analysis\Silvia_coculture_experiments\230520_analysis\fig6b_data.mat') %loading data structure with measures areas per strain/mix/replica

%% Imports data from loaded structure to variables

AreaTransMix1R1_YFP=res_area.t1_mix1_wt(~isnan(res_area.t1_mix1_wt));
AreaTransMix1R1_CFP=res_area.t1_mix1_mut(~isnan(res_area.t1_mix1_mut));
AreaTransMix1R3_YFP=res_area.t3_mix1_wt(~isnan(res_area.t3_mix1_wt));
AreaTransMix1R3_CFP=res_area.t3_mix1_mut(~isnan(res_area.t3_mix1_mut));

AreaTransMix2R1_CFP=res_area.t1_mix2_wt(~isnan(res_area.t1_mix2_wt));
AreaTransMix2R1_YFP=res_area.t1_mix2_mut(~isnan(res_area.t1_mix2_mut));
AreaTransMix2R2_CFP=res_area.t2_mix2_wt(~isnan(res_area.t2_mix2_wt));
AreaTransMix2R2_YFP=res_area.t2_mix2_mut(~isnan(res_area.t2_mix2_mut));
AreaTransMix2R3_CFP=res_area.t3_mix2_wt(~isnan(res_area.t3_mix2_wt));
AreaTransMix2R3_YFP=res_area.t3_mix2_mut(~isnan(res_area.t3_mix2_mut));

AreaExpoMix1R1_YFP=res_area_exp.t1_mix1_wt(~isnan(res_area_exp.t1_mix1_wt));
AreaExpoMix1R1_CFP=res_area_exp.t1_mix1_mut(~isnan(res_area_exp.t1_mix1_mut));
AreaExpoMix1R3_YFP=res_area_exp.t3_mix1_wt(~isnan(res_area_exp.t3_mix1_wt));
AreaExpoMix1R3_CFP=res_area_exp.t3_mix1_mut(~isnan(res_area_exp.t3_mix1_mut));

AreaExpoMix2R1_CFP=res_area_exp.t1_mix2_wt(~isnan(res_area_exp.t1_mix2_wt));
AreaExpoMix2R1_YFP=res_area_exp.t1_mix2_mut(~isnan(res_area_exp.t1_mix2_mut));
AreaExpoMix2R2_CFP=res_area_exp.t2_mix2_wt(~isnan(res_area_exp.t2_mix2_wt));
AreaExpoMix2R2_YFP=res_area_exp.t2_mix2_mut(~isnan(res_area_exp.t2_mix2_mut));
AreaExpoMix2R3_CFP=res_area_exp.t3_mix2_wt(~isnan(res_area_exp.t3_mix2_wt));
AreaExpoMix2R3_YFP=res_area_exp.t3_mix2_mut(~isnan(res_area_exp.t3_mix2_mut));

%% Replicas merged

AreaExpoMix1_WT=cat(1,AreaExpoMix1R1_YFP,AreaExpoMix1R3_YFP);
AreaExpoMix1_glg=cat(1,AreaExpoMix1R1_CFP,AreaExpoMix1R3_CFP);
AreaTransMix1_WT=cat(1,AreaTransMix1R1_YFP,AreaTransMix1R3_YFP);
AreaTransMix1_glg=cat(1,AreaTransMix1R1_CFP,AreaTransMix1R3_CFP);

AreaExpoMix2_WT=cat(1,AreaExpoMix2R1_CFP,AreaExpoMix2R2_CFP,AreaExpoMix2R3_CFP);
AreaExpoMix2_glg=cat(1,AreaExpoMix2R1_YFP,AreaExpoMix2R2_YFP,AreaExpoMix2R3_YFP);
AreaTransMix2_WT=cat(1,AreaTransMix2R1_CFP,AreaTransMix2R2_CFP,AreaTransMix2R3_CFP);
AreaTransMix2_glg=cat(1,AreaTransMix2R1_YFP,AreaTransMix2R2_YFP,AreaTransMix2R3_YFP);


%%
%Mix 1 WT is YFP and glg is CFP

meanWT_Mix1R1_Expo=bootstrp(500,@mean,AreaExpoMix1R1_YFP);
meanglg_Mix1R1_Expo=bootstrp(500,@mean,AreaExpoMix1R1_CFP);
diff_Mix1R1_Expo=zeros(500,1);

meanWT_Mix1R1_Trans=bootstrp(500,@mean,AreaTransMix1R1_YFP);
meanglg_Mix1R1_Trans=bootstrp(500,@mean,AreaTransMix1R1_CFP);
diff_Mix1R1_Trans=zeros(500,1);

meanWT_Mix1R3_Expo=bootstrp(500,@mean,AreaExpoMix1R3_YFP);
meanglg_Mix1R3_Expo=bootstrp(500,@mean,AreaExpoMix1R3_CFP);
diff_Mix1R3_Expo=zeros(500,1);
meanWT_Mix1R3_Trans=bootstrp(500,@mean,AreaTransMix1R3_YFP);
meanglg_Mix1R3_Trans=bootstrp(500,@mean,AreaTransMix1R3_CFP);
diff_Mix1R3_Trans=zeros(500,1);

%Mix 2 WT is CFP and glg is YFP

meanWT_Mix2R1_Expo=bootstrp(500,@mean,AreaExpoMix2R1_CFP);
meanglg_Mix2R1_Expo=bootstrp(500,@mean,AreaExpoMix2R1_YFP);
diff_Mix2R1_Expo=zeros(500,1);

meanWT_Mix2R1_Trans=bootstrp(500,@mean,AreaTransMix2R1_CFP);
meanglg_Mix2R1_Trans=bootstrp(500,@mean,AreaTransMix2R1_YFP);
diff_Mix2R1_Trans=zeros(500,1);

meanWT_Mix2R2_Expo=bootstrp(500,@mean,AreaExpoMix2R2_CFP);
meanglg_Mix2R2_Expo=bootstrp(500,@mean,AreaExpoMix2R2_YFP);
diff_Mix2R2_Expo=zeros(500,1);
meanWT_Mix2R2_Trans=bootstrp(500,@mean,AreaTransMix2R2_CFP);
meanglg_Mix2R2_Trans=bootstrp(500,@mean,AreaTransMix2R2_YFP);
diff_Mix2R2_Trans=zeros(500,1);

meanWT_Mix2R3_Expo=bootstrp(500,@mean,AreaExpoMix2R3_CFP);
meanglg_Mix2R3_Expo=bootstrp(500,@mean,AreaExpoMix2R3_YFP);
diff_Mix2R3_Expo=zeros(500,1);
meanWT_Mix2R3_Trans=bootstrp(500,@mean,AreaTransMix2R3_CFP);
meanglg_Mix2R3_Trans=bootstrp(500,@mean,AreaTransMix2R3_YFP);
diff_Mix2R3_Trans=zeros(500,1);

for i=1:500
    
    diff_Mix1R1_Expo(i,1)=100*((meanWT_Mix1R1_Expo(i,1)-meanglg_Mix1R1_Expo(i,1))/meanWT_Mix1R1_Expo(i,1));
    diff_Mix1R1_Trans(i,1)=100*((meanWT_Mix1R1_Trans(i,1)-meanglg_Mix1R1_Trans(i,1))/meanWT_Mix1R1_Trans(i,1));
    diff_Mix1R3_Expo(i,1)=100*((meanWT_Mix1R3_Expo(i,1)-meanglg_Mix1R3_Expo(i,1))/meanWT_Mix1R3_Expo(i,1));
    diff_Mix1R3_Trans(i,1)=100*((meanWT_Mix1R3_Trans(i,1)-meanglg_Mix1R3_Trans(i,1))/meanWT_Mix1R3_Trans(i,1));
    
    diff_Mix2R1_Expo(i,1)=100*((meanWT_Mix2R1_Expo(i,1)-meanglg_Mix2R1_Expo(i,1))/meanWT_Mix2R1_Expo(i,1));
    diff_Mix2R1_Trans(i,1)=100*((meanWT_Mix2R1_Trans(i,1)-meanglg_Mix2R1_Trans(i,1))/meanWT_Mix2R1_Trans(i,1));
    diff_Mix2R2_Expo(i,1)=100*((meanWT_Mix2R2_Expo(i,1)-meanglg_Mix2R2_Expo(i,1))/meanWT_Mix2R2_Expo(i,1));
    diff_Mix2R2_Trans(i,1)=100*((meanWT_Mix2R2_Trans(i,1)-meanglg_Mix2R2_Trans(i,1))/meanWT_Mix2R2_Trans(i,1));
    diff_Mix2R3_Expo(i,1)=100*((meanWT_Mix2R3_Expo(i,1)-meanglg_Mix2R3_Expo(i,1))/meanWT_Mix2R3_Expo(i,1));
    diff_Mix2R3_Trans(i,1)=100*((meanWT_Mix2R3_Trans(i,1)-meanglg_Mix2R3_Trans(i,1))/meanWT_Mix2R3_Trans(i,1));
    
end

mean_DiffMix1R1_Expo= mean(diff_Mix1R1_Expo);
mean_DiffMix1R1_Trans= mean(diff_Mix1R1_Trans);
mean_DiffMix1R3_Expo= mean(diff_Mix1R3_Expo);
mean_DiffMix1R3_Trans= mean(diff_Mix1R3_Trans);

mean_DiffMix2R1_Expo= mean(diff_Mix2R1_Expo);
mean_DiffMix2R1_Trans= mean(diff_Mix2R1_Trans);
mean_DiffMix2R2_Expo= mean(diff_Mix2R2_Expo);
mean_DiffMix2R2_Trans= mean(diff_Mix2R2_Trans);
mean_DiffMix2R3_Expo= mean(diff_Mix2R3_Expo);
mean_DiffMix2R3_Trans= mean(diff_Mix2R3_Trans);

std_DiffMix1R1_Expo= std(diff_Mix1R1_Expo);
std_DiffMix1R1_Trans= std(diff_Mix1R1_Trans);
std_DiffMix1R3_Expo= std(diff_Mix1R3_Expo);
std_DiffMix1R3_Trans= std(diff_Mix1R3_Trans);

std_DiffMix2R1_Expo= std(diff_Mix2R1_Expo);
std_DiffMix2R1_Trans= std(diff_Mix2R1_Trans);
std_DiffMix2R2_Expo= std(diff_Mix2R2_Expo);
std_DiffMix2R2_Trans= std(diff_Mix2R2_Trans);
std_DiffMix2R3_Expo= std(diff_Mix2R3_Expo);
std_DiffMix2R3_Trans= std(diff_Mix2R3_Trans);

%%

meanWT_Mix1_Expo=bootstrp(1000,@mean,AreaExpoMix1_WT);
meanglg_Mix1_Expo=bootstrp(1000,@mean,AreaExpoMix1_glg);
diff_Mix1_Expo=zeros(1000,1);
meanWT_Mix1_Trans=bootstrp(1000,@mean,AreaTransMix1_WT);
meanglg_Mix1_Trans=bootstrp(1000,@mean,AreaTransMix1_glg);
diff_Mix1_Trans=zeros(1000,1);

meanWT_Mix2_Expo=bootstrp(1000,@mean,AreaExpoMix2_WT);
meanglg_Mix2_Expo=bootstrp(1000,@mean,AreaExpoMix2_glg);
diff_Mix2_Expo=zeros(1000,1);
meanWT_Mix2_Trans=bootstrp(1000,@mean,AreaTransMix2_WT);
meanglg_Mix2_Trans=bootstrp(1000,@mean,AreaTransMix2_glg);
diff_Mix2_Trans=zeros(1000,1);


for i=1:1000
    
    diff_Mix1_Expo(i,1)=100*((meanWT_Mix1_Expo(i,1)-meanglg_Mix1_Expo(i,1))/meanWT_Mix1_Expo(i,1));
    diff_Mix1_Trans(i,1)=100*((meanWT_Mix1_Trans(i,1)-meanglg_Mix1_Trans(i,1))/meanWT_Mix1_Trans(i,1));
    
    diff_Mix2_Expo(i,1)=100*((meanWT_Mix2_Expo(i,1)-meanglg_Mix2_Expo(i,1))/meanWT_Mix2_Expo(i,1));
    diff_Mix2_Trans(i,1)=100*((meanWT_Mix2_Trans(i,1)-meanglg_Mix2_Trans(i,1))/meanWT_Mix2_Trans(i,1));
      
end

mean_DiffMix1_Expo= mean(diff_Mix1_Expo);
mean_DiffMix1_Trans= mean(diff_Mix1_Trans);
std_DiffMix1_Expo= std(diff_Mix1_Expo);
std_DiffMix1_Trans= std(diff_Mix1_Trans);

mean_DiffMix2_Expo= mean(diff_Mix2_Expo);
mean_DiffMix2_Trans= mean(diff_Mix2_Trans);
std_DiffMix2_Expo= std(diff_Mix2_Expo);
std_DiffMix2_Trans= std(diff_Mix2_Trans);

%% By replica
Stats=zeros(6,2);

Stats(1,1)=mean_DiffMix1R3_Expo;
Stats(2,1)=mean_DiffMix1R3_Trans;
Stats(3,1)=mean_DiffMix2R2_Expo;
Stats(4,1)=mean_DiffMix2R2_Trans;
Stats(5,1)=mean_DiffMix2R1_Expo;
Stats(6,1)=mean_DiffMix2R1_Trans;

Stats(1,2)=std_DiffMix1R3_Expo;
Stats(2,2)=std_DiffMix1R3_Trans;
Stats(3,2)=std_DiffMix2R2_Expo;
Stats(4,2)=std_DiffMix2R2_Trans;
Stats(5,2)=std_DiffMix2R1_Expo;
Stats(6,2)=std_DiffMix2R1_Trans;

%% Merged replicas

StatsAll=zeros(4,2);

StatsAll(1,1)=mean_DiffMix1_Expo;
StatsAll(2,1)=mean_DiffMix1_Trans;
StatsAll(3,1)=mean_DiffMix2_Expo;
StatsAll(4,1)=mean_DiffMix2_Trans;

StatsAll(1,2)=std_DiffMix1_Expo;
StatsAll(2,2)=std_DiffMix1_Trans;
StatsAll(3,2)=std_DiffMix2_Expo;
StatsAll(4,2)=std_DiffMix2_Trans;

%% R1 complete
StatsR1=zeros(4,2);

StatsR1(1,1)=mean_DiffMix1R1_Expo;
StatsR1(2,1)=mean_DiffMix1R1_Trans;
StatsR1(3,1)=mean_DiffMix2R3_Expo;
StatsR1(4,1)=mean_DiffMix2R3_Trans;

StatsR1(1,2)=std_DiffMix1R1_Expo;
StatsR1(2,2)=std_DiffMix1R1_Trans;
StatsR1(3,2)=std_DiffMix2R3_Expo;
StatsR1(4,2)=std_DiffMix2R3_Trans;

%%
% overlay fitted line
figure('Name','1','WindowState','maximized');
hold on
q3=tiledlayout(2,2); % Requires R2019b or later
 
hold on
x = 1:6;

scatter(x,Stats(:,1))                
xlim([0 7]);
ylim([0 20]);
hold on

er = errorbar(x,Stats(:,1),Stats(:,2));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off
% set(gca,'YScale','log'); %// NEW
% set(gca,'YDir','reverse')
set(gca,'FontSize',12)
set(gca,'FontName','Arial')
ylabel('Size Difference (%)','FontSize',12)
xlim([0 7]);
ylim([0 20]);
hold off
 
%exportgraphics(q3,'CJW7323_GR.png','Resolution',300)

%% Merged replicas

% overlay fitted line
figure('Name','1','WindowState','maximized');
hold on
q3=tiledlayout(2,2); % Requires R2019b or later
 
hold on
x = 1:4;

scatter(x,StatsAll(:,1))                
xlim([0 11]);
ylim([0 20]);
hold on

er = errorbar(x,StatsAll(:,1),StatsAll(:,2));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off
% set(gca,'YScale','log'); %// NEW
% set(gca,'YDir','reverse')
set(gca,'FontSize',12)
set(gca,'FontName','Arial')
ylabel('Size Difference (%)','FontSize',12)
xlim([0 5]);
ylim([0 20]);
hold off
 
%exportgraphics(q3,'CJW7323_GR.png','Resolution',300)

%% R1 complete

% overlay fitted line
figure('Name','3','WindowState','maximized');
hold on
q4=tiledlayout(2,2); % Requires R2019b or later
 
hold on
x = 1:4;

scatter(x,StatsR1(:,1))                
xlim([0 5]);
ylim([0 20]);
hold on

er = errorbar(x,StatsR1(:,1),StatsR1(:,2));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off
% set(gca,'YScale','log'); %// NEW
% set(gca,'YDir','reverse')
set(gca,'FontSize',12)
set(gca,'FontName','Arial')
ylabel('Size Difference (%)','FontSize',12)
xlim([0 5]);
ylim([0 20]);
hold off
 
%exportgraphics(q3,'CJW7323_GR.png','Resolution',300)
%%
figure('Name','4','WindowState','maximized');
hold on
q5=tiledlayout(1,2); % Requires R2019b or later
nexttile
hold on
histogram(AreaTransMix1_WT,'Normalization','probability')
histogram(AreaTransMix1_glg,'Normalization','probability')
nexttile
hold on
histogram(AreaTransMix2_WT,'Normalization','probability')
histogram(AreaTransMix2_glg,'Normalization','probability')
hold off

