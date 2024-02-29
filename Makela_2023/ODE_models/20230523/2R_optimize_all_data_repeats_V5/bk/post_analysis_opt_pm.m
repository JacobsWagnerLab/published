
% Comparing parameter before and after optimization
% Manaully change location to the datasets

param_matrix = NaN(10,4);
dataset_index = 2;

% Rows
% [1] r1, [2] r2, [3] K1, [4] K2, [5] d, [6] Zs, 
% [7] c, [8] cp, [9] xini, [10] yini

% Cols:
% [1,4,5] condition 1 to 3, before optimization
% [2,4,6] condition 1 to 3, after optimization


tag = {'r1','r2','K1','K2','d','Zs','c','cp','xini','yini'};

for r = 1:10    
        
    param_matrix(r,1) = output0.param{dataset_index}.(tag{r});   
    param_matrix(r,2) = output1.param{dataset_index}.(tag{r});
    param_matrix(r,3) = output2.param{dataset_index}.(tag{r});
    param_matrix(r,4) = output3.param{dataset_index}.(tag{r});
        
end

xtag = {'ini', 'opt1', 'opt2', 'opt3'};
% 'Gly, ini', 'Gly, opt', 'Ala, ini', 'Ala, opt'};

figure;

for r = 1:10
    
    subplot(2,5,r);
    
    bar(param_matrix(r,:));
    xticklabels(xtag);
    ylabel(tag{r});

end










