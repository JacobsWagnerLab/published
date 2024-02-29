
function [param] = load_optimized_param(param_opt)

param = {};  
    
param.r1 = param_opt(1);    
param.r2 = param_opt(2); 
param.K1 = param_opt(3);   
param.K2 = param_opt(4);   
    
param.d = param_opt(5);    
param.Zs = param_opt(6); 
param.c = param_opt(7);   
param.cp = param_opt(8);   
    
param.xini = param_opt(9);  
param.yini = param_opt(10); 

param.k1 = param_opt(11);  
param.km1 = param_opt(12); 
param.k2 = param_opt(13);  
param.k3 = param_opt(14); 

end

