%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function [coeff,fit1,YYfit]=fitErf(XX,YY,model,weights)
%@author: Ivan Surovtsev (2013), adapted by Molly Scott
%@date: March 4, 2016
%==========================================================================
%************************Output**********************:
%coeff:                 are best-fit parameters (values) of the model
%fit1:                  is Matlab's fit object (useful for cofidence intevals, residuals etc)
%YYfit:                 returns values of the model with best-fit parameters values at points from XX
%                       (useful for plotting)
%************************Input**********************:
%XX:                    x-data to fit
%YY:                    y-data to fit (size(YY) must be equal size(XX))
%model:                 string specifing model from preset list.
%weights (optional):    array of wieghts for fitting (must the size of XX)
%==========================================================================
% This function performs a fitting of YY(XX) data with the model '2-erf. Allows
% the user to visualize the fit to the data and outputs a matrix. Should be
% used with the function calculate_width.m to determine the widths of fits
% to peaks in fluorescence.
%-------------------------------------------------------------------------- 
%-------------------------------------------------------------------------- 
function [coeff,fit1,YYfit]=fitErf(XX,YY,model,weights)

if nargin < 4, weights=ones(size(XX)); end

ffit=0;

switch model
          case '2_erf'
       model1 = fittype('a*(erf(x-x0)-erf(x-x0-w))/sqrt(2)+b','dependent',{'y'},'independent',{'x'},'coefficients',{'a', 'w','b','x0'});
         %model1 = fittype('a*(erf(x/sigma/sqrt(2))-erf((x-w)/sigma/sqrt(2)))/sqrt(2)','dependent',{'y'},'independent',{'x'},'coefficients',{'a', 'w', 'sigma'});
         x0=(max(XX)-min(XX))/2;
         par0 = [1, 0, 0,x0];
         ffit=2;
         lb=[0,0,0,0]; ub=[Inf,20,Inf,Inf]; 
         
          otherwise
        disp('such model is not included');
        return;
end

switch ffit
    case 2
      fit1 = fit(XX,YY,model1,'Startpoint',par0,'Weights',weights,'Lower',lb,'Upper',ub);
    case 1
      fit1 = fit(XX,YY,model1,'Startpoint',par0,'Weights',weights);
    case 0
      fit1 = fit(XX,YY,model1,'Startpoint',par0);  
end

switch model
        case '2_erf'    
        coeff=coeffvalues(fit1);
        a=coeff(1);
        w=coeff(2);
        b=coeff(3);
        x0=coeff(4);
        YYfit=a*(erf(XX-x0)-erf(XX-x0-w))/sqrt(2)+b; 
        
    otherwise
        disp('such model is not included');
        return;
end
    
   