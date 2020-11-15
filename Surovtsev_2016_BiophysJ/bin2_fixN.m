function [Xout,Yout]=bin2_fixN(X,Y,N,choice)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function bin_fixN
%PURPOSE: bin [x,y] pairs according to x-values and return corresponding X-
%          and Y- mean values for each bin
%@author: Ivan Surovtsev 
%@date:  11.01.2012  
%@copyright 2012-2015 Yale University
%==========================================================================
%**********INPUT********:
% X, Y - column-vectors of the data, must be the same size 
% N    - number of data-points in each bin
%**********OUTPUT********:
%Xout,Yout - column vectors with mean values for each bin
%************************
% function does following:
% 1.sort data [X,Y] in increasing order of X
% 2.bins [X,Y] data to have N datapoints in each bin (tail of [X,Y], i.e last N_last<=N points falls in the last bin)
% 3.returns mean X and mean Y for each bin
%
%=========================================================================
% Publication: specify Journal name and year such as eLife 2014
% Publication Author(s):
% Publication Title:
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

     XY=[X,Y];
      XY=sortrows(XY,1);  
      
       Xout=[]; Yout=[];  
    
  if nargin<4, choice=1; end
       
     switch choice  
       case 1
         for kk=1:N:length(XY)-N+1
           Xout=[Xout;mean(XY(kk:kk+N-1,1))]; Yout=[Yout;mean(XY(kk:kk+N-1,2))];  
         end
         if kk+N-1<length(XY)
           Xout=[Xout;mean(XY(kk+N:end,1))]; Yout=[Yout;mean(XY(kk+N:end,2))];
         end
           
       case 2
         for kk=1:N:length(XY)-N+1
           Xout=[Xout;median(XY(kk:kk+N-1,1))]; Yout=[Yout;median(XY(kk:kk+N-1,2))];  
         end
         if kk+N-1<length(XY)
           Xout=[Xout;median(XY(kk+N:end,1))]; Yout=[Yout;median(XY(kk+N:end,2))];      
         end
     end
        
end