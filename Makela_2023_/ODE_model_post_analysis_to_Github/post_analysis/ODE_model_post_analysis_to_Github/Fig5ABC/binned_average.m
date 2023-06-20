
function [output] = binned_average(x,y,w,xBD)

% Get the moving averaging 

%x = exp_data{1}.area;
%y = exp_data{1}.exp;
%w = 0.5; 
%xBD = [0 10];

edge = (xBD(1):w:xBD(2))';

bnum = length(edge);

dataB = zeros(bnum,3);  % count; sum, mean


for j = 1:length(x)
    
    ind = bnum;
    
    for b = 1:(bnum-1)
        
        if ( x(j) > edge(b) ) && ( x(j) < edge(b+1) )
            ind = b;
            break;
        end
        
    end
     
    dataB(ind,1) = dataB(ind,1) + 1;    
    dataB(ind,2) = dataB(ind,2) + y(j);
    
end

dataB(:,3) = dataB(:,2)./dataB(:,1);

xval = edge(1:end-1) + (w/2);
yval = dataB(1:end-1,3);

%=================================================================

output= {};
output.dataB = dataB;
output.xval = xval;
output.yval = yval;
output.edge = edge;

%

%figure;
%plot(area, Adot, '.'); hold on;
%plot(xval, yval, 'o-');
%xlim([0 10]);
