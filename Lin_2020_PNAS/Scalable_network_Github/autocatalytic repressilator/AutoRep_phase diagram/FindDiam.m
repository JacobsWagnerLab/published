
% Finding diameter of limit cycle
% diameter is defined by sup { d(x1,x2) | x1, x2 in A} 
%                         A

% Input: the multi-dimensional trajectory that is already converged to the 
% omega-limit set.

function [Diam] = FindDiam(TrajLC)

L = size(TrajLC,1);
DistMat = zeros(L,L);

for j1 = 1:L
    
    for j2 = 1:j1
        
        xtemp1 = TrajLC(j1,:);
        xtemp2 = TrajLC(j2,:);
        
        DistMat(j1,j2) = norm(xtemp1 - xtemp2);   
                
    end
    
end

Diam = max(max(DistMat));


