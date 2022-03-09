
function [A] = PM_2D_gaussian(L, sig)

%L = 10;
%sig = 10;

A = zeros(L,L);

for j1 = 1:L    
    for j2 = 1:L
        
        origin = (L+1)/2;
        
        x = j1 - origin;
        y = j2 - origin;
        
        r_sq = x^2 + y^2;
        
        A(j1, j2) = exp((-1)*(r_sq)/(2*sig));
        
    end
        
end

A = A / sum(sum(A));

%figure; imagesc(A);

