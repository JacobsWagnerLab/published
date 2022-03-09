
function [vec] = PM_1D_gaussian(L, sig)

%L = 10;
%sig = 10;

vec = zeros(L,1);

for j = 1:L    
        
    origin = (L+1)/2;
    
    x = j - origin;        
        
	temp = (x^2) / sig;        
        
    vec(j) = exp( (-1/2)*temp );
        
end

vec = vec / sum(sum(vec));

end
