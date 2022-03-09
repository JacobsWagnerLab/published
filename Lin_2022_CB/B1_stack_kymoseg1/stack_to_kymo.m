
% This subscript create 2D kymograph from a segmented image stacks

function [matN] = stack_to_kymo(segmat, param)

Yrange = param.Yrange;   % Range of frames in image stacks
Nrange = param.Nrange;   % Range of y-axis (rows) in image 

%

data = {};

for j = Nrange(1) : Nrange(2)
    
    data{j}.seg = ( segmat(Yrange(1):Yrange(2),:,j) > 0);
    data{j}.col = Columnize(data{j}.seg);
    
end

% (B) Making space-time matrix. Each image become one column

matN = NaN( length(data{Nrange(1)}.col), Nrange(2)-Nrange(1)+1);

for j = Nrange(1) : Nrange(2)
    
    matN(:,j - Nrange(1) + 1 ) = data{j}.col; 
    
end

end

%imagesc(matN)

% ===================================================================== %

function [col] = Columnize(img)

% Collapse the x-direction of each microfluidic chamber;
% hence each image becomes one data column

col = sum(img,2);

end

% ===================================================================== %
