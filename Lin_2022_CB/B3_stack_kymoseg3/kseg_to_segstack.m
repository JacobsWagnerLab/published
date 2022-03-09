
% Convert kymograph segmentation into 2D segmentation stacks.

function [segmatN] = kseg_to_segstack(segmat, kseg, param)

% conver space-time data into segmentation matrix stacks

segmatN = 0*segmat(:,:,param.Nrange(2));

for fn = param.Nrange(1):param.Nrange(2)
    
    maskA = (segmat(:,:,fn) > 0);
    maskB = 0*maskA;
    
    L = size(maskB);
    
    seg_mask = kron( kseg(:,fn), ones(1,L(2)) );    
    maskB(param.Yrange(1):param.Yrange(2),:) = seg_mask;
    
    % Use both the pre-segmentation mask (cell contour) (maskA)
    % and the kymoseg mask (maskB) to generate the final cell outline
    
    segmatN(:,:,fn) = maskA .* maskB; 
    
end
