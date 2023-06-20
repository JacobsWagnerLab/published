%% segmentation
%
% -Purpose-
%   A rudimentary implementation for image segmentation. Looks simply at
%   the image intensity threshold and object size threholds. 
%
% -Input-
%   - im (numeric array): 2D matrix to represent the image
%   - threshold (double): intensity threshold, above which the pixel is
%   zeroed
%   - minsz (double): minimum number of pixels an object should have
%   - maxsz (double): maximum number of pixels an object should have
%
% -Output-
%   - mask (logical array): Image mask that identifies the objects
%
% -Author-
%   Yingjie Xiang, 2019-01-17
%
% -Patch Notes-
%   2019-01-17: created the function

function mask = segmentation(im,threshold,minsz,maxsz)

% Use very simple intensity threshold to segment cells
im(im>threshold) = 0;
bw = imbinarize(im);
cc = bwconncomp(bw);
keep = cellfun(@(x) length(x) > minsz & length(x) < maxsz,cc.PixelIdxList);
cc.PixelIdxList(~keep) = [];
cc.NumObjects = length(cc.PixelIdxList);

%label = labelmatrix(cc);
%imshow(label2rgb(label,'lines','w','shuffle'));

% Generate mask list from the connected components
mask = false(cc.ImageSize);
for ii = 1:cc.NumObjects
    mask_for_one = false(cc.ImageSize);
    mask_for_one(ind2sub(size(mask_for_one),cc.PixelIdxList{ii})) = true;
    mask = mask | mask_for_one;
end

end

