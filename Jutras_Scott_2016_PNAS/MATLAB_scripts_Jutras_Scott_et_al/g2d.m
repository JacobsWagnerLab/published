function [xy] = g2d(image)
% image is any 2D matrix, the center of mass is returned in a 2 element
% vector x,y
%
% Brad Parry, Christine Jacobs-Wagner lab, Yale University
image = double(image);
M0 = sum(image(:));
x = repmat(1:size(image,2),[size(image,1), 1]);
Mx = sum(sum(x.*image));
ay(:,1) = 1:size(image,1);
y = repmat(ay,[1,size(image,2)]);
My = sum(sum(y.*image));
xy = [Mx/M0, My/M0];
end