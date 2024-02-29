function [ROI, xstart, ystart] = getLocalizationROI(point_pos,im, windowSize)

[sizey, sizex] = size(im);
X0=point_pos(1);
Y0=point_pos(2);

%round X0, Y0 to use as matrix locations
X0_int = round(X0);
Y0_int = round(Y0);
windowSize = round(windowSize); %radius should already be an integer anyway

% setup the limits of the cropped image
xstart =  X0_int-windowSize;
xfinish = X0_int+windowSize;
ystart =  Y0_int-windowSize;
yfinish = Y0_int+windowSize;
% check if any of the limits are out of bounds - if so, skip that point
if (xstart<1) || (xstart > sizex) ||  (xfinish<1) || (xfinish > sizex) ...
        || (ystart<1) || (ystart > sizey) ||  (yfinish<1) || (yfinish > sizey)
    
    ROI = [];
    %warning('sub im outide limits');
else
    %crop to a small area around the point
    ROI = im( ystart:yfinish, xstart:xfinish);
end

end