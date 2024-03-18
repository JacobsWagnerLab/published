function [fractionalLength, halfLength, fractionalWidth, halfWidth, LWratio] = Cell_Projection(cellMesh, xys)

%{
-About-
This function uses Cell_Projection.m to extract pixels within dC of the 
cell centerline and greater than poleLen from the pole.

-Inputs-
cellMesh: the cellMesh as returned from Oufti

xys: n by 2 matrix of x coordinates (column 1) and y coordinates (column 2)

-Outputs-
fractionalLength: position along cell length from -1,1

halfLength: half the cell length

fractionalWidth: position along cell width from -1,1

halfWidth: half the cell width

LWratio: ratio of the cell length to width


-Keywords-
Pixel correlation
Ribosome segregation

-Dependencies-
Taylor_Smooth.m

-Author-
Bradley Parry October 10, 2018
%}

cellMesh = double(cellMesh);
% use Taylor_Smooth to remove irregularities from the Oufti cell mesh
centerLine = [(cellMesh(:,1)+cellMesh(:,3))/2 ((cellMesh(:,2)+cellMesh(:,4))/2)];
CENTERLINE = Taylor_Smooth(centerLine);
cellWidth = max(((cellMesh(:,1) - cellMesh(:,3)).^2 + (cellMesh(:,2) - cellMesh(:,4)).^2).^(1/2));
halfWidth = cellWidth/2;
cellLength = sum(sum((CENTERLINE(2:end,1:2) - CENTERLINE(1:end-1,1:2)).^2,2).^(1/2));    
halfLength = cellLength/2;
% useful for reconstructing 2d cell
LWratio = cellLength / cellWidth;

sz1 = size(xys,1);
fractionalLength = zeros(sz1, 1);
fractionalWidth = zeros(sz1, 1);

for q = 1:sz1
    xy = xys(q,:);

    % find the segments that are closest to point xy
    rxy = repmat(xy, [size(CENTERLINE,1),1]);
    [~,ix] = sort(sum((rxy - CENTERLINE).^2,2),'ascend');
    ix(3:end)=[];

    % is the point on the "left" or "right" side [used at the end] L = -1; R = 1
    LR = -1;
    L = min(sum((rxy - cellMesh(:,1:2)).^2, 2).^(1/2));
    R = min(sum((rxy - cellMesh(:,3:4)).^2, 2).^(1/2));
    LR(R>L) = 1;

    % now shift rxy and cellSeg by shiftFactor to center the data on the origin
    cellSeg = CENTERLINE(ix(1:2),:);
    shiftFactor = cellSeg(1,:);
    cellSeg(:,1) = cellSeg(:,1) - shiftFactor(1);
    cellSeg(:,2) = cellSeg(:,2) - shiftFactor(2);
    xy = xy - shiftFactor;

    % get the projection along cell length
    A = [cellSeg(end,1), cellSeg(end,2)]';
    xy = xy(:);
    p = A*((A'*xy) * inv(A'*A));
    p = p(:)' + shiftFactor;

    % construct a new cellLine with the projection, p inserted
    ix = sort(ix);
    centerLine = cat(1,CENTERLINE(1:ix(1),:),p,CENTERLINE(ix(2):end,:));
    insertLocation = ix(2);
    
    if sum( sum( (repmat(p,[length(centerLine(:,1)),1])-centerLine).^2 ,2) == 0 ) > 0
        %something wonky happened witht he cellList
        insertLocation = find( sum( (repmat(p,[length(centerLine(:,1)),1])-centerLine).^2, 2) == 0 );
    end
    
    pole2pointLength = sum(sum((centerLine(2:insertLocation,1:2) - centerLine(1:insertLocation-1,1:2)).^2,2).^(1/2));
    
    % the distance from this point is the position along the cell width
    
    dxy = centerLine(insertLocation(1),:) - (xy(:)' + shiftFactor);
    distFromCenterLine = sum(dxy.^2).^(1/2);

    % get the fraction along the cell length
    fractionalLength(q) = (pole2pointLength - halfLength) / halfLength;

    % now cell width
    fractionalWidth(q) = LR*distFromCenterLine/halfWidth;
end
end