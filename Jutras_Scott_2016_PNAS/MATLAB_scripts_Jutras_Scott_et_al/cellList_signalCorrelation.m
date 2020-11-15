% Assumes that a cellList is available in the workspace with Bocillin as signal1
% and HADA as signal2
bocillinSignal = 'signal1';
hadaSignal = 'signal2';
nm = @(x) x / sum(x);

bocilin = [];
hada = [];
pole = 3;
for F = 1:length(cellList.meshData)
    for C = 1:length(cellList.meshData{F})
        if isempty(cellList.meshData{F}{C}) || ~isfield(cellList.meshData{F}{C},'mesh') || length(cellList.meshData{F}{C}.mesh)<6
            continue
        end
        
        A = getMeshSegArea(cellList.meshData{F}{C}.mesh);
        S1 = nm(cellList.meshData{F}{C}.(bocillinSignal)) ./ A;
        S2 = nm(cellList.meshData{F}{C}.(hadaSignal)) ./ A;
        
        bocilin = [bocilin, S1(pole:end-(pole-1))'];
        hada = [hada, S2(pole:end-(pole-1))'];
         
    end
end

xy = cat(2,bocilin(:),hada(:));
bins = 500;
[n,c] = hist3(fliplr(xy),[bins,bins]);
posx = repmat(c{1},[length(c{1}), 1]);
posy = repmat(c{2}',[1,length(c{1})]);
posx = posx(:);
posy = posy(:);
n = n(:);

keep = n > 0;
n = log(n(keep));
posx = posx(keep);
posy = posy(keep);

[n,ix] = sort(n);
posx = posx(ix);
posy = posy(ix);

nCols = 100;
cBins = [-1, 0 : (max(n) / nCols) : max(n)]; % >, <=
if cBins(end) ~= max(n), cBins(end+1) = max(n); end
cc = (jet(length(cBins)));

figure
hold on
for k = 1:length(n)
    ix = find(n(k) < cBins, 1);
    if isempty(ix),continue,end
    plot(posx(k),posy(k),'o','color',cc(ix,:),'markerfacecolor',cc(ix,:))
end

