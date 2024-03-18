function plot_extreme_cells(data,paths,n,varargin)

if nargin == 4
    channel = varargin{1};
else
    channel = 'GFP';
end

for ii = 1:size(data,1)
    data{ii}(:,9) = repmat(ii,[size(data{ii},1),1]);
end

data = vertcat(data{:});

data(:,10) = data(:,5) ./ data(:,4);
data(:,11) = data(:,6) ./ data(:,4);
data(:,12) = min(data(:,10),data(:,11));
[~,ix] = sort(data(:,12),'ascend');
data = data(ix,:);
top = data(1:n,:);
bottom = data(end-n+1:end,:);

for ii = 1:n
    subplot(2,n,ii);
    im = helper(top,ii,paths);
    showim(im);
    subplot(2,n,ii+n)
    im = helper(bottom,ii,paths);
    showim(im);
end

    function res = helper(x,ii,paths)
        
        root = '\\172.24.22.178\Data_02\William Gray\William Gray PART-III imaging\Stationary Phase Ribosome Organization\CJW4677EcoliTests\';
        cellList = paths{x(ii,9),4};
        frame = x(ii,1);
        cc = x(ii,2);
        
        Phase_dir = [root,paths{x(ii,9),1},'\1ph\'];
        Phase_list = dir([Phase_dir,'*.tif']);
        im_phase = double(imread(fullfile(Phase_list(frame).folder,Phase_list(frame).name)));
        
        GFP_dir = [root,paths{x(ii,9),1},'\2GFP\'];
        GFP_list = dir([GFP_dir,'*.tif']);
        if isempty(GFP_list)
            GFP_list = dir([[root,paths{x(ii,9),1},'\3GFP\'],'*.tif']);
        end
        im_GFP = double(imread(fullfile(GFP_list(frame).folder,GFP_list(frame).name)));
        
        p_phase = Project_Cell_Signal(cellList.meshData{frame}{cc}.mesh,-5:5,2,im_phase);
        p_phase = imresize(p_phase,[40,10]);
        
        p_GFP = Project_Cell_Signal(cellList.meshData{frame}{cc}.mesh,-5:5,2,im_GFP);
        p_GFP = imresize(p_GFP,[40,10]);
        
        a = sum(sum(p_GFP(1:5,:)));
        b = sum(sum(p_GFP(end-5:end,:)));
        if a > b
            p_GFP = flip(p_GFP,1);
            p_phase = flip(p_phase,1);
        end
        
        if strcmp(channel,'phase')
            res = p_phase;
        else
            res = p_GFP;
        end
    end

end