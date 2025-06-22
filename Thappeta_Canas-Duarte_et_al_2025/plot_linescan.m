function plot_linescan(cellList,frame,cc)

fields = fieldnames(cellList.meshData{frame}{cc});

sigs = cell(length(fields),1);
labels = cell(length(fields),1);

for ii = 1:length(fields)
    if ~isempty(regexp(fields{ii},regexptranslate('wildcard','signal*'), 'once'))
        sigs{ii} = cellList.meshData{frame}{cc}.(fields{ii});
        sigs{ii} = sigs{ii} ./ sum(sigs{ii});
        labels{ii} = fields{ii};
    end
end

remove = cellfun(@isempty,sigs);
sigs(remove) = [];
labels(remove) = [];
xs = cumsum(cellList.meshData{frame}{cc}.steplength);
xs = xs ./ xs(end);

cs = lines(length(sigs));
% figure;
hold on;
for ii = 1:length(sigs)
   plot(xs,sigs{ii},'color',cs(ii,:),'linewidth',2);
end
xlabel('Normalized cell length','fontsize',12);
ylabel('Normalized signal','fontsize',12);
legend({cell2mat(labels)},'fontsize',12);
end