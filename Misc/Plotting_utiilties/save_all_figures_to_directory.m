function save_all_figures_to_directory(base_name)
figlist=findobj('type','figure');
for i=1:numel(figlist)
    saveas(figlist(i),fullfile(cd,[base_name '_' num2str(figlist(i).Number) '.fig']));
end
end