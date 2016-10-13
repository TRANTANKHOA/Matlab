%% Resize graph
w = 1;h = 0.35;
f=gcf;undock(f)
w=round(w*1050);h=round(h*800);
sz=[-1000 1000 w h];
set(f,'Units','pixels','Position',sz);drawnow
%% Printing
f=gcf;
% set(f,'Units','pixels','Position',sz);drawnow
r=300;
[name, path] = uiputfile({'*.pdf','PDF';'*.eps','Vector based .eps';'*.jpg;*.tif;*.png;*.gif','All Image Files';'*.*','All Files' },'Save Image','/Users/erictran/Dropbox/PhD/Reports/MCMC4SISO/');
export_fig(fullfile(path,name),f,'-a1',['-r',num2str(r)],'-nocrop');%, '-grey'
saveas(f,fullfile(path,[name(1:end-4) '.fig']))
%% Miscelaneous
set(gca,'ycolor','k','xcolor','k')
set(gca,'ycolor','w','xcolor','w')
%% Moving plot box and legend
[hleg1, hobj1] = legend(gca);
% textobj = findobj(hobj1, 'type', 'text');
% set(textobj, 'Interpreter', 'latex', 'fontsize', font_size);
leg_box = get(hleg1,'position');
plt_box = get(gca,'position');
plt_box(1) = plt_box(1)+.05*plt_box(3);
set(gca,'position',plt_box)

leg_box(1) = leg_box(1)-.05*leg_box(3);
leg_box(3) = 1.05*leg_box(3); 
legend('Location',leg_box)
% legend('Location','NorthOutside')
% legend({'$\tilde{\phi}(\theta)$','$\rho\phi_\circ(\theta)$'},'Interpreter','latex')
% set(hleg1, 'Color', 'white')
% axis tight
% legend('boxon')
%%Delete markers
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
%% Create zoom out window
kids = get(gca,'Children'); % Collect handles to lines data
xd = get(kids,'XData'); % Extract lines data on X
yd = get(kids,'YData'); % Extract lines data on Y
linestyle = get(kids,'linestyle'); % Extract linestyle
color = get(kids,'color'); % Extract line colors
ax1 = gca;axis tight
ax2 = axes('Position',get(ax1,'Position'),...
    'Color','w',...
    'XColor','k','YColor','k',...
    'xticklabel',[],'yticklabel',[]...
    ); hold on
for i=1:length(kids)
plot(xd{i},yd{i},'linestyle',linestyle{i},'color',color{i})
end
axis tight

set(ax2,'xlim',[-.5 .5])% axes limit for zoom out,
set(ax2,'ylim',[.325 .5])% axes limit for zoom out,
set(ax2, 'Units', 'normalized', 'Position', [0.2 0.2 0.35 0.45])
box on