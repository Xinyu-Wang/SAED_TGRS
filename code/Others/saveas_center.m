function saveas_center(h, save_file, width,height)
% saveas_center(h, save_file, width, height)

set(0,'CurrentFigure',h);
set(gcf,'PaperPositionMode','auto');
set(gca,'position',[0,0,1,1]);
set(gcf,'position',[200,200,2*height,2*width]);
% set(gcf,'visible','off')

saveas(h, save_file);