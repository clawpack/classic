colormap(jet);
axis([-1 1 -1 1 -1 1])
caxis([-0.003 0.003])
grid on

hideslices;
showslices('x',6);
showslices('y',6);
showslices('z',4);

% plot dashed line showing where interface is:
hold on;
[xcm,ycm] = meshgrid([-1 1],[-1 1]);
p = patch(surf2patch(xcm,ycm,0*xcm));
set(p,'cdata',0*xcm);
set(p,'facecolor','w');
set(p,'facealpha',0.7);
set(p,'edgecolor','k');

plot3([-1 1],[0 0],[0 0],'w','linewidth',3)
plot3([0 0],[-1 1],[0 0],'w','linewidth',3)
hold off

showpatchborders;
setpatchbordercolor('k');

cv = linspace(-0.005,0.005,31);
cv(cv == 0) = [];
drawcontourlines(cv);
setcontourlineprops('color','k');

colorbar;

shg;

clear afterframe;
