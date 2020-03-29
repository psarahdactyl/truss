clf;
% axis equal;
hold on;
% v = readDMAT("v.dmat");
v = readDMAT("vterm.dmat");
vl = readDMAT("vandlterm.dmat");
l = readDMAT("lterm.dmat");
histogram(vl,20,'FaceColor','blue')
histogram(v,20,'FaceColor','magenta')
histogram(l,20,'FaceColor','cyan')
% histogram(vterm,20,'FaceColor','red')

% a = readDMAT("areas.dmat");
% ax = readDMAT("ax.dmat");
% histogram(a,200,'FaceColor','cyan')

