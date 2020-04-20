
%[V,F] = load_mesh('~/Dropbox/models/bunny-remesh.obj');
%V = V * axisangle2matrix([0 0 1],-pi*0.1);
%V = V * axisangle2matrix([-1 1 0],pi*0.1);
%% rescale to unit sphere
%V = V-0.5*(max(V)+min(V));
%V = V/max(normrow(V));
%% 0.3 m wide
%V = V*0.3*0.5;
%AO = [];
CB = cbrewer('Set1',5);
cbred = CB(1,:);
cbblue = CB(2,:);
cbgreen = CB(3,:);
cborange = CB(5,:);

cen = centroid(V,F);
X = [-0.09475 -0.0489 0.02751;-0.04366 0.006452 -0.1024;0.05819 -0.005395 0.03442];
X = V(snap_points(X,V),:);
nx = size(X,1);
E = [1 nx+1;1 nx+2;2 nx+3;2 nx+4;3 nx+5];
Y = [ ...
  X(1,:)+[-1 0 0.5]
  X(1,:)+[-1 0 -0.5]
  X(2,:)+[-1 0 0.3]
  X(2,:)+[-1 0 -0.6]
  X(3,:)+[0.1 0 1]
  ];

[CV,CF] = edge_cylinders([X;Y],E,'Thickness',0.005,'PolySize',30);

fcen = [0 0 -9.8];


clf;
hold on
tsh = tsurf(F,V, ...
  'FaceVertexCData',repmat(cbgreen,size(V,1),1),fphong,fsoft,falpha(0.8,0));
%sct(X,'o','LineWidth',3,'MarkerFaceColor','w','MarkerEdgeColor','k','SizeData',300);
sct(cen,'o','LineWidth',3,'MarkerFaceColor','w','MarkerEdgeColor','k','SizeData',300);
%plot_edges([X;Y],E);
csh = tsurf(CF,CV, ...
  falpha(1,0),fsoft,'FaceVertexCData',repmat(cborange,size(CV,1),1));
%qvr(cen,fcen);
hold off;
camproj('persp');
axis equal;
axis manual;
axis([     -0.25 0.17297     -0.11831      0.11831     -0.2 0.20438]);
view(0,0);
l = add_lights();
AO = apply_ambient_occlusion(tsh,'AddLights',false,'AO',AO);
set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');
