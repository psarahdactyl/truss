% this file will be using point loads instead of rigid body loads
clf
filename = ...
    "~/Documents/siggraph2020/hidden-supports-application/viewer/data/scenes/canonical.txt";


% read scene and plot the objects
[objs,bb] = read_scene(filename);
% read grid
GV = readDMAT("~/Documents/siggraph2020/truss-opt/grid.dmat");
voxel_width = GV(1,1) - GV(2,1);
GV = GV - abs(repmat(repelem(voxel_width/2,3),size(GV,1),1));
% read visibility values of GV for each object
Vs = readDMAT("~/Documents/siggraph2020/truss-opt/grid_visibilities.dmat");
side = readDMAT("side.dmat");

hold on
axis equal
% plot background wall
tsurf(objs{1}{2},objs{1}{1},...
    'FaceColor',0.5+0.5*blue,falpha(1,0),fsoft)


MC = repmat(0.8*[0.99 1 0.99],size(objs{i}{2},1),1); % colors
%     % plot object
%     tsurf(objs{i}{2},objs{i}{1},...
%         'FaceVertexCData',MC,falpha(0.8,0),fsoft);

DV = objs{1}{1}; % wall mesh

RV = objs{2}{1}; % object mesh
P = centroid(RV,objs{2}{2});

% hidden-ness/visibility
H = Vs;

% this creates the ground structure, (V,E)
[I,J] = find(ones(size(P,1),size(DV,1)));

V = [P;DV];
E = [I' size(P,1)+J'];

% nodal forces, boundary vertices, stress limits
% place nodal force on first V
f = zeros(size(V));
%     f(1:3,2) = -9.8;
f(1,2) = -9.8;
%     V_forces = V(1:3,:);
V_forces = V(1,:);
bf = size(P,1)+(1:size(DV,1));
sC = 5e2;
sT = 5e2;    

%%
% get areas, axial forces, lengths, nodal equilibrium matrix
[a,n,l,h,BT] = groundstructure(V,E,H,GV,f,bf,sC,sT);

% find areas bigger than a certain threshold and place a cylinder there
NZ = find(max(a,0)>1e-7);
num_rods = size(NZ,1)
num_compression = sum(sign(n(NZ))==1)
num_tension = sum(sign(n(NZ))==-1)
[CV,CF,CJ,CI] = edge_cylinders(V,E(NZ,:),...
    'PolySize',10,'Thickness',sqrt(max(a(NZ),0)/pi));

    % plot the cylinders
    color = sign(n(NZ(CJ)));
    tsurf(CF,CV,falpha(1,0),'CData',color,fsoft);

%%

%     % plot visibilities
%     GV_h = GV(H~=0,:);
%     H_h = H(H~=0);
%     scatter3(GV_h(:,1),GV_h(:,2),GV_h(:,3),...
%         '.','CData', abs(H_h),'SizeData',300);

EV = normalizerow(V(E(:,2),:)-V(E(:,1),:));
l = edge_lengths(V,E);
h = edge_visibilities(V,E,GV,H,l);

o = l.*h';

%%
% plot(graph(E(:,1),E(:,2)),...
%     'XData',V(:,1),'YData',V(:,2),'ZData',V(:,3),...
%     'EdgeCData',o, 'LineWidth',2);

%     scatter3(GV(:,1),GV(:,2),GV(:,3),...
%         '.','CData', abs(Vs),'SizeData',30);

quiver3( ...
  V(:,1),V(:,2),V(:,3), ...
  f(:,1),f(:,2),f(:,3),'r','LineWidth',2);

%%
hold off;
axis equal;
% view(26,15);
camlight;
camup([0 1 0]);
% caxis([min(o), max(o)])
cameratoolbar('SetCoordSys','y')
cameratoolbar('setmode','orbit')
camproj('perspective')
caxis([-1 1])
colormap(flipud(cbrewer('RdBu',256)))
title(sprintf('Vol: %g, #bars: %d',l'*a,numel(NZ)),'FontSize',20);
colorbar

