% this file will be using point loads instead of rigid body loads
clf
clear
close all
hold on
axis equal

%%
%%%%%%%%%% FILE IO
filename = ...
    "~/Documents/siggraph2020/hidden-supports-application/viewer/data/scenes/canonical.txt";


% read scene and plot the objects
[objs,bb] = read_scene(filename);
[AV,AF] = list_to_mesh(objs);

%%
%%%%%%%%%%%% SET UP GROUND STRUCTURE

DF = objs{1}{2};
DV = objs{1}{1};
% plot background wall
tsurf(DF,DV,...
    'FaceColor',0.5+0.5*blue,falpha(1,0),fsoft)

MC = repmat(0.8*[0.99 1 0.99],size(objs{2}{2},1),1); % colors

RV = objs{2}{1}; % object mesh
RF = objs{2}{2};
% % plot object
% tsurf(RF,RV,...
%     'FaceVertexCData',MC,falpha(0.8,0),fsoft);

% sample points on surfaces
Rc = RV(3,:);
Dc = DV;
scatter3(Rc(:,1),Rc(:,2),Rc(:,3),'.g','SizeData',300);
scatter3(Dc(:,1),Dc(:,2),Dc(:,3),'.y','SizeData',300);

% this creates the ground structure, (V,E)
[I,J] = find(ones(size(Rc,1),size(Dc,1)));

V = [Rc;Dc];
E = [I' size(Rc,1)+J'];

% nodal forces, boundary vertices
f = zeros(size(V));
bf = size(Rc,1)+(1:size(Dc,1));
f(1,2) = -9.8;

%%
%%%%%%%%%% LENGTHS
% bar lengths
lengths = edge_lengths(V,E);

%%
%%%%%%%%% VISIBILITY

% create voxel grid
[GV,side,w] = voxel_grid(AV,20);

% generate views
[gx,gy] = meshgrid([-3:-.5:-5],[3:.5:5]);
gx = gx(:);
gy = gy(:);
views = [gx, gy, 4*ones(length(gx),1)];
scatter3(views(:,1),views(:,2),views(:,3),'.m','SizeData',100);

% calculate visibilities
src = repmat(views,size(GV,1),1);
dir = repelem(GV,size(views,1),1);
Vs = ray_mesh_intersect(src,dir-src,AV,AF);
Vs(Vs ~= 0) = 1;
Vs = sum(reshape(Vs,size(views,1),size(GV,1)),1);
Vs = Vs ./ size(views,1);

% bar projected areas
% average distances from viewpoints to bar midpoints
midpoints = (V(E(:,1),:) + V(E(:,2),:))/2;
vs = repmat(views,size(E,1),1);
ms = repelem(midpoints,size(views,1),1);
d = vecnorm(ms-vs,2,2);
d = mean(reshape(d,size(views,1),size(E,1)),1)';

% bar thetas
edges = repelem([V(E(:,1),:) V(E(:,2),:)],size(views,1),1);
v1s = edges(:,1:3)-vs;
v2s = edges(:,4:6)-vs;
normcostheta = dot(v1s,v2s,2);
costheta = normcostheta ./ (vecnorm(v1s,2,2).*vecnorm(v2s,2,2));
meantheta = mean(reshape(costheta,size(views,1),size(E,1)),1)';
theta = acos(meantheta) * (180 / pi );

% bar visibilities
[v, segments] = edge_visibilities(V,E,GV,side,w,Vs,lengths);

pva = theta.^2 .* v'.^2 ./ d.^2;


%%
%%%%%%%%%%% BAR PROPERTIES
sC = 1e2;
sT = 1e2;
sB = 1e2;

dim=3;
n = size(V,1);
m = size(E,1);

%%
%%%%%%%%%%% NODAL EQUILIBRIUM MATRICES
EV = normalizerow(V(E(:,2),:)-V(E(:,1),:)); % B matrix
[fr,fc] = find(f~=0);
fsum = sum(f(fr,:),1);
% PV = fsum.*ones(size(EV,1),3);
% return
% PV = normalizerow(PV); % C matrix
% PV = [zeros(size(EV,1),1) ones(size(EV,1),1) zeros(size(EV,1),1)]; % C matrix
PV = [EV(:,2) -EV(:,1) zeros(size(EV,1),1)]; % C matrix

% MV = (V(E(:,2),:)+V(E(:,1),:))/2; % midpoints for plotting

BT = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[EV;-EV],dim*n,m);
CT = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[PV;-PV],dim*n,m);

% remove fixed vertices
bf = bf(:);
BT(bf+n*(0:dim-1),:) = [];
CT(bf+n*(0:dim-1),:) = [];
ff = f;
ff(bf,:,:) = [];

B = [sparse(size(BT,1),m) BT sparse(size(BT,1),m)];
Bnb = [sparse(size(BT,1),m) BT];
C = [sparse(size(CT,1),2*m) CT];

%%
%%%%%%%%%%% INEQUALITY CONSTRAINT MATRICES
II = speye(m,m);
Z = sparse(m,m);
nf = size(f,3);

Anb = [repmat(-sC*II,nf,1), -II; ...
    repmat(-sT*II,nf,1),  II];
bnb = zeros(2*m*nf,1);

A = [repmat(-sC*II,nf,1), -II, Z; ...
    repmat(-sT*II,nf,1) , II, Z;...
    repmat(-sB*II,nf,1), Z, -II;...
    repmat(-sB*II,nf,1), Z, -II];
b = zeros(4*m*nf,1);


%%
%%%%%%%%% OPTIMIZATION
objective = lengths .* pva;

% %%%%%%% without bending %%%%%%%%
% [an,~,flags,output] = linprog([objective;zeros(m*nf,1)],Anb,bnb,...
%   Bnb,ff(:), ...
%   [zeros(m,1);-inf(m,1)], ...
%   [inf*ones(m,1);inf(m,1)]);
% ar = an(1:m);
% ax = reshape(an(m+(1:m*nf)),m,nf);

% %%%%%%% with bending %%%%%%%%
[anb,~,flags,output] = linprog([lengths;zeros(2*m*nf,1)],A,b,...
  B+C,ff(:), ...
  [zeros(m,1);-inf(2*m,1)], ...
  [inf*ones(m,1);inf(2*m,1)]);
ar = anb(1:m);
ax = reshape(anb(m+(1:m*nf)),m,nf);
be = reshape(anb(2*m+(1:m*nf)),m,nf);

%%
%%%%%%%%%%%% PLOT
% find areas bigger than a certain threshold and place a cylinder there
NZ = find(max(ar,0)>1e-7);
num_rods = size(NZ,1)
num_compression = sum(sign(ax(NZ))==1)
num_tension = sum(sign(ax(NZ))==-1)
[CV,CF,CJ,CI] = edge_cylinders(V,E(NZ,:),...
    'PolySize',10,'Thickness',sqrt(max(ar(NZ),0)/pi));

% plot the cylinders
color = sign(ax(NZ(CJ)));
tsurf(CF,CV,falpha(1,0),'CData',color,fsoft);


%%
hold on
[~, II] = sortrows(E);

% plot(graph(E(:,1),E(:,2)),...
%     'XData',V(:,1),'YData',V(:,2),'ZData',V(:,3),...
%     'EdgeCData',lengths, 'LineWidth',2);

% caxis([min(lengths) max(lengths)])
% colormap(flipud(cbrewer('RdBu',256)))
% colorbar

quiver3( ...
  V(:,1),V(:,2),V(:,3), ...
  f(:,1),f(:,2),f(:,3),'r','LineWidth',2);

% quiver3(...
%     MV(:,1),MV(:,2),MV(:,3), ...
%     PV(:,1),PV(:,2),PV(:,3),0.5,'g','LineWidth',3);
% 
% quiver3(...
%     MV(:,1),MV(:,2),MV(:,3), ...
%     EV(:,1),EV(:,2),EV(:,3),0.5,'c','LineWidth',3);

%%
hold off;
axis equal;
camlight;
camup([0 1 0]);
campos([-33 6 20])
cameratoolbar('SetCoordSys','y')
cameratoolbar('setmode','orbit')
camproj('perspective')
caxis([-1 1])
colormap(flipud(cbrewer('RdBu',256)))
title(sprintf('Vol: %g, #bars: %d out of %i edges',lengths'*ar,numel(NZ),size(E,1)),'FontSize',20);
colorbar

