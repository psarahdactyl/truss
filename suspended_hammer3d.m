
% load and rotate mesh
[MV,MF] = load_mesh('hammer.obj');
MV = MV*axisangle2matrix([1 0 0],-pi/2);
MV = MV*axisangle2matrix([0 1 0],-0.4);

% generate colors for each component for visualization
[~,MI] = connected_components(MF);
MC = repmat(0.8*[0.99 1 0.99],size(MF,1),1);
MC((MI==2),1) = 0.6;
MC((MI==2),2) = 0.5;
MC((MI==2),3) = 0.4;


%[WV] = bounding_box(MV(:,[2 3]))
%WV = [repmat(1.5*min(MV(:,1))-0.5*max(MV(:,1)),4,1) WV];
%WF = [1 2 4;1 4 3];

% find bounding box size
bbd = normrow(max(MV)-min(MV));
gr = bbd*0.02;
[Y,Z] = meshgrid( ...
  min(MV(:,2))-bbd*0.0:gr:max(MV(:,2))+bbd*0.0, ...
  min(MV(:,3))-bbd*0.3:gr:max(MV(:,3))+bbd*0.5);
% make wall mesh
[WQ,WV] = surf2patch(0*Y+1.5*min(MV(:,1))-0.5*max(MV(:,1)),Y,Z);
DV = WV;
% find points on side of hammer
% DV is points on wall mesh, repmat is ray directions, MV and MF are
% vertices on hammer
[~,T] = ray_mesh_intersect(DV,repmat([1 0 0],size(DV,1),1),MV,MF);
% RV is the points on hammer
RV = WV(~isinf(T),:)+T(~isinf(T)).*[1 0 0];

%DV = random_points_on_mesh(WV,WF,50);
%RV = random_points_on_mesh(MV,MF,50);

% this creates the ground structure, (V,E), including the center of mass
% and connection points on both wall and hammer
[I,J] = find(ones(size(RV,1),size(DV,1)));
%cen = centroid(MV,MF);
cen =  ...
  0.9*centroid(MV,MF(MI==5,:))+ ...
  0.1*centroid(MV,MF(MI~=5,:));
V = [cen;RV;DV];
E = [ones(size(RV,1),1) 1+(1:size(RV,1))';1+[I size(RV,1)+J]];
%dir = DV - RV;
%[~,T] = ray_mesh_intersect(DV,dir,MV,MF);

% nodal forces, boundary vertices, stress limits
f = zeros(size(V));
f(1,3) = -9.8;
bf = 1+size(RV,1)+(1:size(DV,1));
sC = 5e1;
sT = 5e1;
% get areas, axial forces, lengths, nodal equilibrium matrix
% [a,n,l,BT] = groundstructure(V,E,f,bf,sC,sT,'IgnoredEdges',1:size(RV,1));
% optimization
[B,C] = create_nodal_equilibrium_matrices(V,E,bf);
[A,b,Aeq,beq] = create_constraint_matrices(V,E,f,bf,sC,sT,sB);
lengths = edge_lengths(V,E);
[x,ar,ax,be] = optimize_lp(lengths,A,b,Aeq,beq,'cvx','IgnoredEdges',1:size(RV,1));
% [x,ar,ax,be] = optimize_lp(lengths+EAs,A,b,Aeq,beq,'yalmip','IgnoredEdges',ig);


naa = size(RV,1);
av = 1 + (1:naa)';
ae = (1:naa)';
dim = size(V,2);
fa = reshape(B(av + size(V,1)*[0:dim-1],naa+1:end)*n(naa+1:end),[],dim);
torque = normrow(sum(cross(fa,V(E(ae,2),:)-V(E(ae,1),:),2)))

% find areas bigger than a certain threshold and place a cylinder there
NZ = find(max(a,0)>1e-7);
[CV,CF,CJ,CI] = edge_cylinders(V,E(NZ,:),'PolySize',10,'Thickness',sqrt(max(a(NZ),0)/pi));

fbf = reshape(B(bf(:) + [0:dim-1]*size(V,1),:)*n,[],dim);

clf;
hold on;
% plot wall
tsurf(WQ,WV,'FaceColor',0.5+0.5*blue,falpha(1,0),fsoft);
% plot mesh 
tsurf(MF,MV,'FaceVertexCData',MC,falpha(0.8,0),fsoft);
% plot points on wall (connection points for bipartite graph)
scatter3(DV(:,1),DV(:,2),DV(:,3),   '.b','SizeData',1000);
% plot points on side of hammer (connection points for bipartite graph)
scatter3(RV(:,1),RV(:,2),RV(:,3),   '.k','SizeData',1000);
% plot center of mass
scatter3(cen(:,1),cen(:,2),cen(:,3),'.r','SizeData',1000);

quiver3( ...
  V(bf,1),V(bf,2),V(bf,3), ...
  fbf(:,1),fbf(:,2),fbf(:,3), ...
  0,'.b','LineWidth',3);

quiver3( ...
  V(:,1),V(:,2),V(:,3), ...
  f(:,1),f(:,2),f(:,3),'r','LineWidth',2);
%tsurf(E,V,falpha(0,0.015));

% plot the cylinders
tsurf(CF,CV,falpha(1,0),'CData',sign(n(NZ(CJ))),fsoft);

hold off;
axis equal;
view(26,15);
camlight;
%caxis(max(abs(caxis))*[-1 1])
caxis([-1 1])
colormap(flipud(cbrewer('RdBu',256)))
colorbar
