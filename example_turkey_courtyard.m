rng(0)
addpath utils
V = {};
F = {};
coms = [nan nan nan];
vols = nan;
masses = nan;
[V{end+1},F{end+1}] = load_mesh('data/meshes/turkey.obj');
[coms(end+1,:),vols(end+1)] = shell_centroid(V{end},F{end});
[V{end},F{end}] = selfintersect(V{end},F{end},'StitchAll',true);
masses = vols*5;

CB = cbrewer('Set1',5);
cbred = CB(1,:);
cbblue = CB(2,:);
cbgreen = CB(3,:);
cborange = CB(5,:);

[SV,SF] = load_mesh('data/meshes/courtyard-simple.obj');

[TV,TF] = load_mesh('data/meshes/viewpoints.obj');
Q = blue_noise(300,TV,TF);

% units are meters

% [X,Y] = meshgrid(linspacace(0,1),linspace(0,1));
% PP = normrow([X(:) Y(:)]-mean([X(:) Y(:)]));
% PP = clamp(max(PP)-PP,0.2,1);
% PP = matrixnormalize(PP);
% PP = reshape(PP,size(X));
% PX = X*0.75+0.25*0.5;
% PY = 0*Y-0.75;
% PZ = Y*0.5+0.25;
% [PF,PV] = surf2patch(PX,PY,PZ,'triangles');
% nu = 200;
% %Q = blue_noise(nu*300,PV,PF);
% Q = random_points_on_mesh(PV,PF,nu*30);
% QP = interp2(PX,PZ,PP,Q(:,1),Q(:,3));
% [~,QI] = histc(rand(nu,1),[0;cumsum(QP)]/sum(QP));
% Q = Q(QI,:);

n = 1000;

VV = cell2mat(reshape(V,[],1));
nV = cumsum([0 cellfun(@(V) size(V,1),{V{:}})]);
FF = cell2mat(arrayfun(@(i) nV(i)+F{i},1:numel(V),'UniformOutput',false)');
CC = cell2mat(arrayfun(@(i) repmat(i,size(F{i},1),1),1:numel(F),'UniformOutput',false)');
VV = [SV;VV];
FF = [SF;size(SV,1)+FF];
CC = [ones(size(SF,1),1);1+CC];

%%%%%%%
% NORMAL GROUND STRUCTURE CODE
% [XX,XE,XC,YX,YE,YC] = groundstructure(VV,FF,CC,n);
% XE(XC(XE(:,1))==1&XC(XE(:,2))==3,:) = [];
% r = 0.01;
% Xvis = groundstructure_visibility(Q,VV,FF,XX,XE,'SampleSize',r);

%%%%%%%%
% MAKING A GROUND STRUCTURE 
% BY PROJECTING THE COM ONTO SUPPORT SURFACE
center=coms(2,:);
theta = 0:pi/1000:2*pi;
xs = 40 * cos(theta) + center(1);
ys = 40 * sin(theta) + center(2);
ne=size(theta,2);
srcs=repmat(coms(2,:),ne,1);
dirs=[xs',ys',coms(2,3)*ones(ne,1)];
[flag, t, lambda] = ray_mesh_intersect(srcs, dirs-srcs, VV, FF);
vve=[coms(2,:); dirs];
ee=[ones(ne-1,1) (2:ne)'];

nps=VV(FF(flag,1),:).*lambda(:,1)+...
  VV(FF(flag,2),:).*lambda(:,2)+...
  VV(FF(flag,3),:).*lambda(:,3);

[flagn, tn, lambdan] = ray_mesh_intersect(nps, dirs-nps, SV, SF);

flagnn=flagn(flagn~=0);
lambdann=lambdan(flagn~=0,:);
wps=SV(SF(flagnn,1),:).*lambdann(:,1)+...
  SV(SF(flagnn,2),:).*lambdann(:,2)+...
  SV(SF(flagnn,3),:).*lambdann(:,3);

nnz=find(flagn~=0);
sn=size(nps(nnz,:),1);
sw=size(wps,1);
XX=[nps(nnz,:);wps];
XE=[(1:sn)' ((sn+1):(sn+sw))'];
XEV = XX(XE(:,2),:)-XX(XE(:,1),:);
XC=[ones(sn,1);2*ones(sw,1)];

% discard if intersecting something else in the scene
[H,T] = ray_mesh_intersect(XX(XE(:,1),:),XEV,VV,FF);
H = H & T<0.999999;
XE =   XE(~H,:);

r = 0.01;
Xvis = groundstructure_visibility(Q,VV,FF,XX,XE,'SampleSize',r);
%%%%%%%%%%%%

g = [0 0 -9.8];
wires = false;
wvis = 10000;
%wvis = 0;
if wires
  XE = [XE;XE];
  Xvis = [Xvis;0*Xvis];
  % tension wires are second half
  sC = [ones(size(XE,1)/2,1) zeros(size(XE,1)/2,1)]*1e3;
  sT = [zeros(size(XE,1)/2,1) ones(size(XE,1)/2,1)]*1e3;
else
  sC = ones(size(XE,1),1)*1e3;
  sT = sC;
end
sB = ones(size(XE,1),1)*1e3;

[A,b,Aeq,beq] = create_constraint_matrices(XX,XE,XC,coms,g,sC,sT,sB,'Mass',masses);
objective = edge_lengths(XX,XE) + wvis*Xvis.^2;
[x,ar,ax,be,fval] = optimize_lp(objective,A,b,Aeq,beq,'linprog');

NZ = find(max(ar,0)>1e-7);
[ar(NZ) ax(NZ) be(NZ,:)]
if wires
  RZ = ar>1e-6 & ax>0;
  WZ = ar>1e-6 & ax<0;
else
  RZ = NZ;
  WZ = [];
end
[CV,CF] = edge_cylinders(XX,XE(RZ,:),'Thickness',1.*sqrt(ar(RZ)),'PolySize',30);


clf;
hold on;
scatter3(Q(:,1),Q(:,2),Q(:,3),'.r','SizeData',300);

ssh = tsurf(SF,SV,'FaceVertexCData',repmat(cbblue,size(SV,1),1),fphong,falpha(1,0),fsoft);
red = [1 0 0];
% psh = surf(PX,PY,PZ, ...
%   'CData',PP,'AlphaData',PP,fphong,'EdgeColor','none','FaceAlpha','interp');
% %psh = tsurf(PF,PV,'CData',PP,fphong,falpha(1,0),fsoft);
% colormap(interp1([1 0],[cbred;1 1 1],linspace(0,1,7)));
% psh.DiffuseStrength = 0;
% psh.AmbientStrength = 1;
% psh.SpecularStrength = 0;
tsh = {};
CM = interp1([0 1],[cbgreen;cbgreen*0.6+0.4],linspace(0,1,numel(V))');
for ii = 1:numel(V)
  tsh{end+1} = tsurf( ...
    F{ii},V{ii},'FaceVertexCData',repmat(CM(ii,:),size(V{ii},1),1), ...
    fphong,falpha(0.5,0),fsoft);
end
csh = tsurf(CF,CV, ...
  falpha(1,0),fsoft,'FaceVertexCData',repmat(cborange,size(CF,1),1));
wsh = tsurf(XE(WZ,:),XX,'LineWidth',1,'EdgeColor','k','FaceColor','none');

axis equal;
% view(-58,15);
view(180,0);

l = { ...
 light('Color',0.5*[1 1 1],'Position',(campos-camtarget),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position',-(campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],pi*0.9),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],-pi*0.9),'Style','local') ...
 light('Color',0.25*[1 1 1],'Position', [-1 -1 10],'Style','local') ...
};
apply_ambient_occlusion([tsh{:},ssh],'AddLights',false);apply_ambient_occlusion(csh,'AddLights',false,'Factor',0.5);
% add_shadow({ssh,tsh{:}},l{6});
%apply_ambient_occlusion([],'AddLights',false);
camproj('persp');
%view(0,0);
set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');

%for t = linspace(0,2*pi,60);camtarget([0.5 1 0.5]);camup([0 0 1]);campos(mean(views)+0.05*[((1+pi.*cos(t)).*cos(t)-pi+1)/2 0 0.5*(1+pi.*cos(t)).*sin(t)/2]);camproj('persp');camva(50);drawnow;end
