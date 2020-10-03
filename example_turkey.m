rng(0)
V = {};
F = {};
coms = [nan nan nan];
vols = nan;
masses = nan;
[V{end+1},F{end+1}] = load_mesh('data/meshes/turkey-bending.obj');
[coms(end+1,:),vols(end+1)] = shell_centroid(V{end},F{end});
[V{end},F{end}] = selfintersect(V{end},F{end},'StitchAll',true);
masses = vols*5;

CB = cbrewer('Set1',5);
cbred = CB(1,:);
cbblue = CB(2,:);
cbgreen = CB(3,:);
cborange = CB(5,:);

[SV,SF] = load_mesh('data/meshes/building.obj');
% [SV,SF] = load_mesh('data/meshes/building-both.obj');

% units are meters

% [X,Y] = meshgrid(linspace(0,1),linspace(0,1));
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
% % NORMAL GROUND STRUCTURE CODE
% [XX,XE,XC,YX,YE,YC] = groundstructure(VV,FF,CC,n);
% XE(XC(XE(:,1))==1&XC(XE(:,2))==3,:) = [];
% r = 0.001;
% Xvis = groundstructure_visibility(Q,VV,FF,XX,XE,'SampleSize',r);

% %%%%%%%%
% % MAKING A GROUND STRUCTURE 
% % BY PROJECTING THE COM ONTO SUPPORT SURFACE
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
% OS=wps(1:sw/2,:);
% wps((sw/2+1):sw,:) = [OS(:,1), -OS(:,2), OS(:,3)];

XX=[nps(nnz,:);wps];
XE=[(1:sn)' ((sn+1):(sn+sw))'];
XEV = XX(XE(:,2),:)-XX(XE(:,1),:);
XC=[ones(sn,1);2*ones(sw,1)];
% 
%%%%%%%%%%%%

XE=[XE;XE+size(XX,1)];
XX=[XX;[XX(:,1) -XX(:,2) XX(:,3)]];
XC=[XC;XC];
%%
g = [0 0 -9.8];
sC = ones(size(XE,1),1)*1e5;
sT = ones(size(XE,1),1)*1e5;
sB = ones(size(XE,1),1)*1e2;

[A,b,Aeq,beq] = create_constraint_matrices(XX,XE,XC,coms,g,sC,sT,sB,'Mass',masses);
objective = edge_lengths(XX,XE);% + 10000*Xvis.^2;
[x,ar,ax,be,fval] = optimize_lp(objective,A,b,Aeq,beq,'linprog');

NZ = ar>1e-6;
[CV,CF] = edge_cylinders(XX,XE(NZ,:),'Thickness',sqrt(ar(NZ)),'PolySize',30);

%%
clf;
hold on;
SF=[SF;SF+size(SV,1)];
SV=[SV;[SV(:,1),-SV(:,2),SV(:,3)]];
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
    fphong,falpha(1,0),fsoft);
end
csh = tsurf(CF,CV, ...
  falpha(1,0),fsoft,'FaceVertexCData',repmat(cborange,size(CF,1),1));
hold off;
axis equal;
% view(-58,15);
% view(90,0);
view(90,31);

l = { ...
 light('Color',0.5*[1 1 1],'Position',(campos-camtarget),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position',-(campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],pi*0.9),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],-pi*0.9),'Style','local') ...
 light('Color',0.25*[1 1 1],'Position', [-1 -1 200],'Style','local') ...
};
apply_ambient_occlusion([tsh{:},ssh,csh],'AddLights',false);apply_ambient_occlusion(csh,'AddLights',false,'Factor',0.5);
dsh = add_shadow({csh,ssh,tsh{:}},l{6});
%apply_ambient_occlusion([],'AddLights',false);
camproj('persp');
%view(0,0);
set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');

%for t = linspace(0,2*pi,60);camtarget([0.5 1 0.5]);camup([0 0 1]);campos(mean(views)+0.05*[((1+pi.*cos(t)).*cos(t)-pi+1)/2 0 0.5*(1+pi.*cos(t)).*sin(t)/2]);camproj('persp');camva(50);drawnow;end
