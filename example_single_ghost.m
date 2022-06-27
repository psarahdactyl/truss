rng(0)
clf
V = {};
F = {};

VI = {};
FI = {};
coms = [nan nan nan];
[V{end+1},F{end+1}] = load_mesh('data/meshes/ghost/ghost-head.obj');
V{end} = V{end} *  ...
  axisangle2matrix([0 0 1],-pi/2);
nv=size(V{end},1);
% V{end} = V{end} + repmat([5 25 0], nv, 1);
V{end} = V{end}/100;
coms(end+1,:) = centroid(V{end},F{end});
masses = [nan 0.035];

[VI{end+1},FI{end+1}] = load_mesh('data/meshes/ghost/ghost-tail.obj');
VI{end} = VI{end} *  ...
  axisangle2matrix([0 0 1],-pi/2);
nv=size(VI{end},1);
% VI{end} = VI{end} + repmat([5 25 0], nv, 1);
VI{end} = VI{end}/100;

CB = cbrewer('Set1',5);
cbred = CB(1,:);
cbblue = CB(2,:);
cbgreen = CB(3,:);
cborange = CB(5,:);

[SV,SF] = load_mesh('data/meshes/ghost/wall-corner.obj');
[SV,SF]=clean_mesh(SV,SF);
SV = SV/100;
SV = SV *  ...
  axisangle2matrix([0 0 1],-pi/2);

%%
% units are meters
[X,Y] = meshgrid(linspace(0,1),linspace(0,1));
PP = normrow([X(:) Y(:)]-mean([X(:) Y(:)]));
PP = clamp(max(PP)-PP,0.2,1);
PP = matrixnormalize(PP);
PP = reshape(PP,size(X));
PX = X*1-0.15;
PY = 0*Y-1;
PZ = Y*0.50-0.01;
[PF,PV] = surf2patch(PX,PY,PZ,'triangles');
nu = 200;
%Q = cy_blue_noise(nu*300,PV,PF);
Q = random_points_on_mesh(PV,PF,nu*30);
QP = interp2(PX,PZ,PP,Q(:,1),Q(:,3));
[~,QI] = histc(rand(nu,1),[0;cumsum(QP)]/sum(QP));
Q = Q(QI,:);

%%
VV = cell2mat(reshape(V,[],1));
nV = cumsum([0 cellfun(@(V) size(V,1),{V{:}})]);
FF = cell2mat(arrayfun(@(i) nV(i)+F{i},1:numel(V),'UniformOutput',false)');
CC = cell2mat(arrayfun(@(i) repmat(i,size(F{i},1),1),1:numel(F),'UniformOutput',false)');
VV = [SV;VV];
FF = [SF;size(SV,1)+FF];
CC = [ones(size(SF,1),1);1+CC];

VVI = cell2mat(reshape(VI,[],1));
nVI = cumsum([0 cellfun(@(VI) size(VI,1),{VI{:}})]);
FFI = cell2mat(arrayfun(@(i) nVI(i)+FI{i},1:numel(VI),'UniformOutput',false)');
CCI = cell2mat(arrayfun(@(i) repmat(i,size(FI{i},1),1),1:numel(FI),'UniformOutput',false)');
VVI = [SV;VVI];
FFI = [SF;size(SV,1)+FFI];
CCI = [ones(size(SF,1),1);1+CCI];

n = 1000;
[XX,XE,XC,YX,YE,YC] = groundstructure(VV,FF,CC,n);
r = 0.001;
% % don't prune
% XX=YX;
% XE=YE;
% XC=YC;
% Xvis = groundstructure_visibility(Q,VV,FF,XX,XE,'SampleSize',r);
Xvis = groundstructure_visibility(Q,VVI,FFI,XX,XE,'SampleSize',r);

%%
% [SXE,SEI]=sortrows(XE);
% XX=XX(1:max(max(XE)),:);
% plot(graph(XE(:,1),XE(:,2)),...
%     'XData',XX(:,1),'YData',XX(:,2),'ZData',XX(:,3),...
%     'EdgeCData',Xvis(SEI), 'LineWidth',2);

%%
force = [0 0 -9.8];
sC = ones(size(XE,1),1)*1e4;
sT = ones(size(XE,1),1)*1e4;
sB = ones(size(XE,1),1)*1e4;

[A,b,Aeq,beq] = create_constraint_matrices(XX,XE,XC,coms,force,...
  sC,sT,sB,'Mass',masses);
objective = edge_lengths(XX,XE) + 100*Xvis.^2;
[x,ar,ax,be,fval] = optimize_lp(objective,A,b,Aeq,beq,'linprog');

NZ = ar>1e-6; % threshold
R=1/pi*sqrt(ar(NZ)); % find radii
R(R<0.005) = 0.005; % round up to rod size
[CV,CF] = edge_cylinders(XX,XE(NZ,:),'Thickness',R,'PolySize',30);

%%
clf;
hold on;
ssh = tsurf(SF,SV,'FaceVertexCData',repmat(cbblue,size(SV,1),1),fphong,falpha(1,0),fsoft);
red = [1 0 0];
psh = surf(PX,PY,PZ, ...
  'CData',PP,'AlphaData',PP,fphong,'EdgeColor','none','FaceAlpha','interp');
%psh = tsurf(PF,PV,'CData',PP,fphong,falpha(1,0),fsoft);
colormap(interp1([1 0],[cbred;1 1 1],linspace(0,1,7)));
psh.DiffuseStrength = 0;
psh.AmbientStrength = 1;
psh.SpecularStrength = 0;
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
view(-58,15);

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
% add_shadow({ssh,tsh{:}},l{6});
% apply_ambient_occlusion([],'AddLights',false);
AO = apply_ambient_occlusion([ssh,tsh{:},csh],'AddLights',false);
add_shadow([ssh,tsh{:},csh],l{6});
camproj('persp');
%view(0,0);
view(-35.419,19.828);
set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');

%%
hold on
tshv = {};
for ii = 1:numel(VI)
  tshv{end+1} = tsurf( ...
    FI{ii},VI{ii},'FaceVertexCData',repmat(CM(ii,:),size(VI{ii},1),1), ...
    fphong,falpha(1,0),fsoft);
end

%for t = linspace(0,2*pi,60);camtarget([0.5 1 0.5]);camup([0 0 1]);campos(mean(views)+0.05*[((1+pi.*cos(t)).*cos(t)-pi+1)/2 0 0.5*(1+pi.*cos(t)).*sin(t)/2]);camproj('persp');camva(50);drawnow;end
