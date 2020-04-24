rng(0)
addpath utils
V = {};
F = {};
coms = [nan nan nan];
vols = nan;
masses = nan;

fid = fopen('data/teaser.txt','r');
txt = textscan(fid,'%s','delimiter',"\n");
fclose(fid);
for i=1:size(txt{1},1)
  txt{1}{i}
  [V{end+1},F{end+1}] = load_mesh(txt{1}{i});
  [coms(end+1,:),vols(end+1,:)] = shell_centroid(V{end},F{end});
  [V{end},F{end}] = selfintersect(V{end},F{end},'StitchAll',true);
end
masses = vols*10;

CB = cbrewer('Set1',5);
cbred = CB(1,:);
cbblue = CB(2,:);
cbgreen = CB(3,:);
cborange = CB(5,:);

[SV,SF] = load_mesh('data/meshes/teaser/wall.obj');

% units are meters

[X,Y] = meshgrid(linspace(0,1),linspace(0,1));
PP = normrow([X(:) Y(:)]-mean([X(:) Y(:)]));
PP = clamp(max(PP)-PP,0.2,1);
PP = matrixnormalize(PP);
PP = reshape(PP,size(X));
PX = X*0.9+0.5;
PY = 0*Y+4.5;
PZ = Y*0.5-0.25;
[PF,PV] = surf2patch(PX,PY,PZ,'triangles');
nu = 200;
%Q = blue_noise(nu*300,PV,PF);
Q = random_points_on_mesh(PV,PF,nu*30);
QP = interp2(PX,PZ,PP,Q(:,1),Q(:,3));
[~,QI] = histc(rand(nu,1),[0;cumsum(QP)]/sum(QP));
Q = Q(QI,:);
Q = Q *  ...
  axisangle2matrix([0 0 1],-pi/4);

QX = X*0.75-0.25;
QY = 0*Y-1.5;
QZ = Y*1.5-0.25;
[RF,RV] = surf2patch(QX,QY,QZ,'triangles');
R = random_points_on_mesh(RV,RF,nu*30);
RP = interp2(QX,QZ,PP,R(:,1),R(:,3));
[~,RI] = histc(rand(nu,1),[0;cumsum(RP)]/sum(RP));
R = R(RI,:);

n = 1500;

VV = cell2mat(reshape(V,[],1));
nV = cumsum([0 cellfun(@(V) size(V,1),{V{:}})]);
FF = cell2mat(arrayfun(@(i) nV(i)+F{i},1:numel(V),'UniformOutput',false)');
CC = cell2mat(arrayfun(@(i) repmat(i,size(F{i},1),1),1:numel(F),'UniformOutput',false)');
VV = [SV;VV];
FF = [SF;size(SV,1)+FF];
CC = [ones(size(SF,1),1);1+CC];
[XX,XE,XC,YX,YE,YC] = groundstructure(VV,FF,CC,n);
% % 19 (ring) , 18 (planet)
% % XE(XC(XE(:,1))==1&XC(XE(:,2))==3,:) = [];
s1=find(ismember(XE(:,1),find(XC==11)));
s2=find(ismember(XE(:,2),find(XC==11)));
e1=find(ismember(XE(:,1),find(XC==12)));
e2=find(ismember(XE(:,2),find(XC==12)));
set_diff=setdiff([s1;s2],[e1;e2])
XE(setdiff([s1;s2],[e1;e2]),:) = [];

%%
r = 0.01;
Xvis = groundstructure_visibility([Q;R],VV,FF,XX,XE,'SampleSize',r);
% Xvis = groundstructure_visibility(Q,VV,FF,XX,XE,'SampleSize',r);
% XE=YE;

%%
g = [0 0 -9.8];
wires = true;
wvis = 10000;
%wvis = 0;
if wires
  XE = [XE;XE];
  Xvis = [Xvis;0*Xvis];
  % tension wires are second half
  sC = [ones(size(XE,1)/2,1) zeros(size(XE,1)/2,1)]*1e5;
  sT = [zeros(size(XE,1)/2,1) ones(size(XE,1)/2,1)]*1e5;
  sB = [ones(size(XE,1)/2,1) zeros(size(XE,1)/2,1)]*1e4;
else
  sC = ones(size(XE,1),1)*1e5;
  sT = sC;
  sB = sC*0.1;
end

[A,b,Aeq,beq] = create_constraint_matrices(XX,XE,XC,coms,g,sC,sT,sB,'Mass',masses);
objective = edge_lengths(XX,XE) + wvis*Xvis.^2;
objective(sT>0) = objective(sT>0)*10;
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

%%
clf;
hold on;
ssh = tsurf(SF,SV,'FaceVertexCData',repmat(cbblue,size(SV,1),1),fphong,falpha(1,0),fsoft);
% psh = surf(QX,QY,QZ, ...
%   'CData',PP,'AlphaData',PP,fphong,'EdgeColor','none','FaceAlpha','interp');
% %psh = tsurf(PF,PV,'CData',PP,fphong,falpha(1,0),fsoft);
% colormap(interp1([1 0],[cbred;1 1 1],linspace(0,1,8)));
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
wsh = tsurf(XE(WZ,:),XX,'LineWidth',1,'EdgeColor','k','FaceColor','none');

hold off;
axis equal;
view(-141,-1);

l = { ...
 light('Color',0.5*[1 1 1],'Position',(campos-camtarget),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position',-(campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],pi*0.9),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],-pi*0.9),'Style','local') ...
 light('Color',0.25*[1 1 1],'Position', [-1 -1 10],'Style','local') ...
};
% apply_ambient_occlusion([tsh{:},ssh],'AddLights',false);apply_ambient_occlusion(csh,'AddLights',false,'Factor',0.5);
add_shadow({ssh,tsh{:}},l{6});
%apply_ambient_occlusion([],'AddLights',false);
camproj('persp');
%view(0,0);
set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');

% the loopy animation needs some work :(
% for t = linspace(0,2*pi,60);camtarget(mean(VV));camup([0 0 1]);campos(mean(Q)+0.05*[((1+pi.*cos(t)).*cos(t)-pi+1)/2 0 0.5*(1+pi.*cos(t)).*sin(t)/2]);camproj('persp');camva(50);drawnow;end
% camtarget(mean(VV)-[0.3 0 0]);camup([0 0 1]);campos(mean(R)-[0 .3 0]);camproj('persp');camva(50);
% camtarget(mean(VV)-[0.3 0 0]);camup([0 0 1]);campos(mean(Q)-[0.3 0 0]);camproj('persp');camva(50);

%%
% save('teaser.mat','XX','XE','RZ','WZ','Q','R',...
%       'x','ar','ax','be','fval','sC','sT','sB','Xvis');