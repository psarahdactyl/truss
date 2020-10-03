rng(0);

addpath ../cyCodeBase/
addpath utils/
[AV,AF] = load_mesh('data/meshes/pterosaure-light-solid.obj');
AV = AV-0.5*(max(AV)+min(AV));
AV = AV/max(normrow(AV));

[SV,SF] = create_regular_grid(2,2);
SV = 2*(SV*[1 0 0;0 1 0]+[-0.5 -0.5 0])+[0 0 1];
SF=fliplr(SF);
[KV,KF] = load_mesh('data/meshes/panel-decimated.ply');
KV = KV*axisangle2matrix([1 0 0],pi);
KV(:,1:2) = KV(:,1:2) - 0.5*(max(KV(:,1:2))+min(KV(:,1:2)));
KV(:,3) = KV(:,3)- 0.5*(max(KV(:,3))+min(KV(:,3)));
KV = KV/max(abs(KV(:)));
KV = KV+[0 0 1];

[~,AC] = connected_components(AF);
% messed up pieces to remove
bad = [12 13 43 44 45 46 65]-1;
[AV,~,~,AF] = remove_unreferenced(AV,AF(~ismember(AC,bad),:));
[~,AC] = connected_components(AF);

V = {};
F = {};
coms = [nan nan nan];
%for ii = 1:max(AC)

for ii = 1:max(AC)
  [Vi,~,~,Fi] = remove_unreferenced(AV,AF(AC==ii,:));
  V{ii} = Vi;
  F{ii} = Fi;
  coms = [coms;centroid(Vi,Fi)];
end



VV = cell2mat(reshape(V,[],1));
nV = cumsum([0 cellfun(@(V) size(V,1),{V{:}})]);
FF = cell2mat(arrayfun(@(i) nV(i)+F{i},1:numel(V),'UniformOutput',false)');
CC = cell2mat(arrayfun(@(i) repmat(i,size(F{i},1),1),1:numel(F),'UniformOutput',false)');
VV = [SV;VV];
FF = [SF;size(SV,1)+FF];
CC = [ones(size(SF,1),1);1+CC];
nc = max(CC);
areas = cell2mat(arrayfun(@(ci) ...
  sum(doublearea(VV,FF(CC==ci,:))),(1:nc)','UniformOutput',false));
% n = 40 looked pretty good
n = 80;
ns = 20;
nn = repmat(n,max(CC),1);
nn(1) = ns;
[XX,XE,XC,YX,YE,YC] = groundstructure(VV,FF,CC,nn);
count = accumarray(XC(XE(:)),1);
assert(min(count>6));


r = 0.1;
Xvis = groundstructure_visibility(Q,VV,FF,XX,XE,'SampleSize',r);
Xvis = zeros(size(XE,1),1);
force = [0 0 -9.8];
XE = [XE;XE];
Xvis = [Xvis;Xvis];
sT = [ones(size(XE,1)/2,1)*1e5;ones(size(XE,1)/2,1)*0.0];
sC = [ones(size(XE,1)/2,1)*0.0;ones(size(XE,1)/2,1)*1e4];
sB = ones(size(XE,1),1)*1e3;

l = edge_lengths(XX,XE);
bad = l>0.12 & sC>0;
XE(bad,:) = [];
Xvis(bad,:) = [];
sC(bad,:) = [];
sT(bad,:) = [];
sB(bad,:) = [];


%% symmetry
%sp = [10 0 0];
%sn = [1 0 0];
%XE = [XE;size(XX,1)+XE];
%XX = [XX;XX-2*sum((XX-sp).*sn,2).*sn];
%XC = repmat(XC,2,1);
%sC = repmat(sC,2,1);
%sT = repmat(sT,2,1);
%sB = repmat(sB,2,1);
%%Xvis = repmat(Xvis,2,1);
%
wvis = 0;
[A,b,Aeq,beq] = create_constraint_matrices(XX,XE,XC,coms,force,sC,sT,sB);
l = edge_lengths(XX,XE);
objective = l + wvis*Xvis.^2;
objective(sT>0) = objective(sT>0)*1000;
[x,ar,ax,be,fval] = optimize_lp(objective,A,b,Aeq,beq,'linprog');

NZ = find(max(ar,0)>1e-7);
%[ar(NZ) ax(NZ) be(NZ,:)]
RZ = ar>1e-6 & ax>0;
WZ = ar>1e-6 & ax<0;
[CV,CF] = edge_cylinders(XX,XE(RZ,:),'Thickness',0.15*sqrt(ar(RZ)),'PolySize',30);

CB = cbrewer('Set1',5);
cbred = CB(1,:);
cbblue = CB(2,:);
cbgreen = CB(3,:);
cborange = CB(5,:);

CM = interp1([0 1],[cbgreen;cbgreen*0.8+0.2],linspace(0,1,numel(V))');

clf;
hold on;
%ssh = tsurf(SF,SV,'FaceVertexCData',repmat(cbblue,size(SF,1),1),falpha(1,0),fsoft);
ssh = tsurf(KF,KV,'FaceVertexCData',repmat(cbblue,size(KF,1),1),falpha(1,0),fsoft);
tsh = {};
for ii = 1:numel(V)
  tsh{end+1} = tsurf( ...
    F{ii},V{ii},'FaceVertexCData',repmat(CM(ii,:),size(F{ii},1),1), ...
    falpha(1.0,0),fsoft);
end
%wsh = tsurf(XE(NZ,:),XX,'LineWidth',1,'EdgeColor','k','FaceColor','none');
csh = tsurf(CF,CV, ...
  falpha(1,0),fsoft,'FaceVertexCData',repmat(cborange,size(CF,1),1));
wsh = tsurf(XE(WZ,:),XX,'LineWidth',0.1,'EdgeColor','k',falpha(0,0.2),'FaceColor','none');
hold off;
axis equal;
camproj('persp');
view(15,-21);
l = { ...
 light('Color',0.5*[1 1 1],'Position',(campos-camtarget),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position',-(campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],pi*0.9),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],-pi*0.9),'Style','local') ...
 light('Color',0.25*[1 1 1],'Position', [-1 10 500],'Style','local') ...
};
AOt = apply_ambient_occlusion([tsh{:}],'AddLights',false,'AO',AOt);
AOs = apply_ambient_occlusion([ssh],'AddLights',false,'AO',AOs,'Factor',0.5);
camproj('persp');
title(sprintf('%d:%d',sum(RZ),sum(WZ)),'FontSize',40);
set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');
