addpath ../cyCodeBase/
addpath utils/
rng(0);
[SV,SF] = create_regular_grid(2,2);
SV = 0.25*(SV*[1 0 0;0 1 0]+[-0.5 -0.5 0])+[0 0 2];
SF=fliplr(SF);

V = {};
F = {};
[V{end+1},F{end+1}] = load_mesh('data/meshes/cupped-sgp/cupped-sgp-left.obj');
[V{end+1},F{end+1}] = load_mesh('data/meshes/cupped-sgp/cupped-sgp-S.obj');
[V{end+1},F{end+1}] = load_mesh('data/meshes/cupped-sgp/cupped-sgp-G.obj');
[V{end+1},F{end+1}] = load_mesh('data/meshes/cupped-sgp/cupped-sgp-P.obj');
[V{end+1},F{end+1}] = load_mesh('data/meshes/cupped-sgp/cupped-sgp-right.obj');
F{end} = fliplr(F{end});
T = @(V) V*axisangle2matrix([1 0 0],-pi/2)*0.1+[0 0 1];
V = cellfun(@(V) T(V),V,'UniformOutput',false);
coms = [nan nan nan;cell2mat(arrayfun(@(i) centroid(V{i},F{i}),1:numel(V),'UniformOutput',0)')];


CB = cbrewer('Set1',5);
cbred = CB(1,:);
cbblue = CB(2,:);
cbgreen = CB(3,:);
cborange = CB(5,:);

VV = cell2mat(reshape(V,[],1));
nV = cumsum([0 cellfun(@(V) size(V,1),{V{:}})]);
FF = cell2mat(arrayfun(@(i) nV(i)+F{i},1:numel(V),'UniformOutput',false)');
CC = cell2mat(arrayfun(@(i) repmat(i,size(F{i},1),1),1:numel(F),'UniformOutput',false)');
VV = [SV;VV];
FF = [SF;size(SV,1)+FF];
CC = [ones(size(SF,1),1);1+CC];
n = 100;
[XX,XE,XC,YX,YE,YC] = groundstructure(VV,FF,CC,n);
XE(XC(XE(:,1))==1&XC(XE(:,2))>2,:) = [];
force = [0 0 -9.8];
sC = 0*ones(size(XE,1),1)*1e2;
sT = ones(size(XE,1),1)*1e2;
sB = 0*ones(size(XE,1),1)*1e3;

[A,b,Aeq,beq] = create_constraint_matrices(XX,XE,XC,coms,force,sC,sT,sB);
objective = edge_lengths(XX,XE);
[x,ar,ax,be,fval] = optimize_lp(objective,A,b,Aeq,beq,'linprog');

NZ = ar>1e-6;
sum(NZ)

%[CV,CF] = edge_cylinders(XX,XE(NZ,:),'Thickness',0.015*sqrt(ar(NZ)),'PolySize',30);

%

ssh = [];
clf;
hold on;
ssh = tsurf(SF,SV,'FaceVertexCData',repmat(blue,size(SV,1),1),fphong,falpha(1,0),fsoft);

tsh = {};
CM = interp1([0 1],[cbgreen;cbgreen*0.8+0.2],linspace(0,1,numel(V))');
for ii = 1:numel(V)
  tsh{end+1} = tsurf( ...
    F{ii},V{ii},'FaceVertexCData',repmat(CM(ii,:),size(F{ii},1),1), ...
    fphong,falpha(1.0,0),fsoft);
end
%csh = tsurf(CF,CV, ...
%  falpha(1,0),fsoft,'FaceVertexCData',repmat(cborange,size(CF,1),1));
wsh = tsurf(XE(NZ,:),XX,'LineWidth',1,'EdgeColor','k','FaceColor','none');
%sct(XX,'.y','SizeData',500);
%sct(coms(2,:),'.r','SizeData',1000);
hold off;
axis equal;
view(90,0);

l = { ...
 light('Color',0.5*[1 1 1],'Position',(campos-camtarget),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position',-(campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],pi*0.9),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],-pi*0.9),'Style','local') ...
 light('Color',0.25*[1 1 1],'Position', [-1 -1 10],'Style','local') ...
};
%apply_ambient_occlusion([tsh{:},ssh],'AddLights',false);apply_ambient_occlusion(csh,'AddLights',false,'Factor',0.5);
apply_ambient_occlusion(tsh,'AddLights',false,'Factor',1);
add_shadow({tsh{:}},l{6},'Nudge',0.1);
camproj('persp');
axis vis3d
set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');

ts = linspace(0,360); ts = ts(1:end-1); for t = ts view(90+t,0); drawnow; end
