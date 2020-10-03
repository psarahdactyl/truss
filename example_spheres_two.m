addpath ../cyCodeBase
addpath utils

rng(0)
V = {};
F = {};
coms = [nan nan nan];
% coms = [];

vols = nan;
masses = nan;

hold on

% XX=zeros(4,3);
XX(1,:) = [0.5 1 0.58];
% XX(1,:) = [0.5 1 0.5];
XX(2,:) = [0.5 1 0.42];

[V{end+1},F{end+1}] = subdivided_sphere(3,'Radius',0.08);
V{end} = V{end} + [0.5 0.7 0.5];
[coms(end+1,:),vols(end+1)] = centroid(V{end},F{end});
[M,I]=max(V{end}(:,3));
XX(3,:) = V{end}(I,:);
[M,I]=min(V{end}(:,3));
XX(4,:) = V{end}(I,:);

% XE = [
%       1 3;
%       1 4;
%       2 3;
%       2 4;
%       ];
%  XC = [1 1 2 2];

% [M,I]=min(V{end}(:,2));
% XX(2,:) = V{end}(I,:);
% XE = [1 2];
%     XC=[1 2];

[V{end+1},F{end+1}] = subdivided_sphere(3,'Radius',0.08);
V{end} = V{end} + [0.5 0.4 0.5];
[coms(end+1,:),vols(end+1)] = centroid(V{end},F{end});
[M,I]=max(V{end}(:,3));
XX(5,:) = V{end}(I,:);
[M,I]=min(V{end}(:,3));
XX(6,:) = V{end}(I,:);

CB = cbrewer('Set1',5);
cbred = CB(1,:);
cbblue = CB(2,:);
cbgreen = CB(3,:);
cborange = CB(5,:);

[SV,SF] = cube(2,2);
SV = SV.*[1 0.05 1] + [0 1 0];

% units are meters
[X,Y] = meshgrid(linspace(0,1),linspace(0,1));
PP = normrow([X(:) Y(:)]-mean([X(:) Y(:)]));
PX = X*0.75+0.25*0.5;
PY = 0*Y-0.75;
PZ = Y*0.5+0.25;
PP = clamp(max(PP)-PP,0.2,1);
PP = matrixnormalize(PP);
PP = reshape(PP,size(X));
[PF,PV] = surf2patch(PX,PY,PZ,'triangles');
nu = 100;
%Q = blue_noise(nu*300,PV,PF);
Q = random_points_on_mesh(PV,PF,nu*30);
QP = interp2(PX,PZ,PP,Q(:,1),Q(:,3));
[~,QI] = histc(rand(nu,1),[0;cumsum(QP)]/sum(QP));
Q = Q(QI,:);

n = 400;
VV = cell2mat(reshape(V,[],1));
nV = cumsum([0 cellfun(@(V) size(V,1),{V{:}})]);
FF = cell2mat(arrayfun(@(i) nV(i)+F{i},1:numel(V),'UniformOutput',false)');
CC = cell2mat(arrayfun(@(i) repmat(i,size(F{i},1),1),1:numel(F),'UniformOutput',false)');
VV = [SV;VV];
FF = [SF;size(SV,1)+FF];
CC = [ones(size(SF,1),1);1+CC];
% [XX,XE,XC,YX,YE,YC] = groundstructure(VV,FF,CC,n);
XE = [
%       1 3;
      1 4;
      2 3;
%       2 4;
%       3 5;
      3 6;
      4 5;
%       4 6
      ];

% XE = [
%       1 3;
%       1 4;
%       2 3;
%       2 4;
%       3 5;
%       3 6;
%       4 5;
%       4 6
%       ];
XC = [1 1 2 2 3 3];

r = 0.1;
W = ones(size(Q,1),1);
W = W./sum(W);
vis_method = 'exact-sampling';
Zint = groundstructure_visibility( ...
  Q,VV,FF,XX,XE,'SampleSize',r,'Method',vis_method,'Weight',W);
Zint = matrixnormalize(Zint);

force = [0 0 -9.8];
sC = repmat(1e3,size(XE,1),1);
sT = repmat(1e3,size(XE,1),1);
sB = repmat(1e3,size(XE,1),1);
[A,b,Aeq,beq] = create_constraint_matrices(XX,XE,XC,coms,force,sC,sT,sB);
objective = edge_lengths(XX,XE);%+10000*Zint.^2;
[x,ar,ax,be] = optimize_lp(objective,A,b,Aeq,beq,'linprog'); 
NZ = ar>1e-6;

[CV,CF] = edge_cylinders(XX,XE(NZ,:),'Thickness',0.2*sqrt(ar(NZ)),'PolySize',30);

EVV = [VV;CV];
EFF = [FF;size(VV,1)+CF];
writeOBJ('spheres-example.obj',EVV,EFF);

for pass = 1

clf;
hold on;
ssh = tsurf(SF,SV,'FaceVertexCData',repmat(cbblue,size(SV,1),1),fphong,falpha(1,0),fsoft);
% red = [1 0 0];
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
%sct(Q,'.k');
csh = [];
switch pass
case 4
  tsurf(YE,YX,'FaceColor','k',falpha(0,0.5));
case 3
  tsurf(reshape(1:numel(XE),size(XE)),XX(XE,:),'CData',1*[Zint;Zint], ...
    falpha(0,0.5),'EdgeColor','flat')
case {1,2}
%   csh = tsurf(CF,CV, ...
%     falpha(1,0),fsoft,'FaceVertexCData',repmat(cborange,size(CV,1),1));
end

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
AO = apply_ambient_occlusion([ssh,tsh{:},csh],'AddLights',false);
add_shadow({ssh,tsh{:}},l{6});
camproj('persp');
%view(0,0);
set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');

% if pass ==2
%   camtarget([0.5 1 0.5]);camup([0 0 1]);campos(mean(PV));camproj('persp');camva(50);
% end
% figpng(sprintf('pruned-bunny-%d.png',pass));

end
