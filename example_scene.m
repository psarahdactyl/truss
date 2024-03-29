rng(6)
V = {};
F = {};
coms = [nan nan nan];
[V{end+1},F{end+1}] = load_mesh('data/meshes/bunny-remesh.obj');
V{end} = V{end} *  ...
  axisangle2matrix([0 0 1],pi*0.1);
coms(end+1,:) = centroid(V{end},F{end});
[V{end+1},F{end+1}] = load_mesh('data/meshes/rocker-arm.obj');
V{end} = V{end} *  ...
  axisangle2matrix([1 1 0],-pi*0.4) * ...
  axisangle2matrix([0 0 1],-pi*0.4) * ...
  axisangle2matrix([1 0 0],-pi*0.1);
coms(end+1,:) = centroid(V{end},F{end});
[V{end+1},F{end+1}] = load_mesh('data/meshes/teapot.obj');
V{end} = V{end} *  ...
  axisangle2matrix([1 0 0],-pi*0.4) * axisangle2matrix([0 0 1],pi/5) * ...
    axisangle2matrix([0 1 0],pi/10);
coms(end+1,:) = shell_centroid(V{end},F{end});
[V{end},F{end}] = selfintersect(V{end},F{end},'StitchAll',true);
for ii = 1:numel(V)
  % rescale to unit sphere
  V{ii} = V{ii}-0.5*(max(V{ii})+min(V{ii}));
  V{ii} = V{ii}/max(normrow(V{ii}));
  % 0.3 m wide
  V{ii} = V{ii}*0.3*0.5;
end
V{1} = V{1} + [0.6 0.5 0.8];
V{2} = V{2} + [0.53 0.1 0.6];
V{3} = V{3} + [0.6 0.6 0.5];

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
PP = clamp(max(PP)-PP,0.2,1);
PP = matrixnormalize(PP);
PP = reshape(PP,size(X));
PX = X*0.75+0.25*0.5;
PY = 0*Y-0.75;
PZ = Y*0.5+0.25;
[PF,PV] = surf2patch(PX,PY,PZ,'triangles');
nu = 200;
%Q = cy_blue_noise(nu*300,PV,PF);
Q = random_points_on_mesh(PV,PF,nu*30);
QP = interp2(PX,PZ,PP,Q(:,1),Q(:,3));
[~,QI] = histc(rand(nu,1),[0;cumsum(QP)]/sum(QP));
Q = Q(QI,:);

n = 1000;

VV = cell2mat(reshape(V,[],1));
nV = cumsum([0 cellfun(@(V) size(V,1),{V{:}})]);
FF = cell2mat(arrayfun(@(i) nV(i)+F{i},1:numel(V),'UniformOutput',false)');
CC = cell2mat(arrayfun(@(i) repmat(i,size(F{i},1),1),1:numel(F),'UniformOutput',false)');
VV = [SV;VV];
FF = [SF;size(SV,1)+FF];
CC = [ones(size(SF,1),1);1+CC];
[XX,XE,XC,YX,YE,YC] = groundstructure(VV,FF,CC,n);
% XE(XC(XE(:,1))==1&XC(XE(:,2))==3,:) = [];
YE(YC(YE(:,1))==1&YC(YE(:,2))==3,:) = [];
r = 0.001;
Xvis = groundstructure_visibility(Q,VV,FF,YX,YE,'SampleSize',r);

force = [0 0 -9.8];
sC = ones(size(YE,1),1)*1e3;
sT = ones(size(YE,1),1)*1e3;
sB = ones(size(YE,1),1)*1e3;

[A,b,Aeq,beq] = create_constraint_matrices(YX,YE,YC,coms,force,sC,sT,sB);
objective = edge_lengths(YX,YE) + 10000*Xvis.^2;
[x,ar,ax,be,fval] = optimize_lp(objective,A,b,Aeq,beq,'linprog');

NZ = ar>1e-6;
[CV,CF] = edge_cylinders(YX,YE(NZ,:),'Thickness',0.015*sqrt(ar(NZ)),'PolySize',30);

%%
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
dsh = add_shadow({csh,ssh,tsh{:}},l{6});
%apply_ambient_occlusion([],'AddLights',false);
camproj('persp');
%view(0,0);
set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');

%for t = linspace(0,2*pi,60);camtarget([0.5 1 0.5]);camup([0 0 1]);campos(mean(views)+0.05*[((1+pi.*cos(t)).*cos(t)-pi+1)/2 0 0.5*(1+pi.*cos(t)).*sin(t)/2]);camproj('persp');camva(50);drawnow;end

% camtarget([0.5 1 0.5]); camup([0 0 1]); campos(mean(PV));camva(50);
% vidobj = [];
% vidname = 'video-footage/bunny-teapot-rocker-jiggle.mp4';
% vidobj = figmp4(vidname,vidobj,10);
% jiggle(0.1,120,@() figmp4(vidname,vidobj) )
% vidobj = figmp4(vidname,vidobj,10);
% vidobj.close();
