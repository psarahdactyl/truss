rng(0);
addpath ../cyCodeBase/
addpath utils/
[SV,SF] = load_mesh('data/meshes/hermitage/hermitage.obj');
[uV,uF] = upsample(SV,SF,'Iterations',10,'OnlySelected',@(V,F) find(doublearea(V,F)>1));
%[V,F] = load_mesh('~/Dropbox/models/balloon-dog.obj');
[V,F] = load_mesh('data/meshes/hermitage/balloon-dog-decimated.obj');
V = V*axisangle2matrix([0 0 1],pi/2);
V = V*axisangle2matrix([1 0 0],-0.1*pi);

V = V*0.02+[-35 -65 20];

CB = cbrewer('Set1',5);
cbred = CB(1,:);
cbblue = CB(2,:);
cbgreen = CB(3,:);
cborange = CB(5,:);

[X,Y] = meshgrid(linspace(-1,1),linspace(-1,1));
PP = normrow([X(:) Y(:)]-mean([X(:) Y(:)]));
PX = 10*X-35;
PY = 25*Y-80;
PZ = 0*X+1;
PP = clamp(max(PP)-PP,0.2,1);
PP = matrixnormalize(PP);
PP = reshape(PP,size(X));
[PF,PV] = surf2patch(PX,PY,PZ,'triangles');
nu = 100;
%Q = cy_blue_noise(nu*300,PV,PF);
Q = random_points_on_mesh(PV,PF,nu*30);
QP = interp2(PX,PY,PP,Q(:,1),Q(:,2));
[~,QI] = histc(rand(nu,1),[0;cumsum(QP)]/sum(QP));
Q = Q(QI,:);

n = 800;
VV = [SV;V];
BF = SF(((barycenter(SV,SF)*[0;1;0])<-38) & ((barycenter(SV,SF)*[0;0;1])>10),:);
FF = [BF;size(SV,1)+F];
CC = [ones(size(BF,1),1);1+ones(size(F,1),1)];
[XX,XE,XC,YX,YE,YC] = groundstructure(VV,FF,CC,n);


r = 0.1;
W = ones(size(Q,1),1);
W = W./sum(W);
vis_method = 'exact-sampling';
Xvis = groundstructure_visibility( ...
  Q,VV,FF,YX,YE,'SampleSize',r,'Method',vis_method,'Weight',W);


force = [0 0 -9.8];
wires = true;
wvis = 10000;
%wvis = 0;
if wires
  YE = [YE;YE];
  Xvis = [0*Xvis;Xvis];
  sT = [ones(size(YE,1)/2,1);zeros(size(YE,1)/2,1)]*1e3;
  sC = [zeros(size(YE,1)/2,1);ones(size(YE,1)/2,1)]*1e1;
else
  sC = ones(size(YE,1),1)*1e1;
  sT = sC;
end
sB = ones(size(YE,1),1)*0;

% symmetry
sp = [-35 0 0];
sn = [1 0 0];
YE = [YE;size(YX,1)+YE];
YX = [YX;YX-2*sum((YX-sp).*sn,2).*sn];
YC = repmat(YC,2,1);
sC = repmat(sC,2,1);
sT = repmat(sT,2,1);
sB = repmat(sB,2,1);
Xvis = repmat(Xvis,2,1);

coms = [nan nan nan;centroid(V,F)];
[A,b,Aeq,beq] = create_constraint_matrices(YX,YE,YC,coms,force,sC,sT,sB);
objective = edge_lengths(YX,YE) + wvis*Xvis.^2;
% [x,ar,ax,be,fval] = optimize_lp(objective,A,b,[Aeq;As],[beq;bs],'linprog');
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
[CV,CF] = edge_cylinders(YX,YE(RZ,:),'Thickness',1.*sqrt(ar(RZ)),'PolySize',30);


clf;
hold on;
%ssh = tsurf(SF,SV,'FaceVertexCData',repmat(blue,size(SV,1),1),fphong,falpha(1,0),fsoft);
ssh = tsurf(uF,uV,'FaceVertexCData',repmat(blue,size(uV,1),1),fphong,falpha(1,0),fsoft);
red = [1 0 0];
psh = surf(PX,PY,PZ, ...
  'CData',PP,'AlphaData',PP,fphong,'EdgeColor','none','FaceAlpha','interp');
%psh = tsurf(PF,PV,'CData',PP,fphong,falpha(1,0),fsoft);
colormap(interp1([1 0],[cbred;1 1 1],linspace(0,1,7)));
caxis([0 max(PP(:))]);
psh.DiffuseStrength = 0;
psh.AmbientStrength = 1;
psh.SpecularStrength = 0;
tsh = {};
tsh = tsurf( ...
  F,V,'FaceVertexCData',repmat(cbgreen,size(V,1),1), ...
  fphong,falpha(1,0));
tsh.SpecularExponent = 100;
tsh.SpecularStrength = 1;
csh = tsurf(CF,CV, ...
  falpha(1,0),fsoft,'FaceVertexCData',repmat(cborange,size(CF,1),1));
wsh = tsurf(XE(WZ,:),XX,'LineWidth',1,'EdgeColor','k','FaceColor','none');
%sct(XX,'.k');

hold off;
axis equal;
view(-58,15);

l = { ...
 light('Color',0.5*[1 1 1],'Position',(campos-camtarget),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position',-(campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],pi*0.9),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],-pi*0.9),'Style','local') ...
 light('Color',0.25*[1 1 1],'Position', [-1 10 500],'Style','local') ...
};
%AO = apply_ambient_occlusion([ssh,tsh],'AddLights',false);
%apply_ambient_occlusion(tsh,'AddLights',false);
%apply_ambient_occlusion(ssh,'AddLights',false,'Factor',0.5);
add_shadow([ssh,tsh,csh],l{6});
camproj('persp');
title(sprintf('%d out of %d, $%g',numel(NZ),size(E,1),fval),'FontSize',30);
%view(0,0);
%set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');
%set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');

%for t = -50;camtarget(coms(2,:)+[0 50 0]);camup([0 0 1]);campos(mean(PV)+[0 t 0]);camproj('persp');camva(50);drawnow;end
