addpath ../cyCodeBase/
[V,F] = load_mesh('~/Dropbox/models/cartoon-elephant.obj');
V=V* axisangle2matrix([0 1 1],pi) * axisangle2matrix([0 0 1],-pi*0.3)* ...
  axisangle2matrix([0 1 0],pi*0.1);
% rescale to unit sphere
V = V-0.5*(max(V)+min(V));
V = V/max(normrow(V));
% 0.3 m wide
V = V*0.4*0.5;
V = V + [0.5 0.8 0.5];

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
Q = mean(PV);
%nu = 100;
%%Q = cy_blue_noise(nu*300,PV,PF);
%Q = random_points_on_mesh(PV,PF,nu*30);
%QP = interp2(PX,PZ,PP,Q(:,1),Q(:,3));
%[~,QI] = histc(rand(nu,1),[0;cumsum(QP)]/sum(QP));
%Q = Q(QI,:);

n = 1600;
VV = [SV;V];
FF = [SF;size(SV,1)+F];
CC = [ones(size(SF,1),1);1+ones(size(F,1),1)];
tic;
[XX,XE,XC,YX,YE,YC] = groundstructure(VV,FF,CC,n);
fprintf('groundstructure: %g secs\n',toc);

r = 0.001;
W = ones(size(Q,1),1);
W = W./sum(W);
vis_method = 'exact-sampling';
tic;
Zint = groundstructure_visibility( ...
  Q,VV,FF,XX,XE,'SampleSize',r,'Method',vis_method,'Weight',W);
fprintf('groundstructure_visibility: %g secs\n',toc);
Zint = matrixnormalize(Zint);

coms = [nan nan nan;centroid(V,F)];
force = [0 0 -9.8];
sC = repmat(1e2,size(XE,1),1);
sT = repmat(1e2,size(XE,1),1);
sB = repmat(1e3,size(XE,1),1);
tic;
[A,b,Aeq,beq] = create_constraint_matrices(XX,XE,XC,coms,force,sC,sT,sB,'bending',1);
fprintf('create_constraint_matrices: %g secs\n',toc);
objective = edge_lengths(XX,XE)+100*Zint.^2;
tic;
[x,ar,ax,be] = optimize_lp(objective,A,b,Aeq,beq,'linprog','bending',1); 
fprintf('optimize_lp: %g secs\n',toc);
be = reshape(be,[],2);
NZ = ar>1e-6;
%[ar(NZ) ax(NZ) be(NZ,:)]

[CV,CF] = edge_cylinders(XX,XE(NZ,:),'Thickness',0.1*sqrt(ar(NZ)),'PolySize',30);

clf;
hold on;
ssh = tsurf(SF,SV,'FaceVertexCData',repmat(blue,size(SV,1),1),fphong,falpha(1,0),fsoft);
%red = [1 0 0];
%psh = surf(PX,PY,PZ, ...
%  'CData',PP,'AlphaData',PP,fphong,'EdgeColor','none','FaceAlpha','interp');
%psh = tsurf(PF,PV,'CData',PP,fphong,falpha(1,0),fsoft);
%colormap(interp1([1 0],[cbred;1 1 1],linspace(0,1,7)));
%psh.DiffuseStrength = 0;
%psh.AmbientStrength = 1;
%psh.SpecularStrength = 0;
tsh = {};
CM = interp1([0 1],[cbgreen;cbgreen*0.6+0.4],linspace(0,1,numel(V))');
tsh = tsurf( ...
  F,V,'FaceVertexCData',repmat(CM(ii,:),size(V,1),1), ...
  fphong,falpha(1,0),fsoft);
%sct(Q,'.k');

csh = tsurf(CF,CV, ...
  falpha(1,0),fsoft,'FaceVertexCData',repmat(cborange,size(CV,1),1));

hold off;
axis equal;
view(-70,8);

l = { ...
 light('Color',0.5*[1 1 1],'Position',(campos-camtarget),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position',-(campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],pi*0.9),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],-pi*0.9),'Style','local') ...
 light('Color',0.25*[1 1 1],'Position', [-1 -1 10],'Style','local') ...
};
AO = apply_ambient_occlusion([ssh,tsh,csh],'AddLights',false);
add_shadow([ssh,tsh,csh],l{6});
camproj('persp');
%view(0,0);
%set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');

title(sprintf('obj: %g', objective'*ar),'FontSize',30);
