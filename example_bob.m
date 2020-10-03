rng(0);
addpath ../cyCodeBase/
addpath utils/
[V,F] = load_mesh('data/meshes/bob/bob_tri.obj');
V = V*axisangle2matrix([1 0 0],-pi/2)*axisangle2matrix([0 1 0],-0.05)+[0 0 1.2];
coms = [nan nan nan;centroid(V,F).*[1 0 1]];

[SV,SF] = create_regular_grid(6,3);
%BC = barycenter(SV,SF);
%bi = snap_points([0 0;0.2 0;1 0;0.2 1;1 1;0 1],BC);
%SF = SF(bi,:);
%SF = fliplr(SF);
SV = 0.5*(SV-0.5)*[1 0 0;0 0.5 0]+[0.2 0 0];
XBC = barycenter(SV,SF);

VV = [SV;V];
FF = [SF;size(SV,1)+F];
CC = [ones(size(SF,1),1);1+ones(size(F,1),1)];

N = normalizerow(normals(V,F));
n = 1000;
[X,I,B] = blue_noise(n,V,F,'Seed',rand);
XN = [N(I,:);nan(size(XBC,1),3)];
XX = [X;XBC];
XC = [repmat(2,size(X,1),1);ones(size(XBC,1),1)];
nc = 2;
AC = ones(nc,nc)-eye(nc);
X2C = sparse(1:size(XX,1),XC,1,size(XX,1),nc);
% sample adjacency
A = X2C*AC*X2C';
[AI,AJ] = find(triu(A));
% edge
XE = [AI,AJ];
% edge vector
XEV = XX(XE(:,2),:)-XX(XE(:,1),:);
XEU = normalizerow(XEV);
% unpruned groundstructure
YX = XX;
YE = XE;
YC = XC;
% discard if the angle w.r.t. to endpoint normal is too large
max_angle = 80/180*pi;
valid = ...
  sum(XEU.*XN(XE(:,1),:),2) > cos(max_angle);
XE = XE(valid,:);
XEV = XEV(valid,:);
XEU = XEU(valid,:);

% discard if intersecting something else in the scene
[H,T] = ray_mesh_intersect(XX(XE(:,1),:),XEV,V,F);
H = H & T<0.999999;
XE =   XE(~H,:);
XEV = XEV(~H,:);
XEU = XEU(~H,:);


% symmetry
sp = [0 0 0];
sn = [0 1 0];
YE = [YE;size(XX,1)+YE];
YX = [YX;YX-2*sum((YX-sp).*sn,2).*sn];
force = [0 0 9.8];
sT = ones(size(YE,1),1)*1e3;
sC = ones(size(YE,1),1)*0;
sB = ones(size(YE,1),1)*0;

Xvis = groundstructure_visibility(Q,VV,FF,YX,YE,'SampleSize',r);

%As = [];
As = [speye(size(YE,1)/2) -speye(size(YE,1)/2) sparse(size(YE,1)/2,size(YE,1)*3)];
bs = zeros(size(As,1),1);
[A,b,Aeq,beq] = create_constraint_matrices(YX,YE,YC,coms,force,sC,sT,sB);
objective = edge_lengths(YX,YE);
[x,ar,ax,be,fval] = optimize_lp(objective,A,b,[Aeq;As],[beq;bs],'linprog');
NZ = find(max(ar,0)>1e-7);
[ar(NZ) ax(NZ) be(NZ,:)]


[HV,HF] = load_mesh('data/meshes/bob/Model3-pose-0000.obj');
HV = (HV*axisangle2matrix([1 0 0],-pi/2)-[0.3822 -0.7546 0.778])*0.1*axisangle2matrix([0 0 1],pi/2);
[HV,HF] = repmesh(HV,HF,YX(YE(NZ,2),:));

[RV,RF] = create_regular_grid(6,3);
RV = 0.5*(RV-0.5)*[5 0 0;0 0.5 0]+[0.2 0 min(HV(:,3))];

%[hV,hF,UV,TF,N,NF] = readOBJ('data/meshes/bob/bob-quads-high.obj','Quads',true);
%hV = hV*axisangle2matrix([1 0 0],-pi/2)*axisangle2matrix([0 1 0],-0.05)+[0 0 1.2];
%tex = im2double(imread('data/meshes/bob/bob_diffuse-green.png'));
%[X,Y] = meshgrid(linspace(0,1,size(tex,2)),linspace(0,1,size(tex,1)));
%C = [];
%C(:,1) = interp2(X,Y,flipud(tex(:,:,1)),UV(TF,1),UV(TF,2));
%C(:,2) = interp2(X,Y,flipud(tex(:,:,2)),UV(TF,1),UV(TF,2));
%C(:,3) = interp2(X,Y,flipud(tex(:,:,3)),UV(TF,1),UV(TF,2));
%CV = full(sparse(repmat(hF(:),1,3),repmat([1 2 3],numel(hF),1),C,size(hV,1),3)./ ...
%  sparse(hF,1,1,size(hV,1),1));
%%CU = full(sparse(repmat(TF(:),1,3),repmat([1 2 3],numel(TF),1),C,size(UV,1),3)./ ...
%%  sparse(TF,1,1,size(UV,1),1));
%%[dV,dF] = decimate_libigl([V CV],[F(:,[1 2 3]);F(:,[1 3 4])],20000,'Method','qslim');
%%dC = dV(:,4:6);
%%dV = dV(:,1:3);
%%save('data/meshes/bob/bob.mat','V','F','TF','UV','C','CV','dV','dF','dC');


clf;
hold on;
%ssh = tsurf(RF,RV,'FaceVertexCData',repmat(cbblue*0.05+0.95,size(SF,1),1),falpha(1,0),fsoft);
hsh = tsurf(HF,HV,'FaceVertexCData',repmat(cbblue,size(HF,1),1),falpha(1,0),fsoft);
tsh = tsurf(hF,hV,'Tets',0,'FaceVertexCData',0.5*CV+0.5*cbgreen,falpha(1,0),fphong,fsoft);
tsh.SpecularExponent = 100;
tsh.SpecularStrength = 0.4
wsh = tsurf(XE(NZ,:),XX,'LineWidth',1,'EdgeColor','k','FaceColor','none');
hold off;
axis equal;
camproj('persp');
view(-46,-2);
l = { ...
 light('Color',0.5*[1 1 1],'Position',(campos-camtarget),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position',-(campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],pi*0.9),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],-pi*0.9),'Style','local') ...
 light('Color',0.25*[1 1 1],'Position', [-1 -1 10],'Style','local') ...
};
AO = apply_ambient_occlusion(tsh,'AddLights',false,'AO',AO);
add_shadow({hsh,tsh},l{6});
set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');
figpng('bob-raw.png');
