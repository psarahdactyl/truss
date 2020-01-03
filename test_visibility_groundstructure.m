addpath ../cyCodeBase/
[V,F] = load_mesh('~/Dropbox/models/decimated-knight.off');
rng(0);
[VV,FF,~,CC] = repmesh(V,F,rand(3,3));
[WV,WF] = create_regular_grid(2,2);
WV(:,3) = 0;
WV = WV*2;
FF = [FF;size(VV,1)+WF];
CC = [CC;repmat(max(CC)+1,size(WF,1),1)];
VV = [VV;WV];
R = axisangle2matrix([1 0 0],-pi/2);
VV = VV*R;

% loop over each component gather samples+normals
XX = [];
XN = [];
XC = [];
n = 100;
for ci = 1:max(CC)
  [Vc,~,~,Fc] = remove_unreferenced(VV,FF(CC==ci,:));
  Nc = normalizerow(normals(Vc,Fc));
  [Xc,Ic,Bc] = blue_noise(n*sum(doublearea(Vc,Fc)),Vc,Fc,'Seed',rand);
  Nc = Nc(Ic,:);
  XX = [XX;Xc];
  XN = [XN;Nc];
  XC = [XC;repmat(ci,size(Xc,1),1)];
end

% component adjacency
AC = ones(max(CC),max(CC))-eye(max(CC));
X2C = sparse(1:size(XX,1),XC,1,size(XX,1),max(CC));
% sample adjacency
A = X2C*AC*X2C';
[AI,AJ] = find(triu(A));
% edge
XE = [AI,AJ];
% edge vector
XEV = XX(XE(:,2),:)-XX(XE(:,1),:);
XEU = normalizerow(XEV);

max_angle = 80/180*pi;
valid = ...
  sum(XEU.*XN(XE(:,1),:),2) > cos(max_angle) & ...
  sum(-XEU.*XN(XE(:,2),:),2) > cos(max_angle);
XE = XE(valid,:);
XEV = XEV(valid,:);
XEU = XEU(valid,:);

[H,T] = ray_mesh_intersect(XX(XE(:,1),:),XEV,VV,FF);
H = H & T<0.999999;
XE =   XE(~H,:);
XEV = XEV(~H,:);
XEU = XEU(~H,:);

P = [1 -3 1];
r = 0.1;
[SX,SE,SI] = upsample(XX,XE, ...
  'Iterations',5,'OnlySelected',@(V,E) normrow(V(E(:,2),:)-V(E(:,1),:))>r);
SBC = barycenter(SX,SE);

for iter = 1:2
  switch iter
  case 1
    vis_method = 'approx-sampling';
  case 2
    vis_method = 'exact-sampling';
  end
  tic;
  switch vis_method
  case 'approx-sampling'
    side_x = 10*(max(VV(:,1))-min(VV(:,1)))/r;
    [GX,side] = voxel_grid(VV,side_x,'Pad',1);
    % precompute Z on grid
    [H,T] = ray_mesh_intersect(repmat(P,size(GX,1),1),GX-P,VV,FF);
    H = H & T<0.999999;
    GZ = 1*~H;
    GX = reshape(GX,[side([2 3 1]) 3]);
    Zs = interp3( ...
      GX(:,:,:,1),GX(:,:,:,2),GX(:,:,:,3),reshape(GZ,size(GX(:,:,:,1))), ...
      SBC(:,1),SBC(:,2),SBC(:,3));
  case 'exact-sampling'
    [H,T] = ray_mesh_intersect(repmat(P,size(SBC,1),1),SBC-P,VV,FF);
    H = H & T<0.999999;
    Zs = 1*~H;
  end
  fprintf('%s: %g secs\n',vis_method,toc);
  % integrated viZibility  
  Zavg = sparse(SI,1,Zs,size(XE,1),1)./sparse(SI,1,1,size(XE,1),1);
  Xl = normrow(XX(XE(:,2),:)-XX(XE(:,1),:));
  Zint = Xl.*Zavg;
  switch vis_method
  case 'approx-sampling';
    Zsa = Zs;
    Zinta = Zint;
  case 'exact-sampling';
    Zse = Zs;
    Zinte = Zint;
  end
end
%histogram((Zinta-Zinte))

CM = cbrewer('Set1',max(CC));

clf;
hold on;
tsurf(FF,VV,'FaceVertexCData',CM(CC,:),falpha(1,0),fsoft);

scatter3(XX(:,1),XX(:,2),XX(:,3),'.k','SizeData',100);
scatter3(P(:,1),P(:,2),P(:,3),'.r','SizeData',1000);
%tsurf(XE(XE(:,1)==1,:),XX,'LineWidth',3);
%quiver3(XX(1,1),XX(1,2),XX(1,3),XN(1,1),XN(1,2),XN(1,3),'LineWidth',3);
%tsurf(XE((XC(XE(:,1))==1) & (XC(XE(:,2))==3),:),XX);
tsurf(reshape(1:numel(XE),size(XE)),XX(XE,:),'CData',1*[Zint;Zint],'EdgeColor','flat')
%just = (XC(XE(:,1))==2)&(XC(XE(:,2))==3);
%just = Zint<0.2;
%tsurf(reshape(1:numel(XE(just,:)),size(XE(just,:))),XX(XE(just,:),:), ...
%  'CData',1*repmat(Zint(just),2,1),'EdgeColor','flat')
colormap(flipud(cbrewer('Greys',256)))
%tsurf(XE,XX);
hold off;
view(3);
axis equal;
l = add_lights;
%add_shadow([],l{5});
view(-82,18);
colorbar;
