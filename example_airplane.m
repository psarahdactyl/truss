rng(0);

addpath ../cyCodeBase/
addpath utils/
[V,F] = load_mesh('data/meshes/airplane/wwi-plane-triangle-smooth.obj');
[~,C] = connected_components(F);
V = V*axisangle2matrix([1 0 0],pi/2);
[V,F,uJ] = mesh_boolean(V,F,[],[],'union');
C = C(uJ);
V = V - 0.5*(max(V)+min(V));
V = V/max(abs(V(:)));

[SV,SF] = create_regular_grid(2,2);
SV = 1*(SV*[1 0 0;0 1 0]+[-0.5 -0.5 0])+[0 0 1.5];
SF=fliplr(SF);

Q = [0 0 -1];

coms = [nan nan nan;centroid(V,F)];

force = [0 0 -9.8];

%for pass = 2
for pass = 1:2
  switch pass
  case 1
    % first shoot at SV,SF
    [~,t] = ray_mesh_intersect(coms(2,:),-force,SV,SF);
    XX = [coms(2,:)-t*force];
    % shoot back at 
    [~,t] = ray_mesh_intersect(XX(1,:),force,V,F);
    XX = [XX;XX(1,:)+t*force];
    XE = [1 2];
    NZ = 1;
  case 2
    VV = [SV;V];
    FF = [SF;size(SV,1)+F];
    CC = [ones(size(SF,1),1);1+ones(size(F,1),1)];
    n = 20;
    [XX,XE,XC,YX,YE,YC] = groundstructure(VV,FF,CC,n);
    r = 0.1;
    Xvis = groundstructure_visibility(Q,VV,FF,XX,XE,'SampleSize',r);
    force = [0 0 -9.8];
    sT = ones(size(XE,1),1)*1e3;
    sC = ones(size(XE,1),1)*0;
    sB = ones(size(XE,1),1)*0;
    % symmetry
    sp = [0 0 0];
    sn = [1 0 0];
    XE = [XE;size(XX,1)+XE];
    XX = [XX;XX-2*sum((XX-sp).*sn,2).*sn];
    XC = repmat(XC,2,1);
    sC = repmat(sC,2,1);
    sT = repmat(sT,2,1);
    sB = repmat(sB,2,1);
    Xvis = repmat(Xvis,2,1);
    wvis = 1000;
    [A,b,Aeq,beq] = create_constraint_matrices(XX,XE,XC,coms,force,sC,sT,sB);
    objective = edge_lengths(XX,XE) + wvis*Xvis.^2;
    [x,ar,ax,be,fval] = optimize_lp(objective,A,b,Aeq,beq,'linprog');
    NZ = find(max(ar,0)>1e-7);
    [ar(NZ) ax(NZ) be(NZ,:)]
  end

  
  CB = cbrewer('Set1',5);
  cbred = CB(1,:);
  cbblue = CB(2,:);
  cbgreen = CB(3,:);
  cborange = CB(5,:);
  
  clf;
  hold on;
  ssh = tsurf(SF,SV,'FaceVertexCData',repmat(cbblue,size(SF,1),1),falpha(1,0),fsoft);
  tsh = tsurf( ...
    F,V,'FaceVertexCData',repmat(cbgreen,size(F,1),1), ...
    falpha(1,0),fsoft);
  wsh = tsurf(XE(NZ,:),XX,'LineWidth',1,'EdgeColor','k','FaceColor','none');
  hold off;
  axis equal;
  camproj('persp');
  view(-119,2);
  l = { ...
   light('Color',0.5*[1 1 1],'Position',(campos-camtarget),'Style','local'), ...
   light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
   light('Color',0.25*[1 1 1],'Position',-(campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
   light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],pi*0.9),'Style','local'), ...
   light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],-pi*0.9),'Style','local') ...
   light('Color',0.25*[1 1 1],'Position', [-1 10 500],'Style','local') ...
  };
  add_shadow({ssh,tsh},l{6},'Nudge',0.5);
  apply_ambient_occlusion([ssh,tsh],'AddLights',false);
  set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');
  figpng(sprintf('example-airplane-%d.png',pass));
end
