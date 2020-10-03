rng(0)
% clear
% clf

addpath utils
addpath ../cyCodeBase 
V = {};
F = {};
coms = [nan nan nan];
vols = nan;
masses = nan;

[V{end+1},F{end+1}] = load_mesh('data/meshes/claymation/bear.obj');
V{end}=V{end}/100;
coms(end+1,:) = centroid(V{end},F{end});

[V{end+1},F{end+1}] = load_mesh('data/meshes/claymation/crane.obj');
V{end}=V{end}/100;
coms(end+1,:) = centroid(V{end},F{end});

[V{end+1},F{end+1}] = load_mesh('data/meshes/claymation/flamingo.obj');
V{end}=V{end}/100;
coms(end+1,:) = centroid(V{end},F{end});

[V{end+1},F{end+1}] = load_mesh('data/meshes/claymation/penguin.obj');
V{end}=V{end}/100;
coms(end+1,:) = centroid(V{end},F{end});

[V{end+1},F{end+1}] = load_mesh('data/meshes/claymation/turtle.obj');
V{end}=V{end}/100;
coms(end+1,:) = centroid(V{end},F{end});

masses = [nan 1.9 0.1 1.3 6.4 3.7]';
masses = masses/100;
% paper crane 0.1g
% turtle 3.7g
% penguin 6.4g
% bear 1.9g
% flamingo 1.3g
%%
CB = cbrewer('Set1',5);
cbred = CB(1,:);
cbblue = CB(2,:);
cbgreen = CB(3,:);
cborange = CB(5,:);

[SV,SF] = load_mesh('data/meshes/claymation/wall.obj');
SV=SV/100;
%%
% units are meters
[X,Y] = meshgrid(linspace(0,1),linspace(0,1));
PP = normrow([X(:) Y(:)]-mean([X(:) Y(:)]));
PP = clamp(max(PP)-PP,0.2,1);
PP = matrixnormalize(PP);
PP = reshape(PP,size(X));
PX = X*0.9-0.8;
PY = 0*Y+1.5;
PZ = Y*0.5+1.25;
[PF,PV] = surf2patch(PX,PY,PZ,'triangles');
tsurf(PF,PV);

nu = 200;
%Q = blue_noise(nu*300,PV,PF);
Q = random_points_on_mesh(PV,PF,nu*30);
QP = interp2(PX,PZ,PP,Q(:,1),Q(:,3));
[~,QI] = histc(rand(nu,1),[0;cumsum(QP)]/sum(QP));
Q = Q(QI,:);
% Q = Q *  ...
%   axisangle2matrix([0 0 1],-pi/4);
% AP=[PX(:) PY(:) PZ(:)];
% AP = AP *  ...
%   axisangle2matrix([0 0 1],-pi/4);
% PX = reshape(AP(:,1),size(PX));
% PY = reshape(AP(:,2),size(PY));
% PZ = reshape(AP(:,3),size(PZ));

% QX = X*0.75-0.25;
% QY = 0*Y-1.5;
% QZ = Y*1.5-0.25;
% [RF,RV] = surf2patch(QX,QY,QZ,'triangles');
% R = random_points_on_mesh(RV,RF,nu*30);
% RP = interp2(QX,QZ,PP,R(:,1),R(:,3));
% [~,RI] = histc(rand(nu,1),[0;cumsum(RP)]/sum(RP));
% R = R(RI,:);
%%
n = 1500;

VV = cell2mat(reshape(V,[],1));
nV = cumsum([0 cellfun(@(V) size(V,1),{V{:}})]);
FF = cell2mat(arrayfun(@(i) nV(i)+F{i},1:numel(V),'UniformOutput',false)');
CC = cell2mat(arrayfun(@(i) repmat(i,size(F{i},1),1),1:numel(F),'UniformOutput',false)');
VV = [SV;VV];
FF = [SF;size(SV,1)+FF];
CC = [ones(size(SF,1),1);1+CC];
[XX,XE,XC,YX,YE,YC] = groundstructure(VV,FF,CC,n);
% don't prune
% XX=YX;
% XE=YE;
% XC=YC;

% for wvis = [10000]

  %%
  r = 0.1;
%   if wvis == 0
%     Xvis = zeros(size(XE,1),1);
%   else
%     Xvis = groundstructure_visibility([Q;R],VV,FF,XX,XE,'SampleSize',r);
%   end
%   % Xvis = groundstructure_visibility(Q,VV,FF,XX,XE,'SampleSize',r);
%   % XE=YE;
  Xvis = groundstructure_visibility(Q,VV,FF,XX,XE,'SampleSize',r);

%   
  %%
  g = [0 0 -9.8];
%   wires = wvis>0;
  wires = 0;
  %wvis = 0;
  if wires
    XE = [XE;XE];
    Xvis = [Xvis;0*Xvis];
    % tension wires are second half
    sC = [ones(size(XE,1)/2,1) zeros(size(XE,1)/2,1)]*1e4;
    sT = [zeros(size(XE,1)/2,1) ones(size(XE,1)/2,1)]*1e4;
    sB = [ones(size(XE,1)/2,1) zeros(size(XE,1)/2,1)]*1e4;
  else
    sC = ones(size(XE,1),1)*1e3;
    sT = sC;
    sB = sC;%*0.1;
  end
  
  [A,b,Aeq,beq] = create_constraint_matrices(XX,XE,XC,coms,g,sC,sT,sB,'Mass',masses);
  objective = edge_lengths(XX,XE)+ 100*Xvis.^2;
%   objective(sT>0) = objective(sT>0)*10;
  [x,ar,ax,be,fval] = optimize_lp(objective,A,b,Aeq,beq,'linprog');
  
  NZ = find(max(ar,0)>1e-7);
  if wires
    RZ = ar>1e-6 & ax>0;
    WZ = ar>1e-6 & ax<0;
  else
    RZ = NZ;
    WZ = [];
  end
  
%%
%plot_groundstructure(GSV,E,ar,ax);
NZ = ar>1e-6;
[CV,CF] = edge_cylinders(XX,XE(NZ,:),'Thickness',sqrt(ar(NZ)/pi),'PolySize',30);

% end
%%
clf;
hold on;

ssh = tsurf(SF,SV,'FaceVertexCData',repmat(cbblue,size(SV,1),1),fphong,falpha(1,0),fsoft);
red = [1 0 0];
psh = surf(PX,PY,PZ, ...
  'CData',PP,'AlphaData',PP,fphong,'EdgeColor','none','FaceAlpha','interp');
colormap(interp1([1 0],[cbred;1 1 1],linspace(0,1,7)));
psh.DiffuseStrength = 0;
psh.AmbientStrength = 1;
psh.SpecularStrength = 0;
tsh = {};
CM = interp1([0 1],[cbgreen;cbgreen*0.6+0.4],linspace(0,1,numel(V))');
for ii = 2:numel(V)
  tsh{end+1} = tsurf( ...
    F{ii},V{ii},'FaceVertexCData',repmat(CM(ii,:),size(V{ii},1),1), ...
    fphong,falpha(1,0),fsoft);
end

csh = tsurf(CF,CV, ...
  falpha(1,0),fsoft,'FaceVertexCData',repmat(cborange,size(CV,1),1));

hold off;
axis equal;
view(-58,15);

l = { ...
 light('Color', 0.5*[1 1 1],'Position',(campos-camtarget),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position',-(campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],pi*0.9),'Style','local'), ...
 light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],-pi*0.9),'Style','local') ...
 light('Color',0.25*[1 1 1],'Position', [-1 -1 10],'Style','local') ...
};
AO = apply_ambient_occlusion([ssh,tsh{:},csh],'AddLights',false);
% add_shadow({ssh,tsh{:}},l{6});
% apply_ambient_occlusion({ssh,tsh{:}},'AddLights',false);
camproj('persp');
%view(0,0);
set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');
% title(sprintf('%d out of %d, $%g',sum(NZ),size(E,1),fval),'FontSize',30);


  
%   if usejava('jvm')
%     [CV,CF] = edge_cylinders(XX,XE(RZ,:),'Thickness',1.*sqrt(ar(RZ)),'PolySize',30);
%     [uV,uF] = mesh_boolean(SV,SF,[],[],'union');
%     [uV,uF] = upsample(SV,SF,'Iterations',10,'OnlySelected',@(V,F) find(doublearea(V,F)>0.001));
%     
%     %
%     clf;
%     hold on;
%     ssh = tsurf(uF,uV,'FaceVertexCData',repmat(cbblue,size(uF,1),1),falpha(1,0),fsoft);
%     
%     % viewpoint distributions
%     psh = surf(QX,QY,QZ, ...
%       'CData',PP,'AlphaData',PP,fphong,'EdgeColor','none','FaceAlpha','interp');
%     %psh = tsurf(PF,PV,'CData',PP,fphong,falpha(1,0),fsoft);
%     colormap(interp1([1 0],[cbred;1 1 1],linspace(0,1,8)));
%     psh.DiffuseStrength = 0;
%     psh.AmbientStrength = 1;
%     psh.SpecularStrength = 0;
%     
%     rsh = surf(PX,PY,PZ, ...
%       'CData',PP,'AlphaData',PP,fphong,'EdgeColor','none','FaceAlpha','interp');
%     %rsh = tsurf(PF,PV,'CData',PP,fphong,falpha(1,0),fsoft);
%     colormap(interp1([1 0],[cbred;1 1 1],linspace(0,1,8)));
%     rsh.DiffuseStrength = 0;
%     rsh.AmbientStrength = 1;
%     rsh.SpecularStrength = 0;
%     
%     vV = [];
%     vF = [];
%     vC = [];
%     for ii = 1:numel(V)
%       vF = [vF;size(vV,1)+F{ii}];
%       vC = [vC;repmat(ii,size(F{ii},1),1)];
%       vV = [vV;V{ii}];
%     end
%     
%     vAO = ambient_occlusion(vV,vF,barycenter(vV,vF),normalizerow(normals(vV,vF)),1000);
%     vF = vF(vAO<0.99,:);
%     vC = vC(vAO<0.99,:);
%     vAO = vAO(vAO<0.99);
%     
%     [vV,vF,J] = upsample(vV,vF,'Iterations',10,'OnlySelected',@(V,F) find(doublearea(V,F)>0.001));
%     vC = vC(J);
%     
%     tsh = {};
%     %%% attempting to make more distinct green colors
%     CM = interp1([0 1],[cbgreen;cbgreen*0.8+0.6],linspace(0,1,numel(V))');
%     rng(8);
%     CM = CM(randperm(size(CM,1)),:);
%     tsh{end+1} = tsurf( ...
%       vF,vV,'FaceVertexCData',CM(vC,:), ...
%       falpha(1,0),fsoft);
%     csh = tsurf(CF,CV, ...
%       falpha(1,0),fsoft,'FaceVertexCData',repmat(cborange,size(CF,1),1));
%     wsh = tsurf(XE(WZ,:),XX,'LineWidth',0.4,falpha(0,0.6),'EdgeColor','k','FaceColor','none');
%     
%     hold off;
%     axis equal;
%     view(-141,-1);
%     
%     l = { ...
%      light('Color',0.5*[1 1 1],'Position',(campos-camtarget),'Style','local'), ...
%      light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
%      light('Color',0.25*[1 1 1],'Position',-(campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
%      light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],pi*0.9),'Style','local'), ...
%      light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],-pi*0.9),'Style','local') ...
%      light('Color',0.25*[1 1 1],'Position', [-1 -1 20],'Style','local') ...
%     };
% %     apply_ambient_occlusion(tsh,'AddLights',false);
% %     apply_ambient_occlusion(ssh,'AddLights',false);
% %     % apply_ambient_occlusion([tsh{:},ssh],'AddLights',false);
% %     % apply_ambient_occlusion(csh,'AddLights',false,'Factor',0.5);
% %     add_shadow({ssh,tsh{:},csh},l{6});
% %     %apply_ambient_occlusion([],'AddLights',false);
%     camproj('persp');
%     %%view(0,0);
%     %
%     %% the loopy animation needs some work :(
%     %% for t = linspace(0,2*pi,60);camtarget(mean(VV));camup([0 0 1]);campos(mean(Q)+0.05*[((1+pi.*cos(t)).*cos(t)-pi+1)/2 0 0.5*(1+pi.*cos(t)).*sin(t)/2]);camproj('persp');camva(50);drawnow;end
%     
%     set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');
%     for pass = 1:3
%       switch pass
%       case 1
%         psh.Visible = 'off';
%         rsh.Visible = 'off';
%       %%% outer view point for teaser
%       view(-175,6.0236);camva(7)
%       case 2
%         psh.Visible = 'off';
%         rsh.Visible = 'off';
%       %%% viewpoint for 1 distribution
%       camtarget(mean(VV)-[0.3 0 0]);camup([0 0 1]);campos(mean(R)-[0 .3 0]);camproj('persp');camva(80);
%       case 3
%       %%% viewpoint for 2 distribution
%         psh.Visible = 'off';
%         rsh.Visible = 'off';
%       camtarget(mean(VV)-[0.3 0 0]);camup([0 0 1]);campos(mean(Q)-[0.3 0 0]);camproj('persp');camva(65);
%       end
%       figpng(sprintf('rocket-%d-%d.png',wvis,pass));
%     end
% 
%   end
% end
% 
% if usejava('jvm')
%   csh.Visible='off';
%   psh.Visible = 'on';
%   rsh.Visible = 'on';
%   %%% outer view point for teaser
%   view(-175,6.0236);camva(7)
%   figpng('rocket-input.png');
% end
% 
% %%
% % save('teaser.mat','XX','XE','RZ','WZ','Q','R',...
% %       'x','ar','ax','be','fval','sC','sT','sB','Xvis');
