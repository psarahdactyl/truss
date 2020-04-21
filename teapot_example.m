V = {};
F = {};

[SV,SF] = cube(2,2);
SV = SV.*[1 0.05 1] + [0 1 0];

V{1} = SV;
F{1} = SF;

[V{end+1},F{end+1}] = load_mesh('data/meshes/teapot.obj');
V{end} = V{end} *  ...
 axisangle2matrix([1 0 0],-pi*0.4) * axisangle2matrix([0 0 1],pi/5) * ...
   axisangle2matrix([0 1 0],pi/10);
for ii = 2:numel(V)
 % rescale to unit sphere
 V{ii} = V{ii}-0.5*(max(V{ii})+min(V{ii}));
 V{ii} = V{ii}/max(normrow(V{ii}));
 % 0.3 m wide
 V{ii} = V{ii}*0.8*0.5;
end
V{2} = V{2} + [0.6 0.5 0.5];

CB = cbrewer('Set1',5);
cbred = CB(1,:);
cbblue = CB(2,:);
cbgreen = CB(3,:);
cborange = CB(5,:);

% units are meters

xs = 2;
[X,Y] = meshgrid(linspace(0,1,xs),linspace(0,1,xs));
PP = normrow([X(:) Y(:)]-mean([X(:) Y(:)]));
PP = clamp(max(PP)-PP,0.2,1);
PP = matrixnormalize(PP);
PP = reshape(PP,size(X));
%[PF,PV] = surf2patch(X*0.75+0.25*0.5,0*Y-0.75,Y*0.5+0.25,'triangles');


objs = {V{:};F{:}};
objs=objs';
[AV,AF,ACV,ACF,coms,vols] = list_to_mesh(objs);
[RV,RF,~] = list_to_mesh(objs(2:end,:));

%%
x=X*0.75+0.25*0.5;
y=0*Y-0.75;
z=Y*0.5+0.25;
views = [x(:) y(:) z(:)];

%[Vs,GV,side,w] = scene_visibility(AV,AF,RV,RF,views);
%[GSV,E,VC] = construct_ground_structure(AV,AF,ACV,ACF);
%% bar lengths
%lengths = edge_lengths(GSV,E);
%% bar visibilities
%EVs = edge_visibilities(GSV,E,GV,side,w,Vs,lengths);
%% % function of visibility
%% fv = @(X) X.^3;
%% fv = @(X) exp(X);
%fv = @(X) (X);
%% bar projected areas
%[EAs] = edge_projected_visible_areas(GSV,E,views,EVs,fv);

n = 100;
[GSV,E,VC,~,Eraw] = groundstructure(AV,AF,ACF,n);
r = 0.01;
EAs = groundstructure_visibility(views,AV,AF,GSV,E,'SampleSize',r);

%%
force = [0 0 -9.8];

% yield stresses
sC = ones(size(E,1),1)*1e3;
% sC(1:size(E,1)/2,:) = 0;

sT = ones(size(E,1),1)*1e3;
% sT(size(E,1)/2+1:end,:) = 0;

sB = 0*ones(size(E,1),1)*1e3;
% sB(1:size(E,1)/2,:) = 0;


%%
% optimization
[A,b,Aeq,beq] = create_constraint_matrices(GSV,E,VC,coms,force,sC,sT,sB);

% [x,ar,ax,be] = optimize_lp(lengths,A,b,Aeq,beq,'linprog');
objective = edge_lengths(GSV,E) + 1000*EAs.^2;
[x,ar,ax,be,fval] = optimize_lp(objective,A,b,Aeq,beq,'linprog');
% objective = lengths;
% objective(size(E,1)/2:end) = ...
%           objective(size(E,1)/2:end) + EAs(size(E,1)/2:end);
% [x,ar,ax,be] = optimize_lp(objective,A,b,Aeq,beq,'linprog','bending',1); 

%%
save('bunny-teapot-rocker-arm.mat','AV','AF','ACV','ACF',...
      'ar','ax','be','sC','sT','sB','GSV','E','VC','EAs','views');
%%
clf;
hold on;

ssh = tsurf(SF,SV,'FaceVertexCData',repmat(blue,size(SV,1),1),fphong,falpha(1,0),fsoft);
red = [1 0 0];
psh = surf(X*0.75+0.25*0.5,0*Y-0.75,Y*0.5+0.25, ...
  'CData',PP,'AlphaData',PP,fphong,'EdgeColor','none','FaceAlpha','interp');
%psh = tsurf(PF,PV,'CData',PP,fphong,falpha(1,0),fsoft);
colormap(interp1([1 0],[cbred;1 1 1],linspace(0,1,7)));
psh.DiffuseStrength = 0;
psh.AmbientStrength = 1;
psh.SpecularStrength = 0;
tsh = {};
CM = interp1([0 1],[cbgreen;cbgreen*0.6+0.4],linspace(0,1,numel(V))');
for ii = 2:numel(V)
  tsh{end+1} = tsurf( ...
    F{ii},V{ii},'FaceVertexCData',repmat(CM(ii,:),size(V{ii},1),1), ...
    fphong,falpha(0.5,0),fsoft);
end

%plot_groundstructure(GSV,E,ar,ax);
NZ = ar>1e-6;
[CV,CF] = edge_cylinders(GSV,E(NZ,:),'Thickness',0.2*sqrt(ar(NZ)),'PolySize',30);

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
add_shadow({ssh,tsh{:}},l{6});
%apply_ambient_occlusion({ssh,tsh{:}},'AddLights',false);
camproj('persp');
%view(0,0);
%set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');
title(sprintf('%d out of %d, $%g',sum(NZ),size(E,1),fval),'FontSize',30);

