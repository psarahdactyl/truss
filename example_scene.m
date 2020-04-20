%V = {};
%F = {};
%[V{end+1},F{end+1}] = load_mesh('~/Dropbox/models/bunny-remesh.obj');
%V{end} = V{end} *  ...
%  axisangle2matrix([0 0 1],pi*0.1);
%[V{end+1},F{end+1}] = load_mesh('~/Dropbox/models/rocker-arm.obj');
%V{end} = V{end} *  ...
%  axisangle2matrix([1 1 0],-pi*0.4) * ...
%  axisangle2matrix([0 0 1],-pi*0.4);
%[V{end+1},F{end+1}] = load_mesh('~/Dropbox/models/teapot.obj');
%V{end} = V{end} *  ...
%  axisangle2matrix([1 0 0],-pi*0.4) * axisangle2matrix([0 0 1],pi/5) * ...
%    axisangle2matrix([0 1 0],pi/10);
%for ii = 1:numel(V)
%  % rescale to unit sphere
%  V{ii} = V{ii}-0.5*(max(V{ii})+min(V{ii}));
%  V{ii} = V{ii}/max(normrow(V{ii}));
%  % 0.3 m wide
%  V{ii} = V{ii}*0.3*0.5;
%end
%V{1} = V{1} + [0.6 0.5 0.8];
%V{2} = V{2} + [0.5 0.1 0.7];
%V{3} = V{3} + [0.6 0.6 0.5];

CB = cbrewer('Set1',3);
cbred = CB(1,:);
cbblue = CB(2,:);
cbgreen = CB(3,:);

[SV,SF] = cube(2,2);
SV = SV.*[1 0.05 1] + [0 1 0];

% units are meters

[X,Y] = meshgrid(linspace(0,1),linspace(0,1));
PP = normrow([X(:) Y(:)]-mean([X(:) Y(:)]));
PP = clamp(max(PP)-PP,0.2,1);
PP = matrixnormalize(PP);
PP = reshape(PP,size(X));
%[PF,PV] = surf2patch(X*0.75+0.25*0.5,0*Y-0.75,Y*0.5+0.25,'triangles');


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
for ii = 1:numel(V)
  tsh{end+1} = tsurf( ...
    F{ii},V{ii},'FaceVertexCData',repmat(CM(ii,:),size(V{ii},1),1), ...
    fphong,falpha(1,0),fsoft);
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
add_shadow({ssh,tsh{:}},l{6});
apply_ambient_occlusion([],'AddLights',false);
camproj('persp');
%view(0,0);
set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');

