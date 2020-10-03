rng(0);
clf

filenames = dir('data/meshes/backflip-boards/*.obj');
V = {};
F = {};
hold off
for fi = 0:15
  sprintf('data/meshes/backflip-boards/%1d.obj',fi)
  [Vi,Fi] = load_mesh(sprintf('data/meshes/backflip-boards/%1d.obj',fi));
%   V{end+1} = ...
%     off*clamp((fi-11)/(15-11),0,1) +  ...
%     Vi*axisangle2matrix([1 0 0],-pi/2)/1200*axisangle2matrix([0 0 1],pi/2);
%   F{end+1} = Fi;
  tsurf(Fi,Vi,'EdgeColor',[1 0 0]);
  view(90,0)
  set(gca,'visible','off')
  axis equal
  saveas(gca,sprintf('data/meshes/backflip-boards/%1d.png',fi));
%   figpng(sprintf('data/meshes/backflip-boards/%1d.png',fi));
end