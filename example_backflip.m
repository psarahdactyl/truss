%rng(0);
%
%filenames = dir('data/meshes/backflip/frame-*');
%V = {};
%F = {};
%% footskate issue (he doesn't land where he starts)
%off = [-0.25793  -3.4438e-07   2.5674e-05];
%for fi = 1:16
%  [Vi,Fi] = load_mesh(sprintf('data/meshes/backflip/frame-%02d.stl',fi));
%  [Vi,~,J] = remove_duplicate_vertices(Vi,0);Fi = J(Fi);
%  V{end+1} = ...
%    off*clamp((fi-11)/(15-11),0,1) +  ...
%    Vi*axisangle2matrix([1 0 0],-pi/2)/1200*axisangle2matrix([0 0 1],pi/2);
%  F{end+1} = Fi;
%end
%
%[SV,SF] = cube(2,2);
%SV = SV.*[1.7 0.02 1.0]*1.1 + [-0.9 2 0];
%
%CB = cbrewer('Set1',5);
%cbred = CB(1,:);
%cbblue = CB(2,:);
%cbgreen = CB(3,:);
%cborange = CB(5,:);
%
%[X,Y] = meshgrid(linspace(-1,1),linspace(-1,1));
%PP = normrow([X(:) Y(:)]-mean([X(:) Y(:)]));
%PP = clamp(max(PP)-PP,0.2,1);
%PP = matrixnormalize(PP);
%PP = reshape(PP,size(X));
%PX = X*0.7+0.0;
%PY = 0*Y-2;
%PZ = 0.35*Y+0.5;
%[PF,PV] = surf2patch(PX,PY,PZ,'triangles');
%nu = 200;
%%Q = blue_noise(nu*300,PV,PF);
%Q = random_points_on_mesh(PV,PF,nu*30);
%QP = interp2(PX,PZ,PP,Q(:,1),Q(:,3));
%[~,QI] = histc(rand(nu,1),[0;cumsum(QP)]/sum(QP));
%Q = Q(QI,:);
%
%CV = {};
%CF = {};
%ar = {};
%for ii = 1:numel(V)
%  fprintf('%d...\n',ii);
%  n = 1000;
%  VV = [SV;V{ii}];
%  FF = [SF;size(SV,1)+F{ii}];
%  CC = [ones(size(SF,1),1);1+ones(size(F{ii},1),1)];
%  [XX,XE,XC,YX,YE,YC] = groundstructure(VV,FF,CC,n);
%  r = 0.01;
%  Xvis = groundstructure_visibility(Q,VV,FF,XX,XE,'SampleSize',r);
%  sC = ones(size(XE,1),1)*1e3;
%  sT = sC;
%  sB = ones(size(XE,1),1)*1e2;
%  coms = [nan nan nan;centroid(V{ii},F{ii})];
%  [A,b,Aeq,beq] = create_constraint_matrices(XX,XE,XC,coms,force,sC,sT,sB);
%  wvis = 1000;
%  objective = edge_lengths(XX,XE) + wvis*Xvis.^2;
%  [x,ar{end+1},ax,be,fval] = optimize_lp(objective,A,b,Aeq,beq,'linprog');
%  NZ = find(max(ar{end},0)>1e-7);
%  [CV{end+1},CF{end+1}] = edge_cylinders(XX,XE(NZ,:),'Thickness',0.75.*sqrt(ar{end}(NZ)),'PolySize',30);
%
%%end
%%for ii = 1:numel(V)
%
%  clf;
%  hold on;
%psh = surf(PX,PY,PZ, ...
%  'CData',PP,'AlphaData',PP,fphong,'EdgeColor','none','FaceAlpha','interp');
%%psh = tsurf(PF,PV,'CData',PP,fphong,falpha(1,0),fsoft);
%colormap(interp1([1 0],[cbred;1 1 1],linspace(0,1,7)));
%psh.DiffuseStrength = 0;
%psh.AmbientStrength = 1;
%psh.SpecularStrength = 0;
%  ssh = tsurf(SF,SV,'FaceVertexCData',repmat(blue,size(SV,1),1),fphong,falpha(1,0),fsoft);
%  tsh = tsurf( ...
%    F{ii},V{ii},'FaceVertexCData',repmat(cbgreen,size(V{ii},1),1), ...
%    fphong,falpha(1,0),fsoft);
%  csh = tsurf(CF{ii},CV{ii}, ...
%    falpha(1,0),fsoft,'FaceVertexCData',repmat(cborange,size(CF{ii},1),1));
%  hold off;
%  axis equal;
%  axis([-0.9         0.97     -2.48651        2.022    -0.029161          1.1]);
%  view(-82,11);
%  l = { ...
%  light('Color',0.5*[1 1 1],'Position',(campos-camtarget),'Style','local'), ...
%  light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
%  light('Color',0.25*[1 1 1],'Position',-(campos-camtarget)*axisangle2matrix([0 0 1],pi/2),'Style','local'), ...
%  light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],pi*0.9),'Style','local'), ...
%  light('Color',0.25*[1 1 1],'Position', (campos-camtarget)*axisangle2matrix([1 0 0],-pi*0.9),'Style','local') ...
%  light('Color',0.25*[1 1 1],'Position', [-1 10 100],'Style','local') ...
%  };
%  add_shadow({ssh,tsh,csh},l{6});
%  %apply_ambient_occlusion([tsh,ssh,csh],'AddLights',false,'Factor',0.5);
%  camproj('persp');
%  %title(sprintf('%d',ii),'FontSize',30);
%  set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');
%  %figgif('backflip-side.gif');
%  %camtarget(coms(2,:));camup([0 0 1]);campos(mean(PV));camproj('persp');camva(50);
%  %figgif('backflip-view.gif');
%  drawnow
%
%end
%
%VV = [];
%FF = [];
%CC = [];
%th = linspace(0,2*pi,16+1);th = th(1:end-1)';
%BV = [cos(th) sin(th)];
%BE = [1:size(BV,1);2:size(BV,1) 1]';
%BC = barycenter(BV,BE);
%Bl = edge_lengths(BV,BE);
%Bd = normrow(BC);
%f = (max(SV(:,1))-min(SV(:,1)))/mean(Bl)*mean(Bd);
%
%NF = zeros(numel(V)+1,1);
%for ii = 1:numel(V)
%  off = [0 -2 0]+[0 -f 0];
%  T = @(V) (V+off)*axisangle2matrix([0 0 1],th(ii));
%  FF = [FF;size(VV,1)+F{ii}];
%  CC = [CC;3*ones(size(F{ii},1),1)];
%  VV = [VV;T(V{ii})];
%  FF = [FF;size(VV,1)+CF{ii}];
%  CC = [CC;5*ones(size(CF{ii},1),1)];
%  VV = [VV;T(CV{ii})];
%  FF = [FF;size(VV,1)+SF];
%  CC = [CC;2*ones(size(SF,1),1)];
%  VV = [VV;T(SV)];
%  NF(ii+1) = size(FF,1);
%end
%
%load('backflip-1000-0.01-1000.mat');
tsh = tsurf(FF,VV,'FaceVertexCData',CB(CC,:),fsoft,falpha(1,0));
AO = apply_ambient_occlusion([],'AddLights',false,'AO',AO,'Factor',0.5);
camlight;
axis equal;
axis([-7.1 7.1 -7.1 7.1 -0.1 1.2]);
set(gca,'Position',[0 0 1 1],'Visible','off');
set(gcf,'Color','w');
view(0,0);

camera_spf = 1/16/2;
shutter_speed = camera_spf/16;
subsamples = 10;
max_t = 3;
% radians per second of zoetrope
rps = 2*pi;
start = [0 -0.1 20];
mid = [0 -15 5];
finish = mean(PV)+off;

filename = 'backflip-camera.gif';
% lights out
tlo = 1;
for touter = 0:camera_spf:max_t-camera_spf
  if touter>tlo
    subsamples = 1;
    set(gcf,'Color',0.1*[1 1 1]);
  end
  blur = 0;
  if subsamples == 1
    ts = touter;
  else
    ts = linspace(touter-shutter_speed/2,touter+shutter_speed/2,subsamples);
  end
  ii = mod(floor(touter*16),16)+1;
  for t = ts
    tsh.Vertices = VV*axisangle2matrix([0 0 1],-rps*t);
    if touter>tlo
      df = 0.1;
      Ct = CB(CC,:)*df;
      if abs(mod(t,camera_spf*2))<1e-8
        Ct(NF(ii)+1:NF(ii+1),:) = Ct(NF(ii)+1:NF(ii+1),:)*(1/df);
      end
      tsh.FaceVertexCData = Ct.*(1-0.5*AO);
    end
    %view(t*360,0);
    cp = interp1([0 0.1 0.5 1 max_t],[start;start;mid;finish;finish],t,'pchip');
    camtarget([0 0 0]);camup([0 0 1]);campos(cp);camproj('persp');camva(50);drawnow;
    drawnow;
    frame = getframe(gcf);
    blur = blur + im2double(frame.cdata)/subsamples;
  end
  [SIf,cm] = rgb2ind(blur,256);
  f = exist(filename,'file');
  if ~f
    imwrite(SIf,cm,filename,'Loop',Inf,'Delay',0);
  else
    imwrite(SIf,cm, filename,'WriteMode','append','Delay',0);
  end
end

%view(112,40);
%l = light('Color',0.5*[1 1 1],'Position',[1 -1 40],'Style','local');
%ssh = add_shadow([],l);
%camtarget([0 0 0]);camup([0 0 1]);campos(finish);camproj('persp');camva(50);drawnow;
%tsh.Vertices = VV*axisangle2matrix([0 0 1],-rps*(11/16));
