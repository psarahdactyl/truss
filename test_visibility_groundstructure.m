%addpath ../cyCodeBase/
%[V,F] = load_mesh('~/Dropbox/models/decimated-knight.off');
%rng(0);
%[VV,FF,~,CC] = repmesh(V,F,rand(3,3));
%[WV,WF] = create_regular_grid(2,2);
%WV(:,3) = 0;
%WV = WV*2;
%FF = [FF;size(VV,1)+WF];
%CC = [CC;repmat(max(CC)+1,size(WF,1),1)];
%VV = [VV;WV];
%R = axisangle2matrix([1 0 0],-pi/2);
%VV = VV*R;
%
%n = 400;
%[XX,XE,XC] = groundstructure(VV,FF,CC,n);
%
%% edge vector
%XEV = XX(XE(:,2),:)-XX(XE(:,1),:);
%XEU = normalizerow(XEV);
%
%P = [1 -3 1];
%r = 0.1;
%[PV,PF] = create_regular_grid(20,2);
%PV = PV-0.5;
%w = 2;
%PV = [w*[1 1] w/3].*[[sin(pi/8*(PV(:,1))) -cos(pi/8*(PV(:,1)))] PV(:,2)] + [1 -3+w 1];
%P = blue_noise(100,PV,PF);
%W = exp(-10*normrow(P-mean(P)));
%W = W-min(W);
%W = W./sum(W);
%
%for iter = 2
%  switch iter
%  case 1
%    vis_method = 'approx-sampling';
%  case 2
%    vis_method = 'exact-sampling';
%  end
%  tic;
%  Zint = groundstructure_visibility( ...
%    P,VV,FF,XX,XE,'SampleSize',r,'Method',vis_method,'Weight',W);
%  fprintf('%s: %g secs\n',vis_method,toc);
%  switch vis_method
%  case 'approx-sampling';
%    Zinta = Zint;
%  case 'exact-sampling';
%    Zinte = Zint;
%  end
%end
%%histogram((Zinta-Zinte))


CM = cbrewer('Set1',max(CC));

clf;
hold on;
tsurf(FF,VV,'FaceVertexCData',CM(CC,:),falpha(1,0),fsoft);

scatter3(XX(:,1),XX(:,2),XX(:,3),'.k','SizeData',100);
scatter3(P(:,1),P(:,2),P(:,3),'.r','SizeData',1000*(W/max(W))+1);
tsurf(PF,PV,'FaceColor',[1 0.8 0.8],falpha(0.5,0));
%tsurf(XE(XE(:,1)==1,:),XX,'LineWidth',3);
%quiver3(XX(1,1),XX(1,2),XX(1,3),XN(1,1),XN(1,2),XN(1,3),'LineWidth',3);
%tsurf(XE((XC(XE(:,1))==1) & (XC(XE(:,2))==3),:),XX);
tsurf(reshape(1:numel(XE),size(XE)),XX(XE,:),'CData',1*[Zint;Zint],'EdgeColor','flat')
%just = (XC(XE(:,1))==2)&(XC(XE(:,2))==3);
%just = Zint<0.1;
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
view(0,0);
