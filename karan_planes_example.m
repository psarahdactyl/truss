addpath ./utils/
% clear
% close all
clf
hold on

filename = ...
    "data/planes.txt";

% read scene and plot the objects
[objs,bb,views] = read_scene(filename);
[AV,AF,ACV,ACF,coms,vols] = list_to_mesh(objs);
[RV,RF,~] = list_to_mesh(objs(2:end,:));

% % % generate views
% % [gx,gy] = meshgrid([-8:-.5:-10],[1:.5:3]);
% [gx,gy] = meshgrid([-4:-.5:-6],[4:.5:6]);
% % [gx,gy] = meshgrid([-1:.5:1],[-1:.5:1]);
% 
% % [gx,gy] = meshgrid(-1:.2:1,0:.2:1);
% 
% gx = gx(:);
% gy = gy(:);
% views = [gx, gy, 5.5*ones(length(gx),1)];
% % views = [-4,4,4] % for debugging

% compute visibility of scene
[Vs,GV,side,w] = scene_visibility(AV,AF,RV,RF,views);

%%

force = [0 -9.8 0];
[V1,E1,VC1] = construct_ground_structure(AV,AF,ACV,ACF);
[V2,E2,VC2] = construct_ground_structure(AV,AF,ACV,ACF);

V = [V1;V2];
E = [E1;E2];
VC = [VC1;VC2];

% V =  V1;
% E = E1;
% VC = VC1;

% yield stresses
% sC = 1e2;
% sT = 1e2;
% sB = 1e3;

sC = ones(size(E,1),1)*1e2;
sC(1:size(E,1)/2,:) = 0;

sT = ones(size(E,1),1)*1e4;
sT(size(E,1)/2+1:end,:) = 0;

sB = ones(size(E,1),1)*1e3;
sB(1:size(E,1)/2,:) = 0;

%%
% bar lengths
lengths = edge_lengths(V,E);

% bar visibilities
EVs = edge_visibilities(V,E,GV,side,w,Vs,lengths);

% % function of visibility
% fv = @(X) X.^3;
% fv = @(X) exp(X);
fv = @(X) (X);

% bar projected areas
[EAs,d,theta] = edge_projected_visible_areas(V,E,views,EVs,fv);

%%
% optimization
[A,b,Aeq,beq] = create_constraint_matrices(V,E,VC,coms,force,sC,sT,sB,'bending',1);

% [x,ar,ax,be] = optimize_lp(lengths,A,b,Aeq,beq,'linprog');
% [x,ar,ax,be] = optimize_lp(lengths+EAs,A,b,Aeq,beq,'linprog');
objective = lengths;
objective(size(E,1)/2:end) = ...
          objective(size(E,1)/2:end) + EAs(size(E,1)/2:end);
[x,ar,ax,be] = optimize_lp(objective,A,b,Aeq,beq,'linprog','bending',1); 

NZ = find(max(ar,0)>1e-5);
%%
save('karan-planes.mat','AV','AF','ACV','ACF',...
      'ar','ax','be','sC','sT','sB','V','E','VC','EAs','views');

%%
% plotting
clf
hold on
% % plot background wall
% 
% tsurf(DF,DV,...
%     'FaceColor',0.5+0.5*blue,falpha(1,0),fsoft)
% 
% MC = repmat(0.8*[0.99 1 0.99],size(objs{2,2},1),1); % colors
% 
% % plot object
% tsurf(RF,RV,...
%     'FaceVertexCData',MC,falpha(0.8,0),fsoft);

% CM = cbrewer('Set1',max(AC));
CM = cbrewer('Blues',max(ACV));

tsurf(AF,AV,'FaceVertexCData',CM(ACF,:),falpha(0.6,0),fsoft);
  
plot_groundstructure(V,E,ar,ax);
% plot_edges(V,E,'Color','#848484','LineWidth',0.001);


% % plot force
% f = zeros(size(coms(2:end,:)));
% f(:,2) = -9.8;
% quiver3(coms(2:end,1),coms(2:end,2),coms(2:end,3),f(:,1)*0.1,f(:,2)*0.1,f(:,3)*0.1,'r','LineWidth',2,'AutoScale','off')
%     
% % plot other vertices
% scatter3(V(:,1),V(:,2),V(:,3),'.k','SizeData',50);

% plot viewpoints
scatter3(views(:,1),views(:,2),views(:,3),'.m','SizeData',100);

% % plot visibility values
% scatter3(GV(:,1),GV(:,2),GV(:,3),'.','CData', Vs,'SizeData',150);

% % plot edges with color
% [~, I] = sortrows(E);

% % nans=find(isnan(EVs));
% % 
% % % EVs(~isnan(EVs)) = min(EVs);
% % % EVs(isnan(EVs)) = max(EVs);
% % sevs=EVs(I);
% % 
% % plot(graph(E(nans',1),E(nans',2)),...
% %     'XData',V(:,1),'YData',V(:,2),'ZData',V(:,3),...
% %     'EdgeCData',sevs(nans'), 'LineWidth',.1);

% plot(graph(E(:,1),E(:,2)),...
%     'XData',V(:,1),'YData',V(:,2),'ZData',V(:,3),...
%     'EdgeCData',EVs(I), 'LineWidth',.1);

% plot(graph(E(:,1),E(:,2)),...
%     'XData',V(:,1),'YData',V(:,2),'ZData',V(:,3),...
%     'EdgeCData',lengths(I), 'LineWidth',2);
% 
% % caxis([5 max(lengths)])
% colormap(flipud(cbrewer('RdBu',256)))
% colorbar
% return

hold off;
axis equal;

% camlight;
camup([0 1 0]);
campos([-33 6 20])
cameratoolbar('SetCoordSys','y')
cameratoolbar('setmode','orbit')
camproj('perspective')
caxis([-1 1])
% colormap(flipud(cbrewer('Blues',256)))
% colormap(flipud(cbrewer('Dark2',256)))
colormap(flipud(cbrewer('RdBu',256)))
title(...
  sprintf('Vol: %g, #bars: %d out of %i edges',lengths'*ar,numel(NZ),size(E,1)),...
  'FontSize',20);
% colorbar

% view(180,0)
% view(-177,11)

