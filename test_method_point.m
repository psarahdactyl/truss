addpath ./utils/
clear
close all
clf
hold on

[px,pz] = meshgrid(1:1:5,1:1:5);
px = px(:);
pz = pz(:);
% plane = [px, pz, zeros(length(px),1)];
plane = [px, zeros(length(px),1), pz];

% V = [3 3 .5;
%   plane];

V = [3 .5 3;
  plane];
f = zeros(size(V,1),3);
f(1,2) = -9.8;

% bf = find(V(:,3) == 0);
bf = find(V(:,2) == 0);

% generate views
% % x between -1 and 1, y between 0 and 1, z at 5.5
% [gx,gy] = meshgrid([-2.5:-.5:-5.5],[3:.5:5]);
% [gx,gy] = meshgrid([-3:-.5:-5],[3:.5:5]);
% [gx,gy] = meshgrid(-1:.2:1,0:.2:1);
% gx = gx(:);
% gy = gy(:);
% views = [gx, gy, 5.5*ones(length(gx),1)];
% views = [-4,4,4] % for debugging

% % compute visibility of scene
% [Vs,GV,side,w] = scene_visibility(AV,AF,RV,RF,views);

%%

E = [ones(size(V,1)-1,1) 1+(1:size(plane,1))']

% force = [0 -9.8 0];
% [V,E,f,bf,ig] = construct_ground_structure(AV,AF,ACV,ACF,coms,force);


% yield stresses
sC = 1e2;
sT = 1e2;
sB = 1e2;

%%
% bar lengths
lengths = edge_lengths(V,E);

% % bar visibilities
% EVs = edge_visibilities(V,E,GV,side,w,Vs,lengths);
% 
% % function of visibility
% % fv = @(X) X.^3;
% % fv = @(X) exp(X);
% fv = @(X) (X);
% 
% % bar projected areas
% [EAs,d,theta] = edge_projected_visible_areas(V,E,views,EVs,fv);


%%
% optimization
[B,C,BT,CT] = create_nodal_equilibrium_matrices(V,E,bf);
[A,b,Aeq,beq] = create_constraint_matrices(V,E,f,bf,sC,sT,sB,'Bending',1);
[x,ar,ax,be] = optimize_lp(lengths,A,b,Aeq,beq,'linprog');
% [x,ar,ax,be] = optimize_lp(lengths+EAs,A,b,Aeq,beq,'linprog','IgnoredEdges',ig);

% naa = size(objs{2,1},1);
% av = 1 + (1:naa)';
% ae = (1:naa)';
% dim = size(V,2);
% fa = reshape(BT(av + size(V,1)*[0:dim-1],naa+1:end)*ax(naa+1:end),[],dim);
% torque = normrow(sum(cross(fa,V(E(ae,2),:)-V(E(ae,1),:),2)))

% find areas bigger than a certain threshold and place a cylinder there
NZ = find(max(ar,0)>1e-7);
num_rods = size(NZ,1)
num_compression = sum(sign(ax(NZ))==1)
num_tension = sum(sign(ax(NZ))==-1)
  
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


plot_groundstructure(V,E,ar,ax);
% plot_edges(V,E,'Color','#848484','LineWidth',0.001);

% plot boundary conditions
scatter3(V(bf,1),V(bf,2),V(bf,3),'.b','SizeData',200);

% plot force
quiver3(V(:,1),V(:,2),V(:,3),f(:,1)*0.1,f(:,2)*0.1,f(:,3)*0.1,'r','LineWidth',3,'AutoScale','off')
    
% plot other vertices
scatter3(V(:,1),V(:,2),V(:,3),'.k','SizeData',50);

% % plot viewpoints
% scatter3(views(:,1),views(:,2),views(:,3),'.m','SizeData',100);

% % plot visibility values
% scatter3(GV(:,1),GV(:,2),GV(:,3),'.','CData', Vs,'SizeData',150);

% % plot edges with color
% [~, I] = sortrows(E);
% 
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
colormap(flipud(cbrewer('Dark2',256)))
% colormap(flipud(cbrewer('RdBu',256)))
title(...
  sprintf('Vol: %g, #bars: %d out of %i edges, sC/sT=%d, sB=%d',lengths'*ar,numel(NZ),size(E,1),sT,sB),...
  'FontSize',20);
% colorbar
