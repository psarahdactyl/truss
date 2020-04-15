addpath ./utils/
clear
close all
clf
hold on


objs = cell(2,2);
% [V,F] = readSTL('data/meshes/barbell.stl');
[V,F] = readOBJ('hammer.obj');
[SV,SVI,SVJ] = remove_duplicate_vertices(V,1e-7);
SF = SVJ(F);
objs{2,1} = SV;
objs{2,2} = SF;

[WV,WF] = create_regular_grid(2,2);
WV(:,3) = 0;
WV = WV*7;

objs{1,1} = WV;
objs{1,2} = WF;
size(objs)

[AV,AF,ACV,ACF,coms] = list_to_mesh(objs);
[RV,RF,~] = list_to_mesh(objs(2:end,:));

% generate views
% % x between -1 and 1, y between 0 and 1, z at 5.5
% [gx,gy] = meshgrid([-2.5:-.5:-5.5],[3:.5:5]);
% % [gx,gy] = meshgrid([-3:-.5:-5],[3:.5:5]);
% % [gx,gy] = meshgrid(-1:.2:1,0:.2:1);
% gx = gx(:);
% gy = gy(:);
% views = [gx, gy, 5.5*ones(length(gx),1)];
views = [4,2,-10] % for debugging

% compute visibility of scene
[Vs,GV,side,w] = scene_visibility(AV,AF,RV,RF,views);

%%
% DF = objs{1,2}; % wall mesh is always the first mesh
% DV = objs{1,1};
% 
% RV = objs{2,1}; % object mesh
% RF = objs{2,2};
% 
% % this creates the ground structure, (V,E), including the center of mass
% [I,J] = find(ones(size(RV,1),size(DV,1)));
% 
% % find centroid of object on which to apply force (gravity)
% cm = centroid(RV,RF);
% V = [cm;RV;DV];
% E = [ones(size(RV,1),1) 1+(1:size(RV,1))';1+[I size(RV,1)+J]];
% % size(E)

force = [0 -9.8 0];
[V,E,f,bf,ig] = construct_ground_structure(AV,AF,ACV,ACF,coms,force);
% size(V)
% size(E)
% return 

% % nodal forces, boundary vertices, stress limits
% f = zeros(size(V));
% bl = snap_points(cm,V); % find closest vertex in ground structure to the COM
% f(bl,2) = -9.8;
% bf = 1+size(RV,1)+(1:size(DV,1)); 
% ig = 1:size(RV,1);

% yield stresses
sC = 1e2;
sT = 1e2;
sB = 1e2;

%%
% bar lengths
lengths = edge_lengths(V,E);

% bar visibilities
EVs = edge_visibilities(V,E,GV,side,w,Vs,lengths);

% function of visibility
% fv = @(X) X.^3;
% fv = @(X) exp(X);
fv = @(X) (X);

% bar projected areas
[EAs,d,theta] = edge_projected_visible_areas(V,E,views,EVs,fv);


%%
% optimization
[B,C] = create_nodal_equilibrium_matrices(V,E,bf);
[A,b,Aeq,beq] = create_constraint_matrices(V,E,f,bf,sC,sT,sB);
% [x,ar,ax,be] = optimize_lp(lengths,A,b,Aeq,beq,'cvx','IgnoredEdges',ig);
[x,ar,ax,be] = optimize_lp(lengths+EAs,A,b,Aeq,beq,'yalmip','IgnoredEdges',ig);

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

% CM = cbrewer('Set1',max(AC));
CM = cbrewer('Greens',max(ACV));

tsurf(AF,AV,'FaceVertexCData',CM(ACF,:),falpha(0.6,0),fsoft);
  
plot_groundstructure(V,E,ar,ax);
% plot_edges(V,E,'Color','#848484','LineWidth',0.001);

% plot boundary conditions
scatter3(V(bf,1),V(bf,2),V(bf,3),'.b','SizeData',200);

% plot force
quiver3(V(:,1),V(:,2),V(:,3),f(:,1)*0.1,f(:,2)*0.1,f(:,3)*0.1,'r','LineWidth',3,'AutoScale','off')
    
% plot other vertices
scatter3(V(:,1),V(:,2),V(:,3),'.k','SizeData',50);

% plot viewpoints
scatter3(views(:,1),views(:,2),views(:,3),'.m','SizeData',100);

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
  sprintf('Vol: %g, #bars: %d out of %i edges',lengths'*ar,numel(NZ),size(E,1)),...
  'FontSize',20);
% colorbar
