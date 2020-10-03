addpath ./utils/
clear
close all
clf
hold on

filename = ...
    "data/sphere.txt";

% read scene and plot the objects
[objs,bb] = read_scene(filename);
[AV,AF,ACV,ACF,coms] = list_to_mesh(objs);
[RV,RF,~] = list_to_mesh(objs(2:end,:));

% generate views
% [gx,gy] = meshgrid([-2.5:-.5:-5.5],[3:.5:5]);
% [gx,gy] = meshgrid([-3:-.5:-5],[3:.5:5]);
% [gx,gy] = meshgrid([-6:-.5:-8],[3:.5:5]);

% [gx,gy] = meshgrid([-8:-.5:-10],[1:.5:3]);
[gx,gy] = meshgrid([-4:-.5:-6],[4:.5:6]);
% [gx,gy] = meshgrid([-1:.5:1],[-1:.5:1]);

% [gx,gy] = meshgrid(-1:.2:1,0:.2:1);
gx = gx(:);
gy = gy(:);
views = [gx, gy, 5.5*ones(length(gx),1)];
% views = [-4,4,4] % for debugging

% compute visibility of scene
[Vs,GV,side,w] = scene_visibility(AV,AF,RV,RF,views);

%%

force = [0 -9.8 0];
[V,E,VC,bf] = construct_ground_structure(AV,AF,ACV,ACF);
% E = [E;E];
% size(V)
% size(E)
% return 

%   V = [1 3 1; 4 2 1;4 4 1;8 2 1;8 4 1]
%   E = [1 2;
%        1 3;
%        1 4;
%        1 5;
%        2 5;
%        2 4;
%        3 4;
%        3 5]
%   ACV = [1;2;2;3;3]
%   VC = [1;2;2;3;3]
%   coms = [1 1 1; 4 3 1; 8 3 1];
  
%   V = [1 3 1; 4 2 1;4 4 1;8 3 1]
%   E = [1 2;
%        1 3;
%        2 4;
%        3 4]
%   ACV = [1;2;2;3]
%   VC = [1;2;2;3]
%   coms = [1 1 1; 4 3 1; 8 3 1];

% yield stresses
sC = 1e2;
sT = 1e3;
sB = 1e3;

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
% [S,G,C,B] = create_nodal_equilibrium_matrices(V,E,VC,coms);
[A,b,Aeq,beq] = create_constraint_matrices(V,E,VC,coms,sC,sT,sB);

% m=size(E,1);
% TT = [sparse(size(B,1),m) B sparse(size(B,1),m)];
% Aeq=B;
% beq=[zeros(3,1);[0 -9.8 0]'];

[vx,var,vax,vbe] = optimize_lp(lengths,A,b,Aeq,beq,'linprog'); % visible
[x,ar,ax,be] = optimize_lp(lengths+EAs,A,b,Aeq,beq,'linprog'); % invisible

vNZ = find(max(var,0)>1e-5);
% num_rods = size(vNZ,1)
% num_compression = sum(sign(vax(vNZ))==1)
% num_tension = sum(sign(vax(vNZ))==-1)

NZ = find(max(ar,0)>1e-5);
% num_rods = size(NZ,1)
% num_compression = sum(sign(ax(NZ))==1)
% num_tension = sum(sign(ax(NZ))==-1)

% new edges
nE = [E(vNZ,:); E(NZ,:)];
[A,b,Aeq,beq] = create_constraint_matrices(V,nE,VC,coms,sC,sT,sB);

num_visible_rods = numel(vNZ);
num_invisible_rods = numel(NZ);

l = [lengths(vNZ); lengths(NZ)]; % new objective function

vt=find(vax(vNZ)>0); % visible tension indices
t=find(ax(NZ)>0)+num_visible_rods; % invisible tension indices

f=l;
f(vt)=f(vt)+EAs(vax(vNZ)>0);
f(t)=f(t)+EAs(ax(NZ)>0);

[nx,nar,nax,nbe] = optimize_lp(f,A,b,Aeq,beq,'linprog'); % new ground structure

nNZ = find(max(nar,0)>1e-5);
num_rods = size(nNZ,1)
num_compression = sum(sign(nax(nNZ))==1)
num_tension = sum(sign(nax(nNZ))==-1)

%%
% plotting
tiledlayout(3,1)

gs1.edges = E;
gs1.lengths = lengths;
gs1.NZ = vNZ;
gs1.areas = var;
gs1.axial = vax;
gs1.text = 'no visibility term'

gs2.edges = E;
gs2.lengths = lengths;
gs2.NZ = NZ;
gs2.areas = ar;
gs2.axial = ax;
gs2.text = 'with visibility term';

gs3.edges = nE;
gs3.lengths = l;
gs3.NZ = nNZ;
gs3.areas = nar;
gs3.axial = nax;
gs3.text = 'only compression rods have the term'

structures = {gs1,gs2,gs3};

for i=1:3
  
S = structures{i};
nexttile

% clf
hold on

% CM = cbrewer('Set1',max(AC));
CM = cbrewer('Blues',max(ACV));

tsurf(AF,AV,'FaceVertexCData',CM(ACF,:),falpha(0.6,0),fsoft);
  
if(~isempty(S.areas))
plot_groundstructure(V,S.edges,S.areas,S.axial);
end
% plot_edges(V,E,'Color','#848484','LineWidth',0.001);

% plot boundary conditions
bf = 1:size(V(VC==1));
scatter3(V(bf,1),V(bf,2),V(bf,3),'.b','SizeData',50);

% plot force
f = zeros(size(coms(2:end,:)));
f(:,2) = -9.8;
quiver3(coms(2:end,1),coms(2:end,2),coms(2:end,3),f(:,1)*0.1,f(:,2)*0.1,f(:,3)*0.1,'r','LineWidth',1,'AutoScale','off')
    
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
% colormap(flipud(cbrewer('Dark2',256)))
colormap(flipud(cbrewer('RdBu',256)))

if(~isempty(S.areas))
  string = sprintf('Vol: %g, #bars: %d out of %i edges',S.lengths'*S.areas,numel(S.NZ),size(S.edges,1));
else
  string = sprintf('Vol: %g, #bars: %d out of %i edges',0,numel(S.NZ),size(S.edges,1));
end
string_list = [S.text, newline, string];
title(...
  string_list,...
  'FontSize',20);
% colorbar
end

% view(180,0)
% view(-177,11)

