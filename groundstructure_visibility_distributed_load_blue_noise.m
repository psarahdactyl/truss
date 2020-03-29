% this file will be using point loads instead of rigid body loads
clf
clear
close all
filename = ...
    "~/Documents/siggraph2020/hidden-supports-application/viewer/data/scenes/canonical.txt";


% read scene and plot the objects
[objs,bb] = read_scene(filename);
[AV,AF] = list_to_mesh(objs);

% create voxel grid
[GV,side,w] = voxel_grid(AV,20);

% % read grid
% GV = readDMAT("~/Documents/siggraph2020/truss-opt/grid.dmat");
% voxel_width = GV(1,1) - GV(2,1);
% GV = GV - abs(repmat(repelem(voxel_width/2,3),size(GV,1),1));

% % read visibility values of GV for each object
% Vs = readDMAT("~/Documents/siggraph2020/truss-opt/grid_visibilities.dmat");
% side = readDMAT("side.dmat");

% generate views
s = rng(0);
num_views = 100;
r = normrnd(4,1.5,[num_views,2]);
views = [-r(:,1) r(:,2) 4*ones(size(r,1),1)];
% views = [-4,4,4] % for debugging
% plot views
scatter3(views(:,1),views(:,2),views(:,3),'.m','SizeData',100);

% calculate visibilities
src = repmat(views,size(GV,1),1);
dir = repelem(GV,size(views,1),1);
% ve = [src;dir];
% s = size(src,1)
% d = size(dir,1)
% ee = [(1:s)' ((s+1):(s+d))'];
% plot_edges(ve,ee);
Vs = ray_mesh_intersect(src,dir-src,AV,AF);
Vs(Vs ~= 0) = 1;
Vs = sum(reshape(Vs,size(views,1),size(GV,1)),1);
Vs = Vs ./ size(views,1);

%%
hold on
axis equal

DF = objs{1}{2};
DV = objs{1}{1};
% plot background wall
tsurf(DF,DV,...
    'FaceColor',0.5+0.5*blue,falpha(1,0),fsoft)

MC = repmat(0.8*[0.99 1 0.99],size(objs{2}{2},1),1); % colors

% hidden-ness/visibility
H = Vs;

RV = objs{2}{1}; % object mesh
RF = objs{2}{2};
% plot object
tsurf(RF,RV,...
    'FaceVertexCData',MC,falpha(0.8,0),fsoft);

% sample points on surfaces
addpath ../cyCodeBase/
[Rc,Irc,Brc] = blue_noise(3*sum(doublearea(RV,RF)),RV,RF,'Seed',rand);
[Dc,Idc,Bdc] = blue_noise(3*sum(doublearea(DV,DF)),DV,DF,'Seed',rand);
scatter3(Rc(:,1),Rc(:,2),Rc(:,3),'.g','SizeData',100);
scatter3(Dc(:,1),Dc(:,2),Dc(:,3),'.y','SizeData',100);

% this creates the ground structure, (V,E), including the center of mass
[I,J] = find(ones(size(Rc,1),size(Dc,1)));

% find centroid of object on which to apply force (gravity)
cm = centroid(RV,RF);
V = [cm;Rc;Dc];
% V = [Rc;Dc];
E = [ones(size(Rc,1),1) 1+(1:size(Rc,1))';1+[I size(Rc,1)+J]];
% E = [I size(Rc,1)+J];

% % prune edges
% num_nodes = size(V,1)
% num_edges = size(E,1)
% 
% [flag, t, lambda] = ray_mesh_intersect(V(E(size(RV,1):end,1),:), V(E(size(RV,1):end,2),:), RV, RF);
% num_pruned = sum(t<1 & t>0)
% 
% E(logical(t<0.99999 & t>0),:) = [];

% nodal forces, boundary vertices, stress limits
f = zeros(size(V));
bl = snap_points(cm,V);
f(bl,2) = -9.8;
bf = 1+size(Rc,1)+(1:size(Dc,1)); 
% bf = size(Rc,1)+(1:size(Dc,1));
% f(1,2) = -9.8;

sC = 5e2;
sT = 5e2;

%%
% bar lengths
lengths = edge_lengths(V,E);

% bar visibilities
[h, segments] = edge_visibilities(V,E,GV,side,w,H,lengths);
v = h;

% bar projected areas
% average distances from viewpoints to bar midpoints
midpoints = (V(E(:,1),:) - V(E(:,2),:))/2;
vs = repmat(views,size(E,1),1);
ms = repelem(midpoints,size(views,1),1);
d = vecnorm(ms-vs,2,2);
d = mean(reshape(d,size(views,1),size(E,1)),1)';
% bar thetas
edges = repelem([V(E(:,1),:) V(E(:,2),:)],size(views,1),1);
v1s = edges(:,1:3)-vs;
v2s = edges(:,4:6)-vs;
normcostheta = dot(v1s,v2s,2);
costheta = normcostheta ./ (vecnorm(v1s,2,2).*vecnorm(v2s,2,2));
meantheta = mean(reshape(costheta,size(views,1),size(E,1)),1)';
theta = acos(meantheta) * (180 / pi );
% [~, I] = sortrows(E);
% plot(graph(E(:,1),E(:,2)),...
%     'XData',V(:,1),'YData',V(:,2),'ZData',V(:,3),...
%     'EdgeCData',lengths(I), 'LineWidth',2);
% caxis([min(lengths) max(lengths)])
% colormap(flipud(cbrewer('RdBu',256)))
% colorbar

pva = theta.^2 .* v'.^2 ./ d.^2;

%%
%%%%%%%%%%%
% % get areas, axial forces, lengths, visibililities, nodal equilibrium matrix
% [a,n,BT] = solve_groundstructure_lp(V,E,f,bf,...
%                 sC,sT,v,lengths,'IgnoredEdges',1:size(Rc,1));

%%%%%%%%%%%

%%
%%%%%%%%%%%
    dim = size(V,2);
    nf = size(f,3);
  
    EV = normalizerow(V(E(:,2),:)-V(E(:,1),:));

    obj = lengths;
    objv = lengths.*(v+0.1)';

    % use ignored edges for free
    ignored_edges = 1:size(Rc,1);
    obj(ignored_edges) = 0;
    objv(ignored_edges) = 0;

    nnn = size(V,1);
    m = size(E,1);

    BT = sparse(E(:)+nnn*(0:dim-1),repmat(1:m,dim,2)',[EV;-EV],dim*nnn,m);
    OBT = BT;

    bf = bf(:);
    % remove fixed values
    BT(bf+nnn*(0:dim-1),:) = [];
    ff = f;
    ff(bf,:,:) = [];

    max_a = inf;
  
    BT = repdiag(BT,nf);
    II = speye(m*nf,m*nf);
    I = speye(m,m);
    Z = sparse(m,m);
    A = [repmat(-sC*I,nf,1) , -II; ...
        repmat(-sT*I,nf,1) ,  II];
    b = zeros(2*m*nf,1);

    vec = @(X) X(:);
    ig_ind = vec(ignored_edges(:)+m*(0:1))+m*(0:nf-1);
    A(ig_ind,:) = [];
    b(ig_ind,:) = [];

    As = [Z, repmat(-sC*I,nf,1) , -II; ...
        Z, repmat(-sT*I,nf,1) ,  II];
    bs = zeros(2*m*nf,1);
    
    As(ig_ind,:) = [];
    bs(ig_ind,:) = [];

      
%     P = sparse(1:m,1:m,pi*lengths,m,m);
%     Q = blkdiag(P,sparse(m,m));
%     [ Qsqrt, p, QS  ] = chol( Q );

%     P0 = -sC*I;
% %     rid = find(sign(EV(:,2))==1);
% %     P0(rid,:) = 0;
%     P1 = -sT*I;
% %     rid2 = find(sign(EV(:,2))==1);
% %     P1(rid2,:) = 0;
% 
%     Q0 = blkdiag(P0,sparse(m,m));
%     Q1 = blkdiag(P1,sparse(m,m));
%     
%     q0 = [zeros(m,1); 1*ones(m,1)];
%     q1 = [zeros(m,1); 1*ones(m,1)];
    
%     q0(ig_ind,:) = [];
%     q1(ig_ind,:) = [];
    
%     [ Q0sqrt, p, Q0S  ] = chol( Q0 );
%     Q0sqrt = -Q0sqrt;% * Q0S;
    
%     [ Q1sqrt, p, Q1S  ] = chol( Q1 );
%     Q1sqrt = -Q1sqrt;% * Q1S;

%     P0sqrt = sqrt(sC)*I;
%     P1sqrt = sqrt(sT)*I;
    
%     Q0sqrt = blkdiag(-P0sqrt,sparse(m,m));
%     Q1sqrt = blkdiag(-P1sqrt,sparse(m,m));
    
%     Q(ig_ind,:) = 0;
%     Q0(ig_ind,:) = [];
%     Q1(ig_ind,:) = [];

    form = 'lp-vis';
    switch form
        case 'lp'
        o = [obj; zeros(m,1)];
        o(ignored_edges) = 0;
        B = [sparse(size(BT,1),m) BT];
        cvx_begin
            variable x(2*m);
            minimize( o' * x );
            subject to
                A * x <= b;
                B * x == ff(:);
        cvx_end
        ar = x(1:m);
        a = ar;
        n = reshape(x(m+(1:m*nf)),m,nf);
        
        case 'lp-vis'
        o = [obj; zeros(m,1)];
        ov = [objv; zeros(m,1)];
        p = [pva*100; zeros(m,1)];
        o(ignored_edges) = 0;
        p(ignored_edges) = 0;
        B = [sparse(size(BT,1),m) BT];
        cvx_begin
            variable x(2*m);
            minimize(ov' * x);% + p' * x );%+ norm(x,1));
            subject to
                A * x <= b;
                B * x == ff(:);
        cvx_end
        ar = x(1:m);
        a = ar;
        n = reshape(x(m+(1:m*nf)),m,nf);
        
%         case 'socp'
%         l = [zeros(m,1); obj; zeros(m,1)];
%         q = [pva*100; zeros(2*m,1)];
%         B = [sparse(size(BT,1),2*m) BT];
%         F = [pi*I,Z,Z;Z,Z,Z;Z,Z,Z];
%         G = sqrt(F);
%         c = [zeros(m,1);ones(m,1);zeros(m,1)];
%         lb = [zeros(2*m,1);-inf(m,1)];
%         ub = [inf(3*m,1)];
%         c1 = [ones(m,1);zeros(2*m,1)];
% %         q(ignored_edges) = 0;
%         cvx_begin
%             variable x(3*m);
%             minimize(l' * x + q' * x );%+ norm(x,1));
%             subject to
%                 As * x <= bs;
%                 B * x == ff(:);
% %                 lb <= x <= ub;
% %                 norm(G*x) <= c'*x;
%                 quad_form(x,F) - c'*x <= 0;
%         cvx_end
%         r = x(1:m);
%         a = x(m+1:2*m);
%         n = reshape(x(2*m+(1:m*nf)),m,nf);
%         
%         case 'qcqp'
%         q = [pva; zeros(m,1)];
% %         q(ignored_edges) = 0;
%         cvx_begin
%             variable x(2*m);
%             minimize( quad_form(x,Q) );% + q'*x);
%             subject to
%                 quad_form(x,Q0) + q0'*x >= 0;
%                 quad_form(x,Q1) + q1'*x >= 0;
% %                 norm( Q0sqrt*x ) + q0*x <= 0;
% %                 norm( Q1sqrt*x ) + q1*x <= 0;
%                 B * x == ff(:);
%         cvx_end
%         ar = x(1:m);
%         a = ar.^2*pi;
%         n = reshape(x(m+(1:m*nf)),m,nf);
    end
    


%%%%%%%%%%%
%%

% find areas bigger than a certain threshold and place a cylinder there
NZ = find(max(a,0)>1e-7);
num_rods = size(NZ,1)
num_compression = sum(sign(n(NZ))==1)
num_tension = sum(sign(n(NZ))==-1)
[CV,CF,CJ,CI] = edge_cylinders(V,E(NZ,:),...
    'PolySize',10,'Thickness',sqrt(max(a(NZ),0)/pi));

% obj = lengths .* v';
% f_value = obj.* a;
% % [a(NZ) n(NZ) v(NZ)' f_value(NZ)]

% plot the cylinders
color = sign(n(NZ(CJ)));
tsurf(CF,CV,falpha(1,0),'CData',color,fsoft);

%%
% % this is the bad one
% % bar visibilities
% new_h = new_h ./ segments';
% new_v = (max(new_h)-new_h);
% % get areas, axial forces, lengths, visibililities, nodal equilibrium matrix
% [new_a,new_n,BT] = groundstructure(V,E,f,bf,sC,sT,new_v,l,'IgnoredEdges',1:size(RV,1));
% 
% new_NZ = find(max(new_a,0)>1e-7);
% 
% new_obj = l .* new_v';
% new_f_value = new_obj.* new_a;
% [new_a(new_NZ) new_n(new_NZ) new_v(new_NZ)' new_f_value(new_NZ)]


%%
hold on
[~, I] = sortrows(E);

% plot(graph(E(:,1),E(:,2)),...
%     'XData',V(:,1),'YData',V(:,2),'ZData',V(:,3),...
%     'EdgeCData',theta(I), 'LineWidth',2);

% plot(graph(E(:,1),E(:,2)),...
%     'XData',V(:,1),'YData',V(:,2),'ZData',V(:,3),...
%     'EdgeCData',lengths, 'LineWidth',2);

% caxis([min(lengths) max(lengths)])
% colormap(flipud(cbrewer('RdBu',256)))
% colorbar

quiver3( ...
  V(:,1),V(:,2),V(:,3), ...
  f(:,1),f(:,2),f(:,3),'r','LineWidth',2);
% 
% scatter3(GV(:,1),GV(:,2),GV(:,3),...
%      '.','CData', Vs,'SizeData',30);


%%
hold off;
axis equal;
% view(26,15);
camlight;
camup([0 1 0]);
campos([-33 6 20])
% caxis([min(o), max(o)])
cameratoolbar('SetCoordSys','y')
cameratoolbar('setmode','orbit')
camproj('perspective')
caxis([-1 1])
colormap(flipud(cbrewer('RdBu',256)))
title(sprintf('Vol: %g, #bars: %d out of %i edges',lengths'*a,numel(NZ),size(E,1)),'FontSize',20);
colorbar

