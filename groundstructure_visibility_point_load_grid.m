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
% % plot object
% tsurf(RF,RV,...
%     'FaceVertexCData',MC,falpha(0.8,0),fsoft);

% sample points on surfaces
Rc = RV;
Dc = DV;
scatter3(Rc(:,1),Rc(:,2),Rc(:,3),'.g','SizeData',100);
scatter3(Dc(:,1),Dc(:,2),Dc(:,3),'.y','SizeData',100);

% this creates the ground structure, (V,E), including the center of mass
[I,J] = find(ones(size(Rc,1),size(Dc,1)));

% find centroid of object on which to apply force (gravity)
cm = centroid(RV,RF);
V = [Rc;Dc];
E = [I size(Rc,1)+J];

% nodal forces, boundary vertices, stress limits
f = zeros(size(V));
bf = size(Rc,1)+(1:size(Dc,1));
f(3,2) = -9.8;

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
    dim = size(V,2);
    nf = size(f,3);
  
    EV = normalizerow(V(E(:,2),:)-V(E(:,1),:));

    obj = lengths;

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
    
    As = [Z, repmat(-sC*I,nf,1) , -II; ...
        Z, repmat(-sT*I,nf,1) ,  II];
    bs = zeros(2*m*nf,1);
    
    form = 'socp';
    switch form
        case 'lp'
        o = [obj; zeros(m,1)];
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
        p = [pva*100; zeros(m,1)];
        B = [sparse(size(BT,1),m) BT];
        cvx_begin
            variable x(2*m);
            minimize(o' * x + p' * x );%+ norm(x,1));
            subject to
                A * x <= b;
                B * x == ff(:);
        cvx_end
        ar = x(1:m);
        a = ar;
        n = reshape(x(m+(1:m*nf)),m,nf);
        
        case 'socp'
        l = [zeros(m,1); obj; zeros(m,1)];
        q = [pva*100; zeros(2*m,1)];
        B = [sparse(size(BT,1),2*m) BT];
        F = [pi*I,Z,Z;Z,Z,Z;Z,Z,Z];
        G = sqrt(F);
        c = [zeros(m,1);ones(m,1);zeros(m,1)];

        cvx_begin
            variables x(3*m);
            minimize(l' * x + q' * x );%+ norm(x,1));
            subject to
                As * x <= bs;
                B * x == ff(:);
                {c'*x,1,G*x} <In> rotated_lorentz(3*m+2);
        cvx_end
        r = x(1:m);
        a = x(m+1:2*m);
        n = reshape(x(2*m+(1:m*nf)),m,nf);
%         l = obj;
%         p = pva*100000;
%         B = [sparse(size(BT,1),m) BT];
%         F = pi*I;
%         G = sqrt(F);
%         c = ones(m,1);
%         
%         cvx_begin
%             variables r(m) a(m) n(m);
%             minimize(l' * a + p' * r );%+ norm(x,1));
%             subject to
%                 A * [a;n] <= b;
%                 B * [a;n] == ff(:);
%                 {c'*a,1,G*r} <In> rotated_lorentz(m+2);
%         cvx_end
        
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

