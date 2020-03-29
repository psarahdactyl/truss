clear
close all
clf
hold on

nx = 4;
ny = 2;
nz = 3;

% % example 4
% [V,E] = regular_grid_network([nx,ny,nz],3);
% V = (V-[1 (ny+1)/2 (nz+1)/2])/(nx-1);
% f = zeros(size(V,1),3);
% f(16,2) = -9.8;

% % example 3
% V = [0 0 0;
%     0 1 0;
%     0 0 1;
%     1 0 0;
%     1 0 1;
%     1 1 0;
%     0 1 1;
%     1 1 1];
% f = zeros(size(V,1),3);
% f(8,2) = -9.8;

% % example 1
% V = [0 1 0;
%     1 1 0;
%     0 1 1;
%     1 1 1];
% f = zeros(size(V,1),3);
% f(4,2) = -9.8;

% example 2
V = [0 0 0;
    0 0 1;
    0 0 2;
    1 0 0;
    1 0 1;
    1 0 2];
f = zeros(size(V,1),3);
f(5,2) = -9.8;


E = [repmat((1:size(V,1))',size(V,1),1) repelem((1:size(V,1))',size(V,1),1)];
E = sort(E,2);
E = unique(E,'rows');
E(any(diff(E,[],2)==0,2),:)=[]; % get rid of rows of edges from a vertex to itself


bf = find(V(:,1) == 0);
 
% % plot boundary conditions
% plot_edges(V,E);
% scatter3(V(bf,1),V(bf,2),V(bf,3),'.b','SizeData',1000);
% 
% % plot force
% quiver3(V(:,1),V(:,2),V(:,3),f(:,1),f(:,2),f(:,3),0.5,'r','LineWidth',3);
% 
% 
% axis equal
% camup([0 1 0])
% cameratoolbar('SetCoordSys','y')
% cameratoolbar('setmode','orbit')
% camproj('perspective');
% 
% return

sC = 1e3;
sT = 1e3;
sB = 1e3;

dim=3;
n = size(V,1);
m = size(E,1);

EV = V(E(:,2),:)-V(E(:,1),:); % edge vectors
EV_normalized = normalizerow(EV); % B matrix

% PV = V(E(:,2),:)-V(E(:,1),:); % perpendicular vectors
% PV = normalizerow(PV);
% PV = [PV(:,2) -PV(:,1) zeros(size(PV,1),1)]; % C matrix

[fr,fc] = find(f~=0);
fsum = sum(f(fr,:),1); % sum of forces

MV = cross(EV,repmat(fsum,size(E,1),1)); % moment vectors
CV = cross(MV,EV);
CV_normalized = normalizerow(CV);
% CV = cross(V(E(:,2),:)-V(E(:,1),:),repmat([0 -9.8 0],size(E,1),1));

MP = (V(E(:,2),:)+V(E(:,1),:))/2; % midpoints for plotting

l = edge_lengths(V,E); % objective

BT = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[EV;-EV],dim*n,m);
CT = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[CV;-CV],dim*n,m);
BT_normalized = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',...
                [EV_normalized;-EV_normalized],dim*n,m);
CT_normalized = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',...
                [CV_normalized;-CV_normalized],dim*n,m);

% remove fixed vertices
bf = bf(:);
BT(bf+n*(0:dim-1),:) = [];
CT(bf+n*(0:dim-1),:) = [];
BT_normalized(bf+n*(0:dim-1),:) = [];
CT_normalized(bf+n*(0:dim-1),:) = [];
ff = f;
ff(bf,:,:) = [];

I = speye(m,m);
Z = sparse(m,m);
nf = size(f,3);

A = [repmat(-sC*diag(1./l),nf,1), -I, Z; ...
    repmat(-sT*diag(1./l),nf,1) , I, Z;...
    repmat(-sB*diag(1./l),nf,1), Z, -I;...
    repmat(-sB*diag(1./l),nf,1), Z, I];

A_normalized = [repmat(-sC*I,nf,1), -I, Z; ...
    repmat(-sT*I,nf,1) , I, Z;...
    repmat(-sB*I,nf,1), Z, -I;...
    repmat(-sB*I,nf,1), Z, I];
b = zeros(4*m*nf,1);


B = [sparse(size(BT,1),m) BT sparse(size(BT,1),m)];
C = [sparse(size(CT,1),2*m) CT];

B_normalized = [sparse(size(BT_normalized,1),m) BT_normalized sparse(size(BT_normalized,1),m)];
C_normalized = [sparse(size(CT_normalized,1),2*m) CT_normalized];
%%%%%%%%%%%%%% optimization %%%%%%%%%%%%%%%

%%%%%%% normalizing %%%%%%%%
[anb,~,flags,output] = linprog([l;zeros(2*m*nf,1)],A_normalized,b,...
  B_normalized+C_normalized,ff(:), ...
  [zeros(m,1);-inf(2*m,1)], ...
  [inf*ones(m,1);inf(2*m,1)]);
arnb = anb(1:m);
axnb = reshape(anb(m+(1:m*nf)),m,nf);
benb = reshape(anb(2*m+(1:m*nf)),m,nf);

%%%%%%% without normalizing (using 1/l) %%%%%%%%
[an,~,flags,output] = linprog([l;zeros(2*m*nf,1)],A,b,...
  B+C,ff(:), ...
  [zeros(m,1);-inf(2*m,1)], ...
  [inf*ones(m,1);inf(2*m,1)]);
ar = an(1:m);
ax = reshape(an(m+(1:m*nf)),m,nf);
be = reshape(an(2*m+(1:m*nf)),m,nf);


tiledlayout(2,1)
nexttile
title(sprintf('Normalizing Vectors, Vol: %g, #bars: %d',l'*arnb,sum(log10(arnb)>-7)),'FontSize',20);
hold on;
plot_groundstructure(V,E,arnb,benb);
plot_edges(V,E,'Color','#848484','LineWidth',0.001);
% plot boundary conditions
scatter3(V(bf,1),V(bf,2),V(bf,3),'.b','SizeData',1000);
% plot force
quiver3(V(:,1),V(:,2),V(:,3),f(:,1),f(:,2),f(:,3),0.5,'r','LineWidth',3);
% % plot moment
quiver3(MP(:,1),MP(:,2),MP(:,3),CV(:,1),CV(:,2),CV(:,3),0.5,'c','LineWidth',3);
axis equal
view(115.69,53.858)
colormap(flipud(cbrewer('RdBu',256)))
caxis([-1 1])
camup([0 1 0])
cameratoolbar('SetCoordSys','y')
cameratoolbar('setmode','orbit')
camproj('perspective');

nexttile
title(sprintf('Using (1/l) term, Vol: %g, #bars: %d',l'*ar,sum(max(ar,0)>1e-7)),'FontSize',20);
hold on;
plot_groundstructure(V,E,ar,be);
plot_edges(V,E,'Color','#848484','LineWidth',0.001);
% plot boundary conditions
scatter3(V(bf,1),V(bf,2),V(bf,3),'.b','SizeData',1000);
% plot force
quiver3(V(:,1),V(:,2),V(:,3),f(:,1),f(:,2),f(:,3),0.5,'r','LineWidth',3);
% % plot direction of rotation
% quiver3(MP(:,1),MP(:,2),MP(:,3),PV(:,1),PV(:,2),PV(:,3),0.5,'c','LineWidth',3);
% % plot moment
quiver3(MP(:,1),MP(:,2),MP(:,3),CV(:,1),CV(:,2),CV(:,3),0.5,'c','LineWidth',3);
% ax = gca
view(115.69,53.858)
% axis equal
colormap(flipud(cbrewer('RdBu',256)))
caxis([-1 1])
camup([0 1 0])
cameratoolbar('SetCoordSys','y')
cameratoolbar('setmode','orbit')
camproj('perspective');

% plot_edges(V,E,'--k','LineWidth',3);
% scatter(V(bf,1),V(bf,2),'.b','SizeData',3000);
% quiver(V(:,1),V(:,2),f(:,1),f(:,2),0.5,'r','LineWidth',3);
% quiver(MV(:,1),MV(:,2),PV(:,1),PV(:,2),0.5,'g','LineWidth',3);
% quiver(MV(:,1),MV(:,2),EV(:,1),EV(:,2),0.5,'c','LineWidth',3);
% quiver(V(1+(1:numel(aa)),1),V(1+(1:numel(aa)),2),fa(:,1),fa(:,2),0,'b','LineWidth',3);
hold off;
axis equal
 
