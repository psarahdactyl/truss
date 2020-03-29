clear
close all
clf
hold on

% setup
% V = [0 0 0;
%     0.57735      0.57735      0.57735];

% V = [0 0 0;
%     1 0 0];

% V = [0 0 0;
%     1 -2 3];

E = [1 2];

dim = 3;
n = size(V,1);
m = size(E,1);

f = zeros(size(V,1),3);
bf = 1;

for i=1:10
f(2,2) = -0.98*i;
force = -0.98*i

% vector calculations
EV = V(E(:,2),:)-V(E(:,1),:); % edge vectors

[fr,fc] = find(f~=0);
fsum = sum(f(fr,:),1); % sum of forces

MV = cross(EV,repmat(fsum,size(E,1),1)); % moment vectors
MV = normalizerow(MV);

RV = cross(MV,EV); % rotation direction vector
% RV = normalizerow(RV);

l = edge_lengths(V,E); % length
r = 0.003175; % radius

BT = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[EV;-EV],dim*n,m);
CT = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[RV;-RV],dim*n,m);

% remove fixed vertices
bf = bf(:);
BT(bf+n*(0:dim-1),:) = [];
CT(bf+n*(0:dim-1),:) = [];
ff = f;
ff(bf,:,:) = [];

% sigma * a / l = n 
% sigma = (n * l) / a
% a = pi*r^2
B = full(BT);
ax = linsolve(B,ff');
sigma_tension = (ax * l) / (pi*r^2)

% sigma * a / l = m
C = full(CT);
be = linsolve(C,ff');
sigma_bending = (be * l) / (pi*r^2)

ratio = sigma_bending/sigma_tension
end

% plotting
[CV,CF,CJ,CI] = edge_cylinders(V,E, ...
    'PolySize',10,'Thickness',0.1);
  tsurf(CF,CV,falpha(1,0.05),'CData',1);
  
scatter3(V(bf,1),V(bf,2),V(bf,3),'.b','SizeData',1000);
  
quiver3(V(:,1),V(:,2),V(:,3),f(:,1),f(:,2),f(:,3),...
    0.5,'r','LineWidth',3);

p = V(2,:);
quiver3(p(1),p(2),p(3),MV(:,1),MV(:,2),MV(:,3),...
    0.5,'b','LineWidth',3);

quiver3(p(1),p(2),p(3),RV(:,1),RV(:,2),RV(:,3),...
    0.5,'c','LineWidth',3);


view(-90.351, 67.466)
axis equal
colormap(flipud(cbrewer('RdBu',256)))
caxis([-1 1])
camup([0 1 0])
cameratoolbar('SetCoordSys','y') 
cameratoolbar('setmode','orbit')
camproj('perspective');