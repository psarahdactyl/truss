clear
close all
clf
hold on

% setup
% bowling ball weights
% 6 lb = 2.72 kg
% 7 lb = 3.18 kg
% 8 lb = 3.63 kg
% 9 lb = 4.08 kg
% 10 lb = 4.54 kg
% 12 lb = 5.44 kg
% 14 lb = 6.35 kg
% 16 lb = 7.26 kg


V = [0 0 0;
    1 0 0];

E = [1 2];

dim = 3;
n = size(V,1);
m = size(E,1);

f = zeros(size(V,1),3);
bf = 1;

g = [0 -9.8 0];
weight = 5.44;

% % 100 random 3D vectors with values between -10 and 10
% forces = -10 + (10-(-10)).*rand(1000,3);

f(2,:) = weight*g;

% vector calculations
EV = V(E(:,2),:)-V(E(:,1),:); % edge vectors

[fr,fc] = find(f~=0);
fsum = sum(f(fr,:),1); % sum of forces

MV = cross(EV,repmat(fsum,size(E,1),1)); % moment vectors
MV = normalizerow(MV);

RV = cross(MV,EV); % rotation direction vector
% RV = normalizerow(RV);

l = edge_lengths(V,E); % length = 1 meter
r = 0.003175; % radius
% diameter is 0.25 inch = 0.635 cm

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
area = (pi*r^2)

B = full(BT);
ax = linsolve(B,ff');
sigma_tension = (ax * l) / (pi*r^2)

% sigma * a / l = m
C = full(CT);
be = linsolve(C,ff');
sigma_bending = (be * l) / (pi*r^2)

ratio = sigma_bending/sigma_tension


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

% tiledlayout(3,1)
% nexttile
% histogram(vecnorm(forces,2,2),200);
% title('Force Magnitudes','FontSize',20)
% 
% nexttile
% histogram(s_tensions,200);
% title('Tension Stresses','FontSize',20)
% 
% nexttile
% histogram(s_bendings,200);
% title('Bending Stresses','FontSize',20)

% scatter3(vecnorm(forces,2,2),abs(s_tensions),abs(s_bendings),...
%         50,abs(ratios),'filled')
% xlabel('Force Magnitudes') 
% ylabel('Tension Stress') 
% zlabel('Bending Stress')

view(-90.351, 67.466)
axis equal
datacursormode on
colormap(flipud(cbrewer('RdBu',256)))
caxis([-1 1])
camup([0 1 0])
cameratoolbar('SetCoordSys','y') 
cameratoolbar('setmode','orbit')
camproj('perspective');