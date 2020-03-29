clear
close all
hold on

% rod stuff
V = [0 0 0;
     1 1 1];
E = [1 2];

% plot_edges(V,E);
scatter3(V(:,1),V(:,2),V(:,3),'.b','SizeData', 500);

% force on rod
F = zeros(size(V,1),3);
F(2,2) = -9.8;
% quiver3(V(:,1),V(:,2),V(:,3),F(:,1),F(:,2),F(:,3),norm(F),'r','LineWidth',3);
quiver3(V(:,1),V(:,2),V(:,3),F(:,1),F(:,2),F(:,3),0.5,'r','LineWidth',3);

r = 0.635; % cm
a = pi*r^2; % m^2

% moment of inertia
I = a*r^2/4; % m^4

% max bending moment
% PV = V(E(:,2),:)-V(E(:,1),:);
% PV = normalizerow(PV);
% PV = [-PV(:,2) PV(:,1) zeros(size(PV,1),1)]; 
PV = [zeros(size(E,1),1) ones(size(E,1),1) zeros(size(E,1),1)]; % C matrix
MV = (V(E(:,2),:)+V(E(:,1),:))/2; % midpoints for plotting

M = cross(V(E(:,2),:)-V(E(:,1),:),F(2,:));
M_max = norm(M); % m
% M_max = norm(rod)*norm(F);

% bending stress
sigma_b = M_max * r / I % N * m^2 / m^4 = N / m^2 = Pa

% quiver3(MV(:,1),MV(:,2),MV(:,3),M(:,1),M(:,2),M(:,3),dot([1 0 0],M),'g','LineWidth',3);
quiver3(MV(:,1),MV(:,2),MV(:,3),M(:,1),M(:,2),M(:,3),1/norm(M),'b','LineWidth',3);
quiver3(MV(:,1),MV(:,2),MV(:,3),PV(:,1),PV(:,2),PV(:,3),'g','LineWidth',3);

[CV,CF,CJ,CI] = edge_cylinders(V,E,...
    'PolySize',10,'Thickness',r);

% plot the cylinders
color = sign(M_max);
tsurf(CF,CV,falpha(1,0),'EdgeColor',[0.2 0.2 0.2])%,'CData',color);

axis equal
camup([0 1 0])
cameratoolbar('SetCoordSys','y')
cameratoolbar('setmode','orbit')
camproj('perspective');