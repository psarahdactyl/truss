% starting in 2D

%                       |
%                       v
% 3. --- 6. --- 9. --- 12. 
%  |      |      |      |
%  |      |      |      |
% 2. --- 5. --- 8. --- 11.
%  |      |      |      |
%  |      |      |      |
% 1. --- 4. --- 7. --- 10.

xRes = 7; yRes = 5; dim=2;

% yield stresses (tension and compression) in pascals, i guess
sigma_t = 4.6%e7;
sigma_c = 6.5%e7;

Ndof = xRes*yRes*dim;

xSpace = linspace(1,xRes,xRes);
ySpace = linspace(1,yRes,yRes);
[X, Y] = meshgrid(xSpace, ySpace);
UV = [X(:), Y(:)];
unit_vecs = eye(dim);

E = nchoosek(1:size(UV,1),dim);
l = vecnorm(UV(E(:,1),:)-UV(E(:,2),:),2,2);

% nodal equilibrium matrix (Ndof by Nb)
Nb = size(E,1);
B = zeros(Ndof,Nb);

for e = 1:Nb
    v1 = UV(E(e,1),:);
    v2 = UV(E(e,2),:);
    edge = v1 - v2;
    B((E(e,1)-1)*2+1,e) = dot(unit_vecs(1,:),edge/norm(edge));
    B((E(e,1)-1)*2+2,e) = dot(unit_vecs(2,:),edge/norm(edge));
    B((E(e,2)-1)*2+1,e) = -dot(unit_vecs(1,:),edge/norm(edge));
    B((E(e,2)-1)*2+2,e) = -dot(unit_vecs(2,:),edge/norm(edge));
end

ax1 = subplot(1,2,1)
plot_edges(UV,E,'LineWidth',3);

axis([1 xRes 1 yRes]);
grid on

% nodal forces (x and y) in newtons
vertex = snap_points([xRes ceil(yRes/2)],UV);
f = zeros(Ndof,1);
f((vertex-1)*2+2) = -9.8;

l = [l; zeros(Nb,1)];

A = [-diag(sigma_c*ones(Nb,1)) -diag(ones(Nb,1)); 
     -diag(sigma_t*ones(Nb,1))  diag(ones(Nb,1))];
b = zeros(Nb*2,1);

% fixed points
I = find(UV(:,1)==1);
I = [(I-1)*2+1; (I-1)*2+2];
B(I,:) = [];
Nfree = Ndof - size(I,1);

Aeq = [zeros(Nfree,Nb) B];
f(I,:) = [];
beq = f;

% % max area (a in [0,1]) in cm or m, i guess
% max_a = 5.;
% max_n = -9.8;
lb = [zeros(Nb,1); -inf*ones(Nb,1)];
% ub = [max_a*ones(Nb,1); zeros(Nb,1)];

% min_a V = l^T*a
[S,vol,exitflag,output] = linprog(l,A,b,Aeq,beq,lb,[])
Sr = S;
Sr = reshape(Sr,numel(Sr)/2,2);
areas = Sr(:,1);
axial_forces = Sr(:,2);

nE = E;
nE(areas == 0,:) = 0;
nE = nE(any(nE ~= 0,2),:);

ax2 = subplot(1,2,2);
% plot_edges(UV,nE,'LineWidth',3);
tsurf(nE,UV,'FaceIndices',1,'VertexIndices',1)
axis([1 xRes 1 yRes]);
grid on
