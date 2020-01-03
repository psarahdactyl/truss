n = 3;
theta = linspace(0,2*pi,n+1);
theta = theta(1:end-1)';
BV = [cos(theta) sin(theta)];
BF = delaunay(BV);

[V,F] = extrude(BV,BF);
phi = (1/n+1/(2*n))*2*pi;
phi = (1/n+1/(2*n))*2*pi-0.1;
V(V(:,3)>0.5,:) = V(V(:,3)>0.5,:)*axisangle2matrix([0 0 1],phi);
r = 31.25;
h = 2*r;
V = V*diag([r r h]);
% let's use meters
V = V/100;
E = edges(F);
%E(4,:) = [];

% 1" diameter
a = ((0.0254*0.5)^2*pi)*ones(size(E,1),1);
EV = V(E(:,2),:)-V(E(:,1),:);
l = normrow(EV);
EV = EV./l;
n = size(V,1);
m = size(E,1);
dim = 3;
% Freund's notation
A = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[EV;-EV],dim*n,m)';
% "volume"
t = l.*a;
% self weight
F = sparse( ...
  repmat(E(:),1,dim), ...
  repmat(1:dim,numel(E),1), ...
  repmat(t,2,1).*[0 0 -9.8], ...
  n,dim);
youngs = 10000;
B = diag(sparse( l.^2 ./ (t.* youngs)));
K = A' * (B\A);

bf = find(V(:,3)==0);
%b = reshape(bf+(0:dim-1)*n,[],1);
b = 2*n+bf;

u = min_quad_with_fixed(0.5*K,-F(:),b,zeros(numel(b),1));
u = reshape(u,size(V));
f = B\(A*u(:));

clf;
hold on;
plot_groundstructure(V,E,a,f);
quiver3(V(:,1),V(:,2),V(:,3),u(:,1),u(:,2),u(:,3),0,'LineWidth',3);
hold off;
axis equal;
caxis(max(abs(caxis))*[-1 1])
colormap(cbrewer(flipud('RdBu'),256));
colorbar
view(15,15)

max(abs(u))
