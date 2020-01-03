
nx =   8+0;
ny = 2*3+1;
[V,F] = create_regular_grid(nx,ny);
  V = (V.*([nx ny]-1)-[0 ny-1]/2)./(nx-1);
% Create over connected graph
I = rangesearch(V,V,max(max(edge_lengths(V,F)))*4.1);
E = cell2mat(arrayfun( ...
  @(a) [repmat(a,1,numel(I{a}));I{a}],1:numel(I),'UniformOutput',false))';
E = E(E(:,1)<E(:,2),:);
tic;
[V,E] = prune_edges(V,E);
toc
EV = normalizerow(V(E(:,2),:)-V(E(:,1),:));
% loaded vertices
bl = snap_points([1 0],V);
% smaller --> thicker
sC = 300;
sT = 300;
% external forces
f = zeros(size(V));
f(bl,end,1) = -9.8;
% fixed vertices
bf = find(V(:,1)==0);
[a,n,l] = optimize_groundstructure(V,E,f,bf,sC,sT);
NZ = find(max(a,0)>1e-7);
%clf;
%hold on;
%plot_groundstructure(V,E,a,n(:,1));
%quiver(V(:,1),V(:,2),f(:,1),f(:,2),'r','LineWidth',2);
%scatter(V(bl,1),V(bl,2),'.r','SizeData',1000);
%scatter(V(bf,1),V(bf,2),'.b','SizeData',1000);
%hold off;
%axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%

E = E(NZ,:);
a = a(NZ);
F = f;
f = n(NZ);


[V,I,J] = remove_unreferenced(V,E);
bf = intersect(I(bf),1:size(V,1));
F = F(J,:);
E = I(E);

% Augh... seems I need to remove vertices with parallel incoming/outgoing edges.
% Why do I have to do this? 
EV = normalizerow(V(E(:,2),:)-V(E(:,1),:));
A = sparse(E,repmat(1:size(E,1),2,1)',1,size(V,1),size(E,1));
bad = [];
for v = 1:size(V,1)
  % should be sorted but to be sure
  vei = sort(find(A(v,:)));
  if numel(vei) == 2
    if abs(sum(EV(vei(1),:).*EV(vei(2),:))) >= 1-1e-7

      E(vei(1),1+(E(vei(1),2)==v)) = E(vei(2),1+(E(vei(2),2)~=v));
      a(vei(1)) = 0.5*(a(vei(1))+ a(vei(2)));
      n(vei(1)) = n(vei(1))+ f(vei(2));
      E(vei(2),:) = v;
      A(v,vei(1)) = 0;
      A(:,vei(1)) = A(:,vei(1)) + A(:,vei(2));
      A(:,vei(2)) = 0;

    end
  end
end

keep = E(:,1)~=E(:,2);
E = E(keep,:);
a = a(keep);
%a(:) = max(a);
f = f(keep);

[V,I,J] = remove_unreferenced(V,E);
bf = intersect(I(bf),1:size(V,1));
F = F(J,:);
E = I(E);


EV = V(E(:,2),:)-V(E(:,1),:);
l = normrow(EV);
EV = EV./l;

n = size(V,1);
m = size(E,1);
dim = 2;
% Freund's notation
A = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[EV;-EV],dim*n,m)';
% "volume"
t = l.*a;
% self weight
F = sparse(repmat(E(:),1,2),repmat([1 2],numel(E),1),repmat(t,2,1).*[0 -9.8],size(V,1),size(V,2));

youngs = 10000;
B = diag(sparse( l.^2 ./ (t.* youngs)));
K = A' * (B\A);
b = reshape(bf+(0:dim-1)*n,[],1);
u = min_quad_with_fixed(0.5*K,-F(:),b,zeros(numel(b),1));
u = reshape(u,size(V));


clf;
hold on;
plot_groundstructure(V,E,a,B\(A*u(:)));
quiver(V(:,1),V(:,2),u(:,1),u(:,2),'LineWidth',3);
hold off;
axis equal;
colormap(cbrewer(flipud('RdBu'),256))
