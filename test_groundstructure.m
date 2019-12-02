% "GRAND3 — Ground structure based topology optimization"
dim = 2;
switch dim
case 2
  nx =   7+0;
  ny = 2*2+1;
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
case 3
  tic;
  nx = 8+0;
  ny = 2*2+1;
  nz = 2*2+1;
  [V,E] = regular_grid_network([nx,ny,nz],3.7);
  V = (V-[1 (ny+1)/2 (nz+1)/2])/(nx-1);
  %[V,F] = regular_tetrahedral_mesh(nx,ny,nz);
  %V = (V.*([nx ny nz]-1)-[0 ny-1 nz-1]/2)./(nx-1);
  %% Create over connected graph
  %I = rangesearch(V,V,max(max(edge_lengths(V,F)))*2.1);
  %E = cell2mat(arrayfun( ...
  %  @(a) [repmat(a,1,numel(I{a}));I{a}],1:numel(I),'UniformOutput',false))';
  %E = E(E(:,1)<E(:,2),:);
  %[V,E] = prune_edges(V,E);
  toc
end

EV = normalizerow(V(E(:,2),:)-V(E(:,1),:));

switch dim
case 2
  % loaded vertices
  bl = snap_points([1 0],V);
  % smaller --> thicker
  sC = 300;
  sT = 300;
  % external forces
  f = zeros(size(V));
  f(bl,end,1) = -9.8;
  f(bl,end,2) =  9.8;
case 3
  % loaded vertices
  bl = snap_points([1 0 0],V);
  %bl = find(normrow(V(:,[1 3])-[1 0])<1e-5);
  % smaller --> thicker
  sC = 2000;
  sT = 2000;
  % external forces
  f = zeros(size(V));
  f(bl,end) = -9.8;
  %f(bl,2) = -1;

  %f = zeros([size(V),2]);
  %f(bl,3,1) = -9.8;
  %f(bl,3,2) =  9.8;

  f = zeros([size(V),4]);
  f(bl,3,1) = -9.8;
  f(bl,3,2) =  9.8;
  f(bl,2,3) = -9.8;
  f(bl,2,4) =  9.8;
  f = cat(3, ...
    f, ...
    f(:,:,1)*axisangle2matrix([1 0 0],pi/4), ...
    f(:,:,2)*axisangle2matrix([1 0 0],pi/4), ...
    f(:,:,3)*axisangle2matrix([1 0 0],pi/4), ...
    f(:,:,4)*axisangle2matrix([1 0 0],pi/4));
    
end
% fixed vertices
bf = find(V(:,1)==0);


tic
method = 'sum';
method = 'max';
method = 'columns';
switch method
case 'sum'
  [a,n,l] = groundstructure(V,E,sum(f,3),bf,sC,sT);
case 'max'
  a = zeros(size(E,1),1);
  n = zeros(size(E,1),size(f,3));
  for i = 1:size(f,3)
    [ai,n(:,i),l] = groundstructure(V,E,f(:,:,i),bf,sC,sT);
    a = max(a,ai);
  end
case 'columns'
  [a,n,l] = groundstructure(V,E,f,bf,sC,sT);
end
toc

NZ = find(max(a,0)>1e-7);

clf;
switch dim
case 2
  
  hold on;
  plot_groundstructure(V,E,a,n(:,1));
  quiver(V(:,1),V(:,2),f(:,1),f(:,2),'r','LineWidth',2);
  scatter(V(bl,1),V(bl,2),'.r','SizeData',1000);
  scatter(V(bf,1),V(bf,2),'.b','SizeData',1000);
  hold off;
case 3
  [CV,CF,CJ,CI] = edge_cylinders(V,E(NZ,:), ...
    'PolySize',10,'Thickness',sqrt(max(a(NZ),0)/pi));
  hold on;
  nf = size(f,3);
  quiver3( ...
    repmat(V(:,1),[1,1,nf]),repmat(V(:,2),[1,1,nf]),repmat(V(:,3),[1,1,nf]), ...
    f(:,1,:),f(:,2,:),f(:,3,:),4,'r','LineWidth',2);
  scatter3(V(bl,1),V(bl,2),V(bl,3),'.r','SizeData',1000);
  scatter3(V(bf,1),V(bf,2),V(bf,3),'.b','SizeData',1000);
  tsurf(CF,CV,falpha(1,0.05),'CData',n(NZ(CJ),1));
  hold off;
  view(19,28)
  drawnow;
  colormap(flipud(cbrewer('RdBu',256)))
  caxis(max(abs(caxis))*[-1 1])
end

axis equal
title(sprintf('Vol: %g, #bars: %d',l'*a,numel(NZ)),'FontSize',20);

  

