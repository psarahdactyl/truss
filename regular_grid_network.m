function [V,E] = regular_grid_network(N,h)
% [V,E] = regular_grid_network(N,h)
  dim = numel(N);
  switch dim
  case 2
  case 3
    nx = N(1);
    ny = N(2);
    nz = N(3);
    [X,Y,Z] = meshgrid(1:nx,1:ny,1:nz);
    V = [X(:) Y(:) Z(:)];
    E = [1 2];
  end
  I = find(normrow(V(2:end,:)-V(1,:))<h)+1;
  E = [ones(numel(I),1) I];
  EV = V(E(:,2),:)-V(E(:,1),:);
  l = normrow(EV);
  EV = normalizerow(EV);
  [l,I] = sort(l,'ascend');
  E = E(I,:);
  EV = EV(I,:);
  A = sum(permute(EV,[1 3 2]).*permute(EV,[3 1 2]),3)>1-1e-7;
  [~,C] = conncomp(A);
  [~,I] = max(sparse(C,1:size(E,1),1,max(C),size(E,1)),[],2);
  l = l(I,:);
  E = E(I,:);
  EV = EV(I,:);
  %plot_groundstructure(V,E,0.01*1./l,l)
  %error

  sizes = N([2 1 3]);
  [eI,eJ,eK] = ind2sub(sizes,E(:,2));
  eI = eI-1;
  eJ = eJ-1;
  eK = eK-1;
  vec = @(X) X(:);
  eI = vec(eI.*[1 -1  1  1]);
  eJ = vec(eJ.*[1 -1 -1 -1]);
  eK = vec(eK.*[1  1  1 -1]);
  [eIJK] = unique([eI eJ eK],'rows');
  eI = eIJK(:,1);
  eJ = eIJK(:,2);
  eK = eIJK(:,3);
  [I,J,K] = ind2sub(sizes,1:size(V,1));
  dI = I+eI;
  dJ = J+eJ;
  dK = K+eK;
  sI = I+0*eI;
  sJ = J+0*eJ;
  sK = K+0*eK;
  good = ...
    (dI>=1 & dI<=sizes(1)) & ...
    (dJ>=1 & dJ<=sizes(2)) & ...
    (dK>=1 & dK<=sizes(3));
  dI = dI(good);
  dJ = dJ(good);
  dK = dK(good);
  sI = sI(good);
  sJ = sJ(good);
  sK = sK(good);
  E = [sub2ind(sizes,sI,sJ,sK) sub2ind(sizes,dI,dJ,dK)];
  E = unique(sort(E,2),'rows');


  %clf;
  %hold on;
  %scatter3(V(:,1),V(:,2),V(:,3));
  %%plot_edges(V,E(any(E==63,2),:),'-','LineWidth',3);
  %plot_edges(V,E,'-','LineWidth',3);
  %hold off;
  %view(-13,28);
  %axis equal;
end
