function [EAs, d, theta] = edge_projected_visible_areas(V, E, views, EVs, func)
  % EDGE_PROJECTED_VISIBLE_AREAS
  % 
  %
  % Inputs:
  %   the groundstructure:
  %   V     #V by 3 matrix of vertex positions of the union of every object in the scene
  %   E     #E by 2 matrix of edge indices into V
  %   views #views by 3 matrix of viewpoints
  %   EVs   #E by 1 list of visibility values for each edge in the ground

  % Outputs:
  %   EAs  #E by 1 list of projected visible area values for each edge in 
  %     the ground structure

  
  % average distances from viewpoints to bar midpoints
  MP = (V(E(:,1),:) + V(E(:,2),:))/2;
  vs = repmat(views,size(E,1),1);
  ms = repelem(MP,size(views,1),1);
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

  EAs = theta.^2 .* func(EVs').^2 ./ d.^2;
  
end