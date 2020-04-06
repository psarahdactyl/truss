function [VV,EE] = prune_edges(V,E,varargin)
  % PRUNE_EDGES Given a network of edges, split long edges that geometrically
  % pass _through_ (but combinatorially are not incident on) existing vertices.
  % Optionally, remesh the network to include vertices at intersections between
  % segments (not at original vertices). This is a non-robust form of
  % snap-rounding. In 2D, one should probably call snap_rounding
  %
  % [VV,EE] = prune_edges(V,E)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   E  #E by 2 list of edges indices into V
  %   Optional:
  %     'SplitIntersections' followed by whether to split edges at intersections
  %       (by adding new vertices) {false}
  % Outputs:
  %   VV  #VV by dim list of vertex positions
  %   EE  #EE by 2 list of edges indices into VV
  %
  % See also:
  %   snap_rounding
  %

  split_intersections = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'SplitIntersections'}, ...
    {'split_intersections'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end


[EB1,EB2] = box_each_element(V,E);
I = box_intersect(EB1,EB2);
% distance between each potential match
[D,~,~,~,T] = segment_segment_squared_distance( ...
  V(E(I(:,1),2),:), ...
  V(E(I(:,1),1),:), ...
  V(E(I(:,2),2),:), ...
  V(E(I(:,2),1),:));
% remove false positives
I = I(D<1e-7,:);
T = T(D<1e-7,:);

% aggregate list of splits
knots = repmat({[0;1]},size(E,1),1);

% points of intersection
if split_intersections
  C = 0.5*( ...
    (T(:,1)).*V(E(I(:,1),2),:)+(1-T(:,1)).*V(E(I(:,1),1),:) + ...
    (T(:,2)).*V(E(I(:,2),2),:)+(1-T(:,2)).*V(E(I(:,2),1),:));
  % whether edges share a vertex
  EV = normalizerow(V(E(:,2),:)-V(E(:,1),:));
  phi = acos(sum(EV(I(:,1),:).*EV(I(:,2),:),2));
  parallel = phi<1e-4;
  both_at_endpoints = all((abs(T - 0.5)-0.5).^2<1e-7,2);
  A = facet_adjacency_matrix(E);
  crossing = ~parallel&~both_at_endpoints&~full(A(sub2ind(size(A),I(:,1),I(:,2))));
  % edges that cross each other without sharing a vertex: to resolve will need to
  % add one new vertex and split both edges
  cI = I(crossing,:);
  cT = T(crossing,:);
  for ci = 1:size(cI,1)
    for p = 1:2
      knots{cI(ci,p)}(end+1) = cT(ci,p);
    end
  end
end

% snap existing vertices at edges
I = box_intersect(EB1,EB2,V,V);
P = V(I(:,2),:);
A = V(E(I(:,1),1),:);
B = V(E(I(:,1),2),:);
UN = (B-A);
N = normalizerow(UN);
T = min(max(sum((P-A).*N,2)./normrow(UN),0),1);
C = (A+T.*UN);
D = normrow(P-C);
C = C(D<1e-7,:);
pT = T(D<1e-7,:);
pI = I(D<1e-7,1);
for pi = 1:size(pI,1)
  knots{pI(pi)}(end+1) = pT(pi);
end

reps = 1e-7;
% remove duplicate knots
knots = cellfun(@(l) sort(unique(round(l/reps))*reps),knots,'UniformOutput',false);

%VV = [V];
%EE = [E];
%for ei = 1:size(E,1)
%  nk = numel(knots{ei});
%  if nk > 2
%    % minus 2 for endpoints that already exist
%    Ei = [1:nk-1;2:nk]'-2;
%    Ei(1) = E(ei,1);
%    Ei(end,2) = E(ei,2);
%    EE(ei,:) = Ei(1,:);
%    EE = [EE;size(VV,1)+Ei(2:end,:)];
%    mids = knots{ei}(2:end-1);
%    VV = [VV;(mids).*V(E(ei,2),:)+(1-mids).*V(E(ei,1),:)];
%  end
%end
% fuck it. we'll do it live.
VV = cell2mat(arrayfun(@(ei) (knots{ei}).*V(E(ei,2),:)+(1-knots{ei}).*V(E(ei,1),:),1:size(E,1),'UniformOutput',false)');
lens = cellfun(@(l) numel(l),knots);
cs = [0;cumsum(lens)];
EE = cell2mat(arrayfun(@(ei) cs(ei)+[1:lens(ei)-1;2:lens(ei)]',1:size(E,1),'UniformOutput',false)');
[VV,~,J] = remove_duplicate_vertices(VV,1e-5);
EE = J(EE);
EE = unique(sort(EE,2),'rows');

end
