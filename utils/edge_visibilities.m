function [EVs] = edge_visibilities(V, E, GV, side, w, Vs, l)
  % EDGE_VISIBILITIES
  % 
  %
  % Inputs:
  %   the groundstructure:
  %   V    #V by 3 matrix of vertex positions of the union of every object in the scene
  %   E    #E by 2 matrix of edge indices into V
  %   GV   #GV by 3 list of voxel bottom left grid cell corners
  %   side 3 vector of voxel grid dimensions
  %   w    scalar representing the width of a voxel
  %   Vs   #GV by 1 list of visibility values in [0,1] for each voxel
  %   l    #E by 1 list of edge lengths

  % Outputs:
  %   EVs  #E by 1 list of visibility values for each edge in the ground
  %     structure


  X = reshape(GV(:,1),side([2 1 3]));
  Y = reshape(GV(:,2),side([2 1 3]));
  Z = reshape(GV(:,3),side([2 1 3]));

  origins = V(E(:,1),:);
  dests = V(E(:,2),:);
  % number of segments for each edge depending on lengths
  segments = ceil(l/(w(1)));

  % edge start and end points
  rep_origins = repelem(origins,max(segments),1);
  rep_dests = repelem(dests,max(segments),1);

  rep_slopes = rep_dests - rep_origins;

  % t values along edge
  ts = zeros(size(segments,1),max(segments));

  % yuck (not vectorized but idk how)
  for i=1:size(segments,1)
      p = linspace(0,1,segments(i));
      ts(i,1:numel(p)) = p;
  end

  tss = reshape(ts',size(ts,1)*size(ts,2),1);
  sts = rep_origins + rep_slopes.*tss;
  scs = interp3(X,Y,Z,reshape(Vs,side([2 1 3])),sts(:,1),sts(:,2),sts(:,3));

  EVs = sum(reshape(1-scs,max(segments),size(E,1)));

end