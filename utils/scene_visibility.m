function [Vs,GV,side,w] = scene_visibility(V,F,OV,OF,views)
  % SCENE_VISIBILITY
  % 
  %
  % Inputs:
  %   (m = #E)
  %   V  #V by 3 matrix of vertex positions of the union of every object in the scene
  %   F  #F by 3 matrix of triangle indices into V
  %   views #views by 3 matrix of viewpoints
  

  % Outputs:
  %   GV   #GV by 3 list of voxel grid cell corners
  %   side 3 vector of voxel grid dimensions
  %   w    scalar representing the width of a voxel
  %   Vs   #GV by 1 list of visibility values in [0,1] for each voxel

  
  % create voxel grid
  [GV,side,w] = voxel_grid(V,20);
  
  % calculate visibilities
  src = repmat(views,size(GV,1),1);
%   size(src)
  dir = repelem(GV,size(views,1),1);
%   size(dir)
  [Vs,~,~] = ray_mesh_intersect(src,dir-src,OV,OF);
%   size(Vs)
  Vs(Vs ~= 0) = 1;
  Vs = sum(reshape(Vs,size(views,1),size(GV,1)),1);
%   size(Vs)
  Vs = Vs ./ size(views,1);
  
end
  