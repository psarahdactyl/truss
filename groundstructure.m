function [XX,XE,XC,YX,YE,YC] = groundstructure(VV,FF,CC,n)
  % [XX,XE,XC] = groundstructure(VV,FF,CC,n)
  %
  % Inputs:
  %   VV  #VV by 3 list of mesh vertex positions
  %   FF  #FF by 3 list of triangle indices into VV
  %   CC  #FF list of mesh "component"/"object" indices
  %   n   number of samples on the largest surface area object
  %       or 
  %       max(CC) list of samples per component
  % Outputs:
  %   XX  #XX by 3 list of ground structure positions
  %   XE  #XE by 2 list of indices into XX
  %   XC  #XC list of mesh component indices
  % 

  % loop over each component gather samples+normals
  XX = [];
  XN = [];
  XC = [];
  nc = max(CC);
  if numel(n) == 1
    areas = cell2mat(arrayfun(@(ci) ...
      sum(doublearea(VV,FF(CC==ci,:))),(1:nc)','UniformOutput',false));
    max_area = max(areas);
    nn = n*areas/max_area;
  else
    nn = n;
  end
  assert(numel(nn) == nc);

  % blue noise on each object for endpoints
  for ci = 1:nc
    [Vc,~,~,Fc] = remove_unreferenced(VV,FF(CC==ci,:));
    Nc = normalizerow(normals(Vc,Fc));
    [Xc,Ic,Bc] = blue_noise(nn(ci),Vc,Fc,'Seed',rand);
    Nc = Nc(Ic,:);
    XX = [XX;Xc];
    XN = [XN;Nc];
    XC = [XC;repmat(ci,size(Xc,1),1)];
  end
  
  % intercomponent all-pairs 
  % component adjacency
  AC = ones(nc,nc)-eye(nc);
  X2C = sparse(1:size(XX,1),XC,1,size(XX,1),nc);
  % sample adjacency
  A = X2C*AC*X2C';
  [AI,AJ] = find(triu(A));
  % edge
  XE = [AI,AJ];
  % edge vector
  XEV = XX(XE(:,2),:)-XX(XE(:,1),:);
  XEU = normalizerow(XEV);

  % unpruned groundstructure
  YX = XX;
  YE = XE;
  YC = XC;
  
  % discard if the angle w.r.t. to endpoint normal is too large
  max_angle = 80/180*pi;
  valid = ...
    sum(XEU.*XN(XE(:,1),:),2) > cos(max_angle) & ...
    sum(-XEU.*XN(XE(:,2),:),2) > cos(max_angle);
  XE = XE(valid,:);
  XEV = XEV(valid,:);
  XEU = XEU(valid,:);
  
  % discard if intersecting something else in the scene
  [H,T] = ray_mesh_intersect(XX(XE(:,1),:),XEV,VV,FF);
  H = H & T<0.999999;
  XE =   XE(~H,:);
  XEV = XEV(~H,:);
  XEU = XEU(~H,:);
  fprintf('K: %d\n',max(CC)-1);
  fprintf('YE: %d\n',size(YE,1));
  fprintf('XE: %d\n',size(XE,1));
end
