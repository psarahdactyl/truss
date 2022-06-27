function [V,E,f,bf,ig] = construct_ground_structure(AV,AF,ACV,ACF,coms,force)
  % CONSTRUCT_GROUND_STRUCTURE
  % 
  %
  % Inputs:
  %   AV #AV by 3 list of vertices of each object in scene
  %   coms #length(objs) by 3 list of each object's center of mass
  %   force 1 by dim vector of force to place on each com

  % Outputs:
  %   V  #V by 3 list of ground structure vertices
  %   E  #E by 2 list of ground structure edges
  %   bf #bf by 1 list of fixed vertex indices into V
  %   ig #com*#AV by 1 list of ignored edges indices into E
  
  
  % create ground structure edges
  % loop over each component gather samples+normals
  XX = [];
  XC = [];
  n = 5;
%   for ci = 1:max(ACV)
%     [Vc,~,~,Fc] = remove_unreferenced(AV,AF(ACF==ci,:));
%     [Xc,~,~] = cy_blue_noise(n*sum(doublearea(Vc,Fc)),Vc,Fc,'Seed',rand);
%     XX = [XX;Xc];
%     XC = [XC;repmat(ci,size(Xc,1),1)];
%   end
  
  XX = AV;
  XC = ACV;

  % component adjacency
  C = ones(max(ACV),max(ACV))-eye(max(ACV));
  X2C = sparse(1:size(XX,1),XC,1,size(XX,1),max(ACV));

  % sample adjacency
  A = X2C*C*X2C';
  [AI,AJ] = find(triu(A));
  
  % edges
  E = [AI,AJ];
  
  % now add the center of masses to the list of vertices
  V = [XX;coms(2:end,:)];
  
  % forces
  f = zeros(size(V,1),3);
  f((size(XX,1)+1):end,:) = repmat(force,size(coms,1)-1,1);
  
  OE = E;
  
  [N,~,~] = histcounts(XC,'BinMethod','integers');
  for i=2:max(ACV)
    s = N(i);
    NE = [repelem(size(XX,1)+(i-1),s)' ((1:s)+sum(N(1:(i-1))))'];
    E = [E;NE];
  end

  % the first mesh will be the fixed one
  bf = find(XC==1);
  
  ig = (size(OE,1)+1):size(E,1);
  
end