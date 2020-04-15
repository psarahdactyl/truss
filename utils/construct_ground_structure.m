function [V,E,VC,bf] = construct_ground_structure(AV,AF,ACV,ACF)
  % CONSTRUCT_GROUND_STRUCTURE
  % 
  %
  % Inputs:
  %   AV  #AV by 3 list of vertices of each object in scene
  %   AF  #AF by 3 list of face indices into AV
  %   ACV #AV by 1 list of object indices into rows of AV
  %   ACF #AF by 1 list of object indices into rows of AF
  %   

  % Outputs:
  %   V  #V by 3 list of ground structure vertices
  %   E  #E by 2 list of ground structure edges
  %   VC #VC by 1 list of object indices into V, the ground structure
  
  
  % create ground structure edges
  % loop over each object and gather samples into V
  
  bf=1:size(find(ACV==1),1);
  
%   XX = [AV(ACV==1,:)];
%   XC = [ACV(ACV==1)];

  XX = [];
  XC = [];

  n = 7;
  for ci = 1:max(ACV)
    [Vc,~,~,Fc] = remove_unreferenced(AV,AF(ACF==ci,:));
    [Xc,~,~] = blue_noise(n*sum(doublearea(Vc,Fc)),Vc,Fc,'Seed',rand);
    XX = [XX;Xc];
    XC = [XC;repmat(ci,size(Xc,1),1)];
  end
  
%   XX = AV;
%   XC = ACV;

  % component adjacency
  C = ones(max(ACV),max(ACV))-eye(max(ACV));
  X2C = sparse(1:size(XX,1),XC,1,size(XX,1),max(ACV));

  % sample adjacency
  A = X2C*C*X2C';
  [AI,AJ] = find(triu(A));
  
  % edges
  E = [AI,AJ];
  
  V=XX;
  VC=XC;
    
end