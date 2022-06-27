function [V,E,VC] = construct_ground_structure(AV,AF,ACV,ACF)
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
    
%   XX = [AV(ACV==1,:)];
%   XC = [ACV(ACV==1)];

  XX = [];
  XC = [];

  n = 100;
  tic % hacky way to randomize the blue noise sampling
  for ci = 1:max(ACV)
    [Vc,~,~,Fc] = remove_unreferenced(AV,AF(ACF==ci,:));
    if(n*sum(doublearea(Vc,Fc))<6)
      num=6;
    else
      num=n*sum(doublearea(Vc,Fc));
    end
    [Xc,~,~] = cy_blue_noise(num,Vc,Fc,'Seed',toc);
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
  
%   [V,E] =  prune_edges(V,E);
    
end
