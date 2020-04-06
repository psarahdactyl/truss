function [B,C] = create_nodal_equilibrium_matrices(V,E,bf)
  % CREATE_NODAL_EQUILIBRIUM_MATRICES 
  % 
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   E  #E by 2 list of edges indices into V
  %   bf #bf by 1 list of vertex indices which will be "fixed" (can
  %   withstand infinite force

  % Outputs:
  %   B  #V*dim by #E nodal equilibrium matrices for tension/compression
  %   C  #V*dim by #E nodal equilibrium matrices for bending
  %
  
  dim = size(V,2);
  n = size(V,1);
  m = size(E,1);
  
  EV = V(E(:,2),:)-V(E(:,1),:); % edge vectors for B matrix

  BV = zeros(size(EV));
  for i=1:size(E,1)
    b = null(EV(i,:));
    BV(i,:) = norm(EV(i,:))*b(:,1)'; % vectors for C matrix
  end
  
  B = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[EV;-EV],dim*n,m);
  C = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[BV;-BV],dim*n,m);

  % remove fixed vertices
  bf = bf(:);
  B(bf+n*(0:dim-1),:) = [];
  C(bf+n*(0:dim-1),:) = [];
  
  

end
