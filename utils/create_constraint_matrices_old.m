function [A,b,Aeq,beq] = create_constraint_matrices_old(V,E,f,bf,sC,sT,sB,varargin)
  % CREATE_CONSTRAINT_MATRICES 
  % 
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   E  #E by 2 list of edges indices into V
  %   f  #V by dim list of forces on each vertex
  %   bf #bf by 1 list of vertex indices which will be "fixed" (can
  %     withstand infinite force
  %   sC,sT,sB scalars representing the compression, tension, bending yield
  %     stress, respectively

  % Outputs:
  %   (m = #E)
  %   A 4*m by 3*m inequality matrix (stress bounds)
  %   b 4*m by 1 vector of 0's
  %   Aeq #V*dim by 3*m equality matrix (force balance)
  %   beq #V*dim by 1 vectorized force matrix
  
  bending=1;
  
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Bending'}, ...
    {'bending'});
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
  
  n = size(V,1);
  m = size(E,1);
  nf = size(f,3);

  I = speye(m,m); % m x m identity matrix 
  Z = sparse(m,m); % m x m zero matrix
  
  l = edge_lengths(V,E); % lengths
  L = spdiags(1./l(:),0,m,m);

  A = [repmat(-sC*L,nf,1), -I,  Z;...
       repmat(-sT*L,nf,1),  I,  Z;...
       repmat(-sB*L,nf,1),  Z, -I;...
       repmat(-sB*L,nf,1),  Z,  I];

  b = zeros(4*m*nf,1);
  
  [BT,CT] = create_nodal_equilibrium_matrices_old(V,E,bf); % already here the fixed vertices have been removed
  f(bf,:,:) = []; % must remove relevant rows in the force vector too
  
  B = [sparse(size(BT,1),m) BT sparse(size(BT,1),m)];
  C = [sparse(size(CT,1),2*m) CT];
  
  if bending==1
    Aeq = B+C;
  else
    Aeq = B;  
  end
   
  beq = f(:);
  

end
