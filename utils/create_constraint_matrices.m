function [A,b,Aeq,beq] = create_constraint_matrices(V,E,ACV,coms,force,sC,sT,sB,varargin)
  % CREATE_CONSTRAINT_MATRICES 
  % 
  %
  % Inputs:
  %   V   #V by dim list of vertex positions
  %   E   #E by 2 list of edges indices into V
  %   ACV #AV by 1 list of object indices into rows of AV
  %   coms #max(ACV) by 3 list of each object's center of mass
  %   sC,sT,sB scalars representing the compression, tension, bending yield
  %     stress, respectively

  % Outputs:
  %   (m = #E)
  %   A 6*m by 3*m inequality matrix (stress bounds)
  %   b 6*m by 1 vector of 0's
  %   Aeq #k*dim by 3*m equality matrix (force and torque balance)
  %   beq #k*dim by 1 vectorized force matrix
  
  bending=1;
  
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Bending'},...
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
  
  bending
%   n = size(V,1);
  m = size(E,1);
  nf = 1;
  

  I = speye(m,m); % m x m identity matrix 
  Z = sparse(m,m); % m x m zero matrix
  
  l = edge_lengths(V,E); % lengths

  sCI = spdiags(sC(:),0,m,m);
  sTI = spdiags(sT(:),0,m,m);
  sBL = spdiags(sB.*l(:),0,m,m);
  L = spdiags(1./l(:),0,m,m);
  
  if(bending==0)
  A = [repmat(-sTI,nf,1), -I;...
       repmat(-sCI,nf,1),  I];
  b = zeros(2*m*nf,1);
  else
  A = [repmat(-sTI,nf,1), -I,  Z,  Z;...
       repmat(-sCI,nf,1),  I,  Z,  Z;...
       repmat(-sBL,nf,1),  Z, -I,  Z;...
       repmat(-sBL,nf,1),  Z,  I,  Z;...
       repmat(-sBL,nf,1),  Z,  Z, -I;...
       repmat(-sBL,nf,1),  Z,  Z,  I];
  b = zeros(6*m*nf,1);
%   A = [repmat(-sC.*I,nf,1), -I,  Z,  Z;...
%        repmat(-sT.*I,nf,1),  I,  Z,  Z;...
%        repmat(-sB.*L,nf,1),  Z, -I,  Z;...
%        repmat(-sB.*L,nf,1),  Z,  I,  Z;...
%        repmat(-sB.*L,nf,1),  Z,  Z, -I;...
%        repmat(-sB.*L,nf,1),  Z,  Z,  I];
%   b = zeros(6*m*nf,1);
  end
  
  [S,G,C,B] = create_nodal_equilibrium_matrices(V,E,ACV,coms);
  
  % S*C is 3*k x m
  % S*B is 3*k x 2*m
  
  if(bending==0)
    FE = [sparse(3*(max(ACV)-1),m) S*C];% S*B];
    TE = [sparse(3*(max(ACV)-1),m) G*C];% G*B];
  else
    FE = [sparse(3*(max(ACV)-1),m) S*C S*B];
    TE = [sparse(3*(max(ACV)-1),m) G*C G*B];
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % bending+tension/compression+torque
  
  Aeq = [FE;TE];
  beq = [repmat(force',max(ACV)-1,1);zeros(3*(max(ACV)-1),1)];

%   Aeq = [sparse(size(T,1)+size(F,1),m) [F;T]]; % 6 by m [0 t;f 0]*[a;n] = [0;-f]
%   beq = [repmat([0 -9.8 0]',max(ACV)-1,1);zeros(3*(max(ACV)-1),1)];
       
%   size(Aeq)
%   size(beq)
  
%   % weights
%   beq(1:3,:) = beq(1:3,:)*7.1;
%   beq(4:6,:) = beq(4:6,:)*0.9;
%   beq(7:9,:) = beq(7:9,:)*2.5;
  

%   if bending==1
%     Aeq = FE+TE;
%     beq = [zeros(3,1);[0 -9.8 0]'];
%   else
%     Aeq = FE;  
%     beq = [0 -9.8 0]';
%   end
   
  

end
