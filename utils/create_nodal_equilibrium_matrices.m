function [S,G,C,B] = create_nodal_equilibrium_matrices(V,E,ACV,coms)
  % CREATE_NODAL_EQUILIBRIUM_MATRICES 
  % 
  %
  % Inputs:
  %   V   #V by dim list of vertex positions
  %   E   #E by 2 list of edges indices into V
  %   ACV #AV by 1 list of object indices into rows of AV
  %   coms #max(ACV) by 3 list of each object's center of mass

  % Outputs:
  %   B  dim by #E nodal equilibrium matrix for bending
  %   T  dim by #E nodal equilibrium matrix for torque
  %   F  dim by #E nodal equilibrium matrix for tension/compression
  
  
%   V = [1 1 1; 4 2 1;4 3 1;4 4 1;8 5 1]
%   E = [1 2;
%        1 3;
%        1 4;
%        2 5;
%        3 5;
%        4 5]
%   ACV = [1;2;2;2;3]
  
  
  dim = size(V,2);
  n = size(V,1);
  m = size(E,1);
  
%   "selection" matrix 3*k by 3*V, with 1s to "select" vertices incident on
%                        rigid object $i$, so that $\S_i$ has 3 rows so that $\S_i
%                       \f$ gathers just the total forces on object $i$

%   [0 \S\C \S\B] [\a] = 0
%   [0 \G\C \G\B] [\c]
%                 [\b]
% 


% % S here is a selection matrix
%   S = sparse(1:q,find(~b),ones(size(1:q)),q,size(V,1));
%   S = blkdiag(S,S);

  S = sparse(dim*(max(ACV)-1),dim*size(V,1));
  for i=2:max(ACV)
    vs=find(ACV==i);
    rows=(1:3)+3*(i-2);
    % rows for object i
    cols=[vs;vs+size(V,1);vs+2*size(V,1)]; 
    % columns in the selection matrix for vertices belonging to object 1
%     S(rows,cols)=1;
    S(1+3*(i-2),vs)=1;
    S(2+3*(i-2),vs+size(V,1))=1;
    S(3+3*(i-2),vs+2*size(V,1))=1;
  end
    
  CP = @(x) ([0    -x(3)  x(2); 
              x(3)  0    -x(1); 
             -x(2)  x(1)  0 ]); % cross product matrix
  G = sparse(dim*(max(ACV)-1),dim*size(V,1));
  for i=2:max(ACV)
    vs=find(ACV==i);
    ds=V(vs,:)-repmat(coms(i,:),size(vs,1),1);
    rows=(1:3)+3*(i-2);
    for j=1:size(ds,1)
      cpm=CP(ds(j,:));
      G(rows,vs(j))=cpm(:,1);
      G(rows,vs(j)+size(V,1))=cpm(:,2);
      G(rows,vs(j)+2*size(V,1))=cpm(:,3);
    end
  end
  
  EV = V(E(:,2),:)-V(E(:,1),:); % tangent edge vectors for C (compression)
  BV1 = zeros(size(EV));
  BV2 = zeros(size(EV));
  for i=1:size(E,1)
    b = null(EV(i,:));
    BV1(i,:) = norm(EV(i,:))*b(:,1)'; % normal vectors for B (bending)
    BV2(i,:) = norm(EV(i,:))*b(:,2)'; % binormal vectors for B (bending)
  end

  NEV=normalizerow(EV);
  C = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[NEV;-NEV],dim*n,m);
  OB1 = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[BV1;-BV1],dim*n,m);
  OB2 = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[BV2;-BV2],dim*n,m);
  B = [OB1 OB2];

%   
%   T = zeros((max(ACV)-1)*dim,m); % torque
%   F = zeros((max(ACV)-1)*dim,m); % torque
%   num_obj = max(ACV);
%   
%   % for every object that is not the wall
%   for i=2:num_obj
%     av=find(ACV==num_obj); % all vertices in said object
%     com=coms(i,:); % center of mass
%     Ti = zeros(dim,m); % torque
%     Fi = zeros(dim,m); % forces
%     % for every edge
%     for e=1:m
%       ei = E(e,:);
%       ev = EV(e,:);
%       nev = ev/norm(ev); % normalized edge vector
%       cross2 = @(A,B) A(:,1).*B(:,2)-A(:,2).*B(:,1);
% 
%       if(ismember(ei(1),av) && ismember(ei(2),av))
%         continue
%         
%       elseif(ismember(ei(1),av))
%         if dim==2
%           Ti(:,e) = cross2(nev,V(ei(1),:)-com); % units = m
%         elseif dim==3
%           Ti(:,e) = cross(nev,V(ei(1),:)-com); % units = m
%         end
%         Fi(:,e) = ev;
%         
%       elseif(ismember(ei(2),av))
%         if dim==2
%           Ti(:,e) = cross2(nev,V(ei(1),:)-com); % units = m
%         elseif dim==3
%           Ti(:,e) = cross(nev,V(ei(1),:)-com); % units = m
%         end
%         Fi(:,e) = ev;
%       end
%       
%     end
%     T((1:3)+3*(i-2),:) = Ti;
%     F((1:3)+3*(i-2),:) = Fi;
%   end
% 
%   size(T)
%   size(F)
  size(G)
  size(S)
  size(C)
  size(B)
  
end
