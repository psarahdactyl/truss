function [a,n,l,OBT] = groundstructure(V,E,f,bf,sC,sT,varargin)
% [a,n] = groundstructure(V,E,f,bf,sC,sT,...)
%
% Inputs:
%   V  #V by dim list of node locations
%   E  #E by 2 list of edge indices into V
%   f  #V by dim by #f list of loading conditions
%   bf  #bf list of Dirichlet boundaries; indices into V (nodes where force
%     balance is not required.)
%   sC  maximum compression stress
%   sT  maximum tension stress
% Outputs:
%   a  #E list of cross-sectional widths/areas
%   n  #E by #loads list of axial forces on each edge for each load


  ignored_edges = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'IgnoredEdges'}, ...
    {'ignored_edges'});
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

  dim = size(V,2);
  nf = size(f,3);

  function [X,mqwf] = admm_argmin_aN(l,BT,F,Sp,Sn,rho,Ua,UN,mqwf)
    %energy = @(a,N,Sp,Sn,Ua,UN,rho) ...
    %  l'*a + 
    %  ;
    %% Debug using global solve
    %% min l'a + (ρ/2)‖A*X + B*[Sp(:);Sn(:)] + [Ua(:);UN(:)]‖²
    %% subject to BT*N = F;
    %Aeq = [sparse(size(BT,1)*nf,m) repdiag(BT,nf)];
    %Beq = F(:);
    %Y = B*[Sp(:);Sn(:)] + [Ua(:);UN(:)];
    %% min X'[l;0] + (ρ/2)X'*A'A*X + ρ X' A'Y + constant
    %% subject to BT*N = F;
    %Xmo = quadprog(rho*A'*A,[l;zeros(m*nf,1)]+rho*A'*Y,[],[],Aeq,Beq);
    %amo = Xmo(1:m);
    %Nmo = reshape(Xmo(m+(1:m*nf)),m,nf);

    % Optimize for a and N separately
    %
    % min l'a + (ρ/2) ∑∑ (ai - (s⁺ij/σt + s⁻ij/σc) + uaij )²
    % let yij = (s⁺ij/σt + s⁻ij/σc) - uaij
    % min li ai + (ρ/2) ∑ (ai - yij )²
    % min li ai + (ρ/2) ∑ (ai² - 2ai yij + yij²)
    % let m = #j
    % min li ai + (ρ/2)m ai² - ρ (∑ yij) ai + constant
    % ρm ai = ρ (∑ yij) - li
    % ai = (∑ yij)/m - li/(ρm)
    Y = Sp/sT + Sn/sC - Ua;
    a = sum(Y,2)/nf - l/(rho*nf);
    %
    % min (ρ/2) ‖ N + [-I I][Sp;Sn] + UN‖²  subject to B'N = F
    % min (ρ/2) ‖ N + (-Sp+Sn+UN)‖²  subject to B'N = F
    % let W = Sp-Sn-UN
    % min (ρ/2) ‖N - W‖²  subject to B'N = F
    % Cool/suspicious: ρ doesn't matter.
    % min ‖N - W‖²  subject to B'N = F
    W = Sp-Sn-UN;
    I = speye(size(W,1));
    %% Augh... apparently the constraints are not li
    %mqwf.force_Aeq_li = false;
    %[N,mqwf] = min_quad_with_fixed(I,-2*W,[],[],BT,-F,mqwf);
    % Deal with constraints the lazy way
    wc = 1e10;
    [N,mqwf] = min_quad_with_fixed(I+wc*BT'*BT,-2*W-2*wc*BT'*F,[],[],[],[],mqwf);
    X = [a(:);N(:)];
  end

  function [Z,data] = admm_argmin_SpSn(a,N,Ua,UN)
    %% Debug using global solve
    %% min (ρ/2)‖A*X+B*Z+U‖²
    %% Z≥0
    %% min (ρ/2)Z'*B'B*Z + ρ Z'*B'*(A*X+U) + constant
    %% Z≥0
    %Zmo = quadprog( ...
    %  B'*B,B'*(A*[a(:);N(:)]+[Ua(:);UN(:)]),[],[],[],[],zeros(m*nf*2,1),[]);
    %Spmo = reshape(Zmo(1:m*nf),m,nf);
    %Snmo = reshape(Zmo(m*nf+(1:m*nf)),m,nf);

    % min (ρ/2) (ai - s⁺ij/σt - s⁻ij/σc + uaij)² + (ρ/2)(nij - s⁺ij + s⁻ij + unij)²
    % subject to s⁺ij ≥ 0
    %            s⁻ij ≥ 0
    %
    % Again, ρ has disappeared...
    %
    % min (ai - s⁺ij/σt - s⁻ij/σc + uaij)² + (nij - s⁺ij + s⁻ij + unij)²
    % subject to s⁺ij ≥ 0
    %
    %  min  [s⁺ij s⁻ij] / 1/σt²+1    1/(σtσc)-1 \ /s⁺ij\  + [s⁺ij s⁻ij] / -ai/σt-uaij/σt-nij-unij \ + constant
    %                   \ 1/(σtσc)-1 1/σc²+1    / \s⁻ij/                \ -ai/σc-uaij/σc+nij+unij /
    % subject to s⁺ij ≥ 0
    %            s⁻ij ≥ 0
    %
    % let wij = - / -ai/σt-uaij/σt-nij-unij \
    %             \ -ai/σc-uaij/σc+nij+unij /
    % let sij = [s⁺ij s⁻ij]'
    %
    % Again, ρ has disappeared...
    %
    %  min  (1/2) sij' Q sij  - sij' wij 
    % subject to sij ≥ 0
    % 
    % tiny QP... system matrix is even constant w.r.t. ij FWIW
    %
    % Exhaustive Algorithm:
    %   sij_11 = Q\wij
    %   sij_01 = [0 Q(2,2)\wij(2)]
    %   sij_10 = [Q(1,1)\wij(1) 0]
    %   sij_00 = [0 0]
    %   Eij_11 = ½ sij_11'*Q*sij_11 - sij_11'*wij
    %   Eij_01 = ½ sij_01'*Q*sij_01 - sij_01'*wij
    %   Eij_10 = ½ sij_10'*Q*sij_10 - sij_10'*wij
    %   Eij_00 = ½ sij_00'*Q*sij_00 - sij_00'*wij
    %   Eij_11(any(sij_11<0)) = inf;
    %   Eij_01(any(sij_01<0)) = inf;
    %   Eij_10(any(sij_10<0)) = inf;
    %   Eij_00(any(sij_00<0)) = inf;
    %   sij = sij_argmin_uv Eij_uv
    %
    Q = [1/(sT^2)+1 1/(sT*sC)-1;1/(sT*sC)-1 1/(sC^2)+1];
    W = -cat(3,-a/sT - Ua/sT - N - UN,-a/sC-Ua/sC+N+UN);
    % feelin' dirty
    Qinv = inv(Q);

    O = inf(m,nf);
    Sp = nan(m,nf);
    Sn = nan(m,nf);
    for ci = 1:4
      switch ci
      case 1
        SSp = Qinv(1,1)*W(:,:,1) + Qinv(1,2)*W(:,:,2);
        SSn = Qinv(2,1)*W(:,:,1) + Qinv(2,2)*W(:,:,2);
      case 2
        SSp = zeros(m,nf);
        SSn = W(:,:,2)/Q(2,2);
      case 3
        SSp = W(:,:,1)/Q(1,1);
        SSn = zeros(m,nf);
      case 4
        SSp = zeros(m,nf);
        SSn = zeros(m,nf);
      end
      OO = ...
        0.5*( ...
          Q(1,1)*SSp.^2 + ...
          2*Q(1,2)*SSp.*SSn + ...
          Q(2,2)*SSn.^2) - (W(:,:,1).*SSp+W(:,:,2).*SSn);
      OO(SSp<0 | SSn<0) = inf;
      Sp(OO<O) = SSp(OO<O);
      Sn(OO<O) = SSn(OO<O);
      O(OO<O) = OO(OO<O);
    end
    Z = [Sp(:);Sn(:)];

    %% Alternative
    %SS = ...
    %  cat(4, ...
    %    cat(3, ...
    %      Qinv(1,1)*W(:,:,1) + Qinv(1,2)*W(:,:,2), ...
    %      Qinv(2,1)*W(:,:,1) + Qinv(2,2)*W(:,:,2)), ...
    %    cat(3,zeros(size(W,1),size(W,2)),W(:,:,2)/Q(2,2)), ...
    %    cat(3,W(:,:,1)/Q(1,1),zeros(size(W,1),size(W,2))), ...
    %    zeros(size(W,1),size(W,2),2));
    %O = ...
    %  0.5*( ...
    %    Q(1,1)*SS(:,:,1,:).^2 + ...
    %    2*Q(1,2)*SS(:,:,1,:).*SS(:,:,2,:) + ...
    %    Q(2,2)*SS(:,:,2,:).^2) - sum(W.*SS,3);
    %O = permute(O,[1 2 4 3]);
    %O(permute(any(SS<0,3),[1 2 4 3])) = inf;
    %S =  sum(permute(repmat(O == min(O,[],3),[1 1 1 size(SS,3)]),[1 2 4 3]).*SS,4);
    %Sp = S(:,:,1);
    %Sn = S(:,:,2);
    %Z = [Sp(:);Sn(:)];

    data = [];
  end

  EV = normalizerow(V(E(:,2),:)-V(E(:,1),:));
  % bar lengths
  l = edge_lengths(V,E);
  % use ignored edges for free
  %l(ignored_edges) = 0;


  n = size(V,1);
  m = size(E,1);

  BT = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[EV;-EV],dim*n,m);
  OBT = BT;


  bf = bf(:);
  % remove fixed values
  BT(bf+n*(0:dim-1),:) = [];
  f(bf,:,:) = [];


  max_a = inf;
  % -sC a ≤ n ≤ sT a
  %
  % -sC a ≤ n
  %     n ≤ sT a
  %
  % -sC a + -n ≤ 0
  % -sT a +  n ≤ 0

  if(exist('mosekopt','file'))
    % Sometimes this is _very_ slow
    %params = default_quadprog_param;
    params = [];
  else
    params = [];
  end

  % manually reformulate in "standard" lp form (no linear inequalities)
  %method = 'standard';
  method = 'linprog';
  %method = 'admm';
  switch method
  case 'standard'
    assert(isempty(ignored_edges));
    A = [speye(m,m)./sT speye(m,m)./sC];
    A = [repmat(A,nf-1,1) repdiag(-A,nf-1)];
    b = zeros(size(A,1),1);
    s = linprog( ...
       repmat([l*sC;l*sT],nf,1), ...
       [],[], ...
       [repdiag([BT, -BT],nf);A],[f(:);b], ...
       zeros(2*m*nf,1),inf(2*m*nf,1),[],params);
    s = reshape(s,[m,2,nf]);
    sp = permute(s(:,1,:),[1 3 2]);
    sn = permute(s(:,2,:),[1 3 2]);
    a = mean(sp/sT+sn/sC,2);
    n = sp-sn;
  case 'linprog'
    BT = repdiag(BT,nf);
    II = speye(m*nf,m*nf);
    I = speye(m,m);
    A = [repmat(-sC*I,nf,1) , -II; ...
        repmat(-sT*I,nf,1) ,  II];
    b = zeros(2*m*nf,1);
    vec = @(X) X(:);
    ig_ind = vec(ignored_edges(:)+m*(0:1))+m*(0:nf-1);
    A(ig_ind,:) = [];
    b(ig_ind,:) = [];
    % Something weird is going on with the sign...
    [an,~,flags,output] = linprog([l;zeros(m*nf,1)],A,b,[sparse(size(BT,1),m) BT],f(:), ...
      [zeros(m,1);-inf(m,1)], ...
      [max_a*ones(m,1);inf(m,1)], ...
      [], ...
      params);
    a = an(1:m);
    n = reshape(an(m+(1:m*nf)),m,nf);
  case 'admm'
    % 
    %  min l'a 
    %  a,N
    %  subject to BN = F
    %             -σc ai ≤ nij ≤ σt ai ∀ i,j
    %
    %  introduce s⁺,s⁻
    %
    %    min     l'a 
    %  a,N,s⁺,s⁻
    %  subject to   BN = F
    %               ai = s⁺ij/σt + sij⁻/σc ∀ i,j
    %              nij = s⁺ij - s⁻ij
    %             s⁺ij ≥ 0
    %             s⁻ij ≥ 0
    %
    %     min     f(a,N) + g(s⁺,s⁻)
    %  a,N,S⁺,S⁻
    %  subject to  a + /-I/σt -I/σc\ /S⁺\  = /0\
    %              N   \-I     I   / \S⁻/    \0/
    % where
    % f(a,N) = l'a + / 0  if BN=F
    %                \ ∞  otherwise
    % g(s⁺,s⁻) = / 0  if s⁺ij ≥ 0 and s⁻ij ≥ 0
    %            \ ∞  otherwise
    % 
    % Let X = [a;N(:)]
    % Let Z = [Sp(:);Sn(:)]
    % Let U = [Ua(:);UN(:)]
    %
    A = [repmat(speye(m,m+m*nf),nf,1); sparse(m*nf,m) speye(m*nf)];
    B = [-speye(m*nf)/sT -speye(m*nf)/sC;-speye(m*nf) speye(m*nf)];
    c = zeros(size(A,1),1);
    %     min     f(a,N) + g(Sp,Sn)
    %  a,N,Sp,Sn
    %     subject to: A*[a;N(:)] + B*[Sp(:);Sn(:)] - c = 0

    Sp_fun = @(Z) reshape(Z(      1:m*nf),m,nf);
    Sn_fun = @(Z) reshape(Z(m*nf+(1:m*nf)),m,nf);
    Ua_fun = @(U) reshape(U(      1:m*nf),m,nf);
    UN_fun = @(U) reshape(U(m*nf+(1:m*nf)),m,nf);
    F = reshape(f,size(f,1)*size(f,2),size(f,3));
    argmin_X = @(Z,U,rho,data) admm_argmin_aN(l,BT,F,Sp_fun(Z),Sn_fun(Z),rho,Ua_fun(U),UN_fun(U),data);
    a_fun = @(X) X(1:m);
    N_fun = @(X) reshape(X(m+(1:m*nf)),m,nf);
    argmin_Z = @(X,U,rho,data) admm_argmin_SpSn(a_fun(X),N_fun(X),Ua_fun(U),UN_fun(U));
    state = [];

    %
    state.X = zeros(size(A,2),size(c,2));
    state.Z = zeros(size(B,2),size(c,2));
    state.U = zeros(size(c,1),size(c,2));
    % a = s⁺/σt + s⁻/σc
    % n = s⁺ - s⁻
    % n + s⁻ = s⁺
    % a = (n + s⁻)/σt + s⁻/σc
    % a = n/σt + s⁻/σt + s⁻/σc
    % a - n/σt = s⁻/σt + s⁻/σc
    % a - n/σt = s⁻(1/σt + 1/σc)
    % s⁻ = (a - n/σt)/(1/σt + 1/σc)
    % s⁺ = n + s⁻
    %load('an.mat','a','n');
    %Sn = (a - n/sT)./(1/sT + 1/sC);
    %Sp = n + Sn;
    %state.X = [a(:);n(:)];
    %state.Z = [Sp(:);Sn(:)];

    state.rho_prev = nan;
    state.rho = 1e-3;
    state.argmin_X_data = [];
    state.argmin_Z_data = [];
    for iter = 1:10000
      tic;
      [X,Z,state] = admm(argmin_X,argmin_Z,A,B,c,state,'MaxIter',1);
      toc
      fprintf('im alive\n');
      a = a_fun(X);
      n = N_fun(X);
      %if mod(iter,1) == 1
        aeff = a;
        subplot(2,1,1);
      %plot_groundstructure(V,E,aeff,n(:,1));
      tsurf(E(aeff>1e-7,:),V,'CData',aeff(aeff>1e-7,:));
      axis equal;
      title(sprintf('Vol: %g, #bars: %d/%d, iter: %d',l'*aeff,sum(aeff>1e-7),size(E,1),iter-0),'FontSize',20);
        subplot(2,1,2);
        aeff = mean(Sp/sT + Sn/sC,2);
      tsurf(E(aeff>1e-7,:),V,'CData',aeff(aeff>1e-7,:));
      %plot_groundstructure(V,E,aeff,n(:,1));
      axis equal;
      title(sprintf('Vol: %g, #bars: %d, iter: %d',l'*aeff,sum(aeff>1e-7),iter-0),'FontSize',20);
      max(abs(BT*n - F))
      drawnow;
      %figgif('admm-2d-gs.gif');
      %pause
      %end
    end


  end
  
end
