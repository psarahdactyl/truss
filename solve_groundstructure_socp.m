function [a,n,BT] = solve_groundstructure_socp(V,E,f,bf,sC,sT,h,l,varargin)
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

  EV = normalizerow(V(E(:,2),:)-V(E(:,1),:));

%   test_h = h;
%   test_h(ignored_edges) = 0;
%   test_l = l;
%   test_l(ignored_edges) = 0;
%   
%   max_length = max(test_l)
%   min_length = min(test_l)
%   
%   max_vis = max(test_h)
%   min_vis = min(test_h)

%   obj = l+300*h;
  obj = l.*h;
  
%   max_objective = max(obj)
%   min_objective = min(obj)

  % use ignored edges for free
  obj(ignored_edges) = 0;

  nnn = size(V,1);
  m = size(E,1);

  BT = sparse(E(:)+nnn*(0:dim-1),repmat(1:m,dim,2)',[EV;-EV],dim*nnn,m);
  OBT = BT;

  bf = bf(:);
  
  % remove fixed values
  BT(bf+nnn*(0:dim-1),:) = [];
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

    opt = mskoptimset('');
    %     opt = mskoptimset(opt,'Diagnostics','on');
    %     opt = mskoptimset(opt,'MSK_IPAR_INFEAS_REPORT_AUTO', 'MSK_ON');

    [an,~,flags,output] = linprog([obj;zeros(m*nf,1)], A, b,...
      [sparse(size(BT,1),m) BT], f(:), ...
      [zeros(m,1);-inf(m,1)], ...
      [max_a*ones(m,1);inf(m,1)],...
      opt);
  
    a = an(1:m);
    n = reshape(an(m+(1:m*nf)),m,nf);
    max_area = max(a)
%%
% check error
% 
%     err = obj.*a;
%     i = find(err>1e-10);
%     [i err(i)];
% 
%     % check constraints are being satisfied
%     inequality_constraints = sum(A*an-b>1e-7)
%     Bb = [sparse(size(BT,1),m) BT];
%     size(Bb)
%     equality_constraints = norm([sparse(size(BT,1),m) BT]*an-f(:))
% 
%     other_err = BT*n;
%     [other_err f(:) other_err - f(:)];
%     %     i = find(other_err>1e-10)
%     %     other_err(i)
  
end
