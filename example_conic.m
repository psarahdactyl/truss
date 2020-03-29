
function [V,T] = createAMS(V,T,viewpoint,lambdaMin,lambdaMax,fixedVertices,weights,mu,zOrder)
% createAMS implements the method of Appearance-Mimicking Surfaces
%
% Author: Christian Schüller
% Email: schuellc[@]inf.ethz.ch
% Project page: http://igl.ethz.ch/projects/ams/
% 
% Dependencies: libigl - http://igl.ethz.ch/projects/libigl/
%
% createAMS(V,T,viewpoint,lambdaMin,lambdaMax,fixedVertices,weights,mu,zOrder)
%
% Inputs:
%   V             #Vertices by 3 list of the vertex positions of the model
%   T             #Faces by 3 list of triangle indices into V
%   viewpoint     1 by 3 vector of the coordinates of the viewpoint
%   lambdaMin     #Constraints by 2 list of [vertexIdx,lambdaMinValue]
%   lambdaMax     #Constraints by 2 list of [vertexIdx,lambdaMaxValue]
%   fixedVertices vector of size #FixedVertices [vertexIdx]
%   weights       vector of size #Vertices [weightFactor]
%   mu            vector of size #Vertices with correspondinf mu index
%   zOrder        #Constraints x 8 [order,vIdx,vAIdx,vBIdx,vCIdx,u,v,w]
%                 Contains the depth order constraints
%                 for vertex <-> triangle point pairs
%                 order: 1 vertex in front of triangle; -1 behind triangle
%                 vIdx:  index of vertex
%                 vAIdx: index of first vertex (A) of triangle
%                 vBIdx: index of second vertex (B) of triangle
%                 vCIdx: index of third vertex (C) of triangle
%                 u: barycentric coordinate of point of triangle A 
%                 v: barycentric coordinate of point of triangle B
%                 w: barycentric coordinate of point of triangle C
% Outputs:
%   V             #Vertices by 3 list of the vertex positions of the deformed model
%   T             #Faces by 3 list of triangle indices into V
% Comments:
%   If you work with high resolution meshes or irregular triangulations
%   using MOSEK might bee faster and numerically more stable. Set MOSEK to 1.
    
% to use MOSEK set to 1
MOSEK = 0;

n = size(V,1);
alpha = 1e-6;

nMu = max(mu); % # independent thickness constraint regions
nZOrder = size(zOrder,1); % # z-order constraints

% tranform to origin (viewpoint)
V = V - repmat(viewpoint,n,1);

% gather inequality constraint matrix entry tripples [roxIdx,colIdx,value]
nlambdaMin = size(lambdaMin,1);
nlambdaMax = size(lambdaMax,1);
nCon = nlambdaMin + nlambdaMax + nZOrder;

rId = 0;
IC_min = [];
IC_max = [];
IC_z = [];

% lambda_min
if nlambdaMin>0
    idx = [(rId+1):(rId+nlambdaMin)]';
    rId = rId+nlambdaMin;
    rowC0 = idx;
    colC0 = lambdaMin(:,1);
    valC0 = ones(n,1)*-1;
    rowC1 = idx;
    colC1 = n+mu(lambdaMin(:,1));
    valC1 = lambdaMin(:,2);
    IC_min = [[rowC0;rowC1],[colC0;colC1],[valC0;valC1]];
end

% lambda_min
if nlambdaMax>0
    idx = [(rId+1):(rId+nlambdaMax)]';
    rId = rId+nlambdaMax;
    rowC0 = idx;
    colC0 = lambdaMax(:,1);
    valC0 = ones(n,1);
    rowC1 = idx;
    colC1 = n+mu(lambdaMax(:,1));
    valC1 = -lambdaMax(:,2);
    IC_max = [[rowC0;rowC1],[colC0;colC1],[valC0;valC1]];
end

% z-order constraints
if nZOrder>0
    idx = [(rId+1):(rId+nZOrder)]';
    rId = rId+nZOrder;
    rowC0 = idx;
    colC0 = zOrder(:,2);
    valC0 = zOrder(:,1);
    rowC1 = idx;
    colC1 = zOrder(:,3);
    valC1 = -zOrder(:,6);
    rowC2 = idx;
    colC2 = zOrder(:,4);
    valC2 = -zOrder(:,7);
    rowC3 = idx;
    colC3 = zOrder(:,5);
    valC3 = -zOrder(:,8);
    IC_z = [[rowC0;rowC1;rowC2;rowC3],[colC0;colC1;colC2;colC3],[valC0;valC1;valC2;valC3]];
end

IC_t = [IC_min;IC_max;IC_z];

% gather equality constraints entries [vIdx,value] (used to fix a vertex)
EC = [fixedVertices,sqrt(sum(V([fixedVertices],:).^2,2))];

% compute laplacian
L=-cotmatrix(V,T);
M=inv(massmatrix(V,T,'voronoi'));
M=kron(M,spdiags(ones(3),0,3,3));

% compute sqrt of mass matrix
Msqrt=spdiags(sqrt(diag(M)),0,3*n,3*n);

% direction matrix
R = normr(V);

% weighting matrix
W=spdiags(kron(weights,[1;1;1]),0,3*n,3*n);

% multiplier matrix (1 lambda per vertex (x,y,z))
P=kron(speye(n),[1;1;1]);

% selector matrix for lambda/mu
S=[speye(n),sparse(n,nMu)];

% Inequality constraint matrix
C=sparse(IC_t(:,1),IC_t(:,2),IC_t(:,3),nCon,n+nMu);
b=sparse(nCon,1);

% Compose energy matrices
l_0=sqrt(sum(V'.^2))';
v=reshape(V',3*n,1);
d=reshape(R',3*n,1);
D=spdiags(d,0,3*n,3*n);

% compute scale invariant laplacian
L_d=spdiags(diag(L),0,n,n);
L_s=spdiags(1./l_0,0,n,n)*(L-L_d)*spdiags(l_0,0,n,n) + L_d;
L=kron(L,speye(3));
L_s=kron(L_s,speye(3));
H=D*P;
G=spdiags(L_s*d,0,3*n,3*n)*P;
F=W*(L*H-G);
f=sparse(n,1)';
F=Msqrt*F;
Q=F'*F;

% compute small alpha weight factor
Qdiag = diag(Q);
alpha0=alpha*sum(Qdiag)/size(Qdiag,1);

% add mu variables
f=f*S;
f(n+1:end)=ones(nMu,1)*0;

[rows,cols] = size(F);
I=sparse(rows+nMu,cols+nMu);
I(1:rows,1:cols)=F;
I(rows+1:end,cols+1:end)=speye(nMu)*sqrt(alpha0);
F=I;

Q=F'*F;

if(MOSEK) % Cast it to conic program if MOSEK is available
    nLambda = n+nMu;
    [nt,cols] = size(F);

    % slack variable vector t
    t=sparse(size(F,1),1);

    % objective function vector
    obj=[f,t',1,0]';

    % inequality constraint matrix
    con=[[F speye(nt,nt) sparse(nt,2)];[C sparse(size(C,1),nt+2)]];

    % right hand side vectors of inequality constraints
    blc=[sparse(nt,1);ones(size(C,1),1)*-inf;];
    buc=[sparse(nt,1);sparse(size(C,1),1)];

    % bounds on variables
    blx=[ones(nLambda+nt,1)*-inf;0;1];
    bux=[ones(nLambda+nt+1,1)*inf;1];
    blx(EC(:,1),:)=EC(:,2);
    bux(EC(:,1),:)=EC(:,2);

    % cone constraint matrix
    cone=[nLambda+nt+2,nLambda+nt+1,[(nLambda+1):(nLambda+nt)]];

    % initialize problem
    [r,res] = mosekopt('symbcon');
    prob = [];
    prob.c = obj;
    prob.a = con;
    prob.buc = buc;
    prob.blc = blc;
    prob.bux = bux;
    prob.blx = blx;
    prob.cones.type = [res.symbcon.MSK_CT_RQUAD];
    prob.cones.sub = cone;
    prob.cones.subptr = [1];

    % debugging
    mosekopt('anapro',prob);

    % solve problem
    [r,res] = mosekopt('minimize info',prob);

    % get result
    lambdas = res.sol.itr.xx(1:nLambda);
else
    [lambdas,F,Lambda] = min_quad_with_fixed_active_set(Q,2*f',EC(:,1),EC(:,2),[],[],C,b,[],[],[],'MaxIter',1000);
end

% create deformed surface
V = repmat(viewpoint,n,1) + spdiags(lambdas(n+mu),0,n,n)*spdiags(lambdas(1:n),0,n,n)*R;

end