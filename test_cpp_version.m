
E = readDMAT("GSME.dmat");
V = readDMAT("GSMV.dmat");

lc = readDMAT("lc.dmat");
uc = readDMAT("uc.dmat");

m = size(E,1);
n = size(V,1);

BT = mmread('BTmatrix.mtx');
A = mmread('Amatrix.mtx');
f = mmread('f.mtx');
% bf = readDMAT("bf.dmat");
% fe = readDMAT("fe.dmat");

nf = 1;

% vec = @(X) X(:);
% fe = [fe;fe];
% fe = fe+1;
% ig_ind = vec(fe(:)+m*(0:1))+m*(0:nf-1);
% A(ig_ind,:) = [];
% 

b = sparse(size(A,1),1);
% bf = [bf;bf];
% bf = bf+1;

vl = readDMAT("vandlterm.dmat");
vl = [vl;zeros(m,1)];

%%
% dim = 3;
% BT(bf+n*(0:dim-1),:) = [];
% f(bf,:,:) = [];

linear_c_m = mmread('lcm.mtx');
lcm = [A;BT];


max_a = inf;
[an,~,flags,output] = linprog(vl, A, uc(1:size(A,1)),...
    BT, f(:));%, ...
%     [zeros(m,1);-inf(m,1)]);, ...
%     [max_a*ones(m,1);inf(m,1)]);

a = an(1:m);
n = reshape(an(m+(1:m*nf)),m,nf);

%%
% find areas bigger than a certain threshold and place a cylinder there
NZ = find(max(a,0)>1e-4);
num_rods = size(NZ,1)
num_compression = sum(sign(n(NZ))==1)
num_tension = sum(sign(n(NZ))==-1)
[CV,CF,CJ,CI] = edge_cylinders(V,E(NZ,:),...
    'PolySize',10,'Thickness',sqrt(max(a(NZ),0)/pi));

color = sign(n(NZ(CJ)));

cylinders = tsurf(CF,CV,falpha(1,0),'CData',color,fsoft);

%%


acpp = readDMAT("areas.dmat");
ncpp = readDMAT("ax.dmat");


% % find areas bigger than a certain threshold and place a cylinder there
% NZ = find(max(acpp,0)>1e-4);
% num_rods = size(NZ,1)
% num_compression = sum(sign(ncpp(NZ))==1)
% num_tension = sum(sign(ncpp(NZ))==-1)
% [CV,CF,CJ,CI] = edge_cylinders(V,E(NZ,:),...
%     'PolySize',10,'Thickness',sqrt(max(acpp(NZ),0)/pi));
% 
% color = sign(ncpp(NZ(CJ)));
% 
% cylinders = tsurf(CF,CV,falpha(1,0),'CData',color,fsoft);

%%

% inequality test
ineq_s_mosek = (A*an <= b);
ineq_s_osqp = (A*[acpp;ncpp] <= b);

size(ineq_s_mosek)
sum(ineq_s_mosek)

size(ineq_s_osqp)
sum(ineq_s_osqp)

%%

% equality test
eq_s_osqp = ((BT*[acpp;ncpp] - f(:)) < 1e-4);
size(eq_s_osqp)
sum(eq_s_osqp)

eq_s_mosek = ((BT*an - f(:)) < 1e-4);
size(eq_s_mosek)
sum(eq_s_mosek)

%%


% sat_cpp_m_upper = (linear_c_m*an <= uc);
% sat_cpp_m_lower = (linear_c_m*an >= lc);

sat_cpp_cpp_upper = (linear_c_m*[acpp;ncpp] <= uc);
sat_cpp_cpp_lower = (linear_c_m*[acpp;ncpp] >= lc);
size(sat_cpp_cpp_upper)
sum(sat_cpp_cpp_upper)
size(sat_cpp_cpp_lower)
sum(sat_cpp_cpp_lower)

%%

sat_m_m_upper = (lcm*an <= [b;f(:)]);
sat_m_m_lower = (lcm*an >= [-inf(size(b,1),1);f(:)]);
size(sat_m_m_upper)
sum(sat_m_m_upper)
size(sat_m_m_lower)
sum(sat_m_m_lower)

% sat_m_cpp_upper = (lcm*[acpp;ncpp] <= uc);
% sat_m_cpp_lower = (lcm*[acpp;ncpp] >= lc);

%%
isequal(lcm,linear_c_m)
