E = readDMAT("GSME.dmat");
V = readDMAT("GSMV.dmat");
m = size(E,1);
n = size(V,1);

lower_c = readDMAT("lc.dmat");
upper_c = readDMAT("uc.dmat");

A_matlab = mmread('BTmatrix.mtx');
A_cpp = mmread('BTmatrix.mtx');

D_matlab = mmread('Amatrix.mtx');
D_cpp = mmread('Amatrix.mtx');

b_matlab = mmread('f.mtx');
b_matlab = b_matlab(:);
b_cpp = upper_c(size(D_cpp,1)+1:end);

e_matlab = sparse(size(D_matlab,1),1);
e_cpp = upper_c(1:size(D_cpp,1));

L_cpp = readDMAT("vandlterm.dmat");
L_cpp = [L_cpp;zeros(m,1)];

L_matlab = readDMAT("vandlterm.dmat");
L_matlab = [L_matlab;zeros(m,1)];

lcm_cpp = mmread('lcm.mtx');
lcm_matlab = [D_matlab;A_matlab];

assert(isequal(A_cpp,A_matlab))
assert(isequal(D_cpp,D_matlab))
assert(isequal(e_cpp,e_matlab))
assert(isequal(b_cpp,b_matlab))
assert(isequal(L_cpp,L_matlab))
assert(isequal(lcm_cpp,lcm_matlab))

%%
A = [A_matlab;-A_matlab;D_matlab];
b = [b_matlab;-b_matlab;e_matlab];

[u,~,flags,output] = linprog(L_matlab, A, b);

%%
B = [lcm_cpp;-lcm_cpp];
y = [e_cpp;b_cpp];
x = lower_c;
z = [y;-x];

[u,~,flags,output] = linprog(L_cpp, B, z);

%%
u1 = u(1:m);
u2 = u(m+1:end);

NZ = find(max(u1,0)>1e-4);
num_rods = size(NZ,1)
num_compression = sum(sign(u2(NZ))==1)
num_tension = sum(sign(u2(NZ))==-1)
[CV,CF,CJ,CI] = edge_cylinders(V,E(NZ,:),...
    'PolySize',10,'Thickness',sqrt(max(u1(NZ),0)/pi));

color = sign(u2(NZ(CJ)));

cylinders = tsurf(CF,CV,falpha(1,0),'CData',color,fsoft);

