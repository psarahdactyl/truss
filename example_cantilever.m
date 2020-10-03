rng(0)
V = {};
F = {};
coms = [nan nan nan];

clf
hold on

% XX=zeros(3,3);
% XX = [0 0 0;
%     2 1 0;
%     0 1 0];
% XE=[2 3];
% XC=[1 2 3];

XXX = [0 1 0;
     2 1 0;
     4 1 0];
XXE=[1 2
  2 3];
XXC=[1 2 3];
% XXC=[1 2];

CB = cbrewer('Set1',5);
cbred = CB(1,:);
cbblue = CB(2,:);
cbgreen = CB(3,:);
cborange = CB(5,:);

% units are meters
% force = [0 -9.8 0];
force = [0 0 -9.8];
sC = ones(size(XXE,1),1)*1e2;
sT = ones(size(XXE,1),1)*1e2;
sB = ones(size(XXE,1),1)*1e2;

%%
axis equal
% plot3(XX(XE(:),1),XX(XE(:),2),XX(XE(:),3))
plot_edges(XXX,XXE)
scatter3(XXX(:,1),XXX(:,2),XXX(:,3),'.b','SizeData',1000)
quiver3(XXX(3,1),XXX(3,2),XXX(3,3),force(1),force(2),force(3),0.1)

bending=1;
[A,b,Aeq,beq] = create_constraint_matrices(XXX,XXE,XXC,coms,force,sC,sT,sB,...
  'Bending',bending);
objective = edge_lengths(XXX,XXE);

%%
% TROUBLE SHOOTING
EV = XXX(XXE(:,2),:)-XXX(XXE(:,1),:); % tangent edge vectors for C (compression)
BV1 = zeros(size(EV));
BV2 = zeros(size(EV));
for i=1:size(XXE,1)
  bs = null(EV(i,:));
  BV1(i,:) = norm(EV(i,:))*bs(:,1)'; % normal vectors for B (bending)
  BV2(i,:) = norm(EV(i,:))*bs(:,2)'; % binormal vectors for B (bending)
end

n=size(XXX,1)
m=size(XXE,1)
dim=size(XXX,2)
NEV=normalizerow(EV);
C = sparse(XXE(:)+n*(0:dim-1),repmat(1:m,dim,2)',[NEV;-NEV],dim*n,m);
OB1 = sparse(XXE(:)+n*(0:dim-1),repmat(1:m,dim,2)',[BV1;-BV1],dim*n,m);
OB2 = sparse(XXE(:)+n*(0:dim-1),repmat(1:m,dim,2)',[BV2;-BV2],dim*n,m);
B = [OB1 OB2];
% B = OB1;
  
bf = find(XXX(:,1) == 0);
bf = bf(:);
B(bf+n*(0:dim-1),:) = [];
C(bf+n*(0:dim-1),:) = [];
ff = repmat(force',max(XXC)-1,1);
FE = [sparse(size(C,1),m) C B];
Aeq=FE;
beq=ff;
beq(3)=0;
full(Aeq)

%%
% a=full(Aeq);
% quiver3(XX(3,1),XX(3,2),XX(3,3),a(1,2),a(1,3),a(1,4))
% quiver3(XX(3,1),XX(3,2),XX(3,3),a(2,2),a(2,3),a(2,4))
% quiver3(XX(3,1),XX(3,2),XX(3,3),a(3,2),a(3,3),a(3,4))
% [x,ar,ax,be,fval] = optimize_lp(objective,A,b,Aeq(4:end,:),beq(4:end,:),'linprog','Bending',bending);
[x,ar,ax,be,fval] = optimize_lp(objective,A,b,Aeq,beq,'linprog','Bending',bending);

NZ = ar>1e-6;
[CV,CF] = edge_cylinders(XXX,XXE(NZ,:),'Thickness',0.15*sqrt(ar(NZ)),'PolySize',30);

csh = tsurf(CF,CV, ...
  falpha(1,0),fsoft,'FaceVertexCData',repmat(cborange,size(CV,1),1));
