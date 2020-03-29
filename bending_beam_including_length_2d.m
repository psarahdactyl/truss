clear
close all
clf

%         |
%         v
% 2 - 4 - 6 
% |   |   | 
% 1 - 3 - 5
%
% V = [0 0;
%      0 1;
%      1 0;
%      1 1;
%      2 0;
%      2 1];
%  
% E = [1 2;
%      1 3;
%      2 3;
%      1 4;
%      3 6;
%      4 5;
%      2 4;
%      3 4;
%      4 6;
%      3 5;
%      5 6];
% 
% f = zeros(size(V,1),2);
% f(6,2) = -9.8;
% bf = [1 2];

%%%%%%%%%%%%%
% 3 --      |
%      -    v
% 2 ------- 4 (length = 2)
%     -
% 1 --
%
% V = [0 0;
%     0 1;
%     0 2;
%     2 1];
%  
% E = [1 4;
%     2 4;
%     3 4];
% 
% f = zeros(size(V,1),2);
% f(4,2) = -9.8;
% bf = [1 2 3];

%%%%%%%%%%%%%
%           |
%           v
% 1 -- 2 -- 3 (lengths = 1,1)
%
V = [0 0;
    1 0;
    2 0];
 
E = [1 2;
    2 3];

f = zeros(size(V,1),2);
f(3,2) = -9.8;
bf = [1];

%%%%%%%%%%%%%
%          |
%          v
% 1 ------ 2 (length = 2)
%
% V = [0 0;
%     2 0];
%  
% E = [1 2];
% 
% f = zeros(size(V,1),2);
% f(2,2) = -9.8;
% bf = [1];

% plot_edges(V,E);
 
sC = 1e2;
sT = 1e2;
sB = 1e3;

dim=2;
n = size(V,1);
m = size(E,1);

EV = V(E(:,2),:)-V(E(:,1),:);
% EV = normalizerow(EV); % B matrix

PV = V(E(:,2),:)-V(E(:,1),:);
% PV = normalizerow(PV);
PV = [PV(:,2) -PV(:,1)]; % C matrix

[fr,fc] = find(f~=0);
fsum = sum(f(fr,:),1);

MV = (V(E(:,2),:)+V(E(:,1),:))/2; % midpoints for plotting
l = edge_lengths(V,E); % objective
  
BT = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[EV;-EV],dim*n,m);
CT = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[PV;-PV],dim*n,m);
% return


% remove fixed vertices
bf = bf(:);
BT(bf+n*(0:dim-1),:) = [];
CT(bf+n*(0:dim-1),:) = [];
ff = f;
ff(bf,:,:) = [];

I = speye(m,m);
Z = sparse(m,m);
nf = size(f,3);

Anb = [repmat(-sC*I,nf,1), -I; ...
    repmat(-sT*I,nf,1),  I];
bnb = zeros(2*m*nf,1);

A = [repmat(-sC*diag(1./l),nf,1), -I, Z; ...
    repmat(-sT*diag(1./l),nf,1) , I, Z;...
    repmat(-sB*diag(1./l),nf,1), Z, -I;...
    repmat(-sB*diag(1./l),nf,1), Z, -I];
b = zeros(4*m*nf,1);


B = [sparse(size(BT,1),m) BT sparse(size(BT,1),m)];
Bnb = [sparse(size(BT,1),m) BT];
C = [sparse(size(CT,1),2*m) CT];

%%%%%%%%%%%%%% optimization %%%%%%%%%%%%%%%
% o = [l; zeros(2*m,1)];
% 
% % Define variables
% x = sdpvar(3*m,1);
% 
% % Define constraints
% Constraints = [A * x <= b, B * x + C * x == fbf(:)];
% 
% % Define an objective
% Objective = o' * x;
% 
% % Set some options for YALMIP and solver
% options = sdpsettings('verbose',1);
% 
% sol = optimize(Constraints,Objective,options);
% 
% % Analyze error flags
% if sol.problem == 0
%  % Extract and display value
%  solution = value(x);
%  ar = solution(1:m);
%  ax = reshape(solution(m+(1:m*nf)),m,nf);
%  be = reshape(solution(2*m+(1:m*nf)),m,nf);
% else
%  display('Hmm, something went wrong!');
%  sol.info
%  yalmiperror(sol.problem)
% end

%%%%%%% without bending %%%%%%%%
[an,~,flags,output] = linprog([l;zeros(m*nf,1)],Anb,bnb,...
  Bnb,ff(:), ...
  [zeros(m,1);-inf(m,1)], ...
  [inf*ones(m,1);inf(m,1)]);
arnb = an(1:m);
axnb = reshape(an(m+(1:m*nf)),m,nf);

%%%%%%% with bending %%%%%%%%
[anb,~,flags,output] = linprog([l;zeros(2*m*nf,1)],A,b,...
  B+C,ff(:), ...
  [zeros(m,1);-inf(2*m,1)], ...
  [inf*ones(m,1);inf(2*m,1)]);
ar = anb(1:m);
ax = reshape(anb(m+(1:m*nf)),m,nf);
be = reshape(anb(2*m+(1:m*nf)),m,nf);


tiledlayout(3,1) 
% nexttile
% title(sprintf('Without Bending, Vol: %g, #bars: %d',l'*arnb,sum(log10(arnb)>-4)),'FontSize',20);
% hold on;
% plot_groundstructure(V,E,arnb,axnb);
% scatter(V(bf,1),V(bf,2),'.b','SizeData',3000);
% quiver(V(:,1),V(:,2),f(:,1),f(:,2),0.5,'r','LineWidth',3);
% % quiver(MV(:,1),MV(:,2),PV(:,1),PV(:,2),0.5,'g','LineWidth',3);
% axis equal

nexttile
title(sprintf('With Bending, Vol: %g, #bars: %d',l'*ar,sum(log10(ar)>-4)),'FontSize',20);
hold on;
plot_groundstructure(V,E,ar,ax);
quiver(V(:,1),V(:,2),f(:,1),f(:,2),0.5,'r','LineWidth',3)
scatter(V(:,1),V(:,2),'.k','SizeData',1000);
scatter(V(bf,1),V(bf,2),'.b','SizeData',3000);
% quiver(MV(:,1),MV(:,2),PV(:,1),PV(:,2),0.5,'g','LineWidth',3);
axis equal

% plot_edges(V,E,'--k','LineWidth',3);
% scatter(V(bf,1),V(bf,2),'.b','SizeData',3000);
% quiver(V(:,1),V(:,2),f(:,1),f(:,2),0.5,'r','LineWidth',3);
% quiver(MV(:,1),MV(:,2),PV(:,1),PV(:,2),0.5,'g','LineWidth',3);
% quiver(MV(:,1),MV(:,2),EV(:,1),EV(:,2),0.5,'c','LineWidth',3);
% quiver(V(1+(1:numel(aa)),1),V(1+(1:numel(aa)),2),fa(:,1),fa(:,2),0,'b','LineWidth',3);
hold off;
axis equal
 
