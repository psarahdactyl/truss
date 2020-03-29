clear
close all
clf

%
%           .
%         .    .
%        .  1   .
%         .    .
%           .
% 
%
% 

% radius = [7 6];
% theta = 0:pi/10:2*pi;
% xs = 1 * cos(theta) + radius(1);
% ys = 1 * sin(theta) + radius(2);
% 
% nw = 9; % number of wall points
% 
% V = [radius;  % com
%     xs' ys'; % circle points (21)
%     zeros(nw,1) 1+(1:nw)']; % wall points
% 
% E = [ones(size(xs,2),1) 1+(1:size(xs,2))';
%     repmat(1+(1:size(xs,2))',nw,1) repelem((2+size(xs,2):size(V,1))',size(xs,2),1)];
% 
% f = zeros(size(V,1),2);
% f(1,2) = -9.8;
% bf = 2+size(xs,2):size(V,1);
% 
% ignored_edges = 1:size(xs,2);

% plot_edges(V,E(1+size(xs,2):end,:));
% return

%%%%%%%%%%%%%
%           |           %
%           v           %
% 1 -- 2 -- 3 -- 4 -- 5 %(lengths = 1,1,1,1)
% % 
% V = [0 0;
%     1 0;
%     2 0;
%     3 0;
%     4 0];
%  
% E = [1 2;
%     2 3;
%     3 4;
%     4 5];
% 
% f = zeros(size(V,1),2);
% f(3,2) = -9.8;
% bf = [1 5];

%%%%%%%%%%%%%
%        |        %
%        v        %
% 1 ---- 2 ---- 3 %(lengths = 2,2)
%
V = [0 0;
    2 0;
    4 0];
 
E = [1 2;
    2 3];

f = zeros(size(V,1),2);
f(2,2) = -9.8;
bf = [1 3];

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
%      --   v
% 2 ------- 4 (length = 2)
%      --
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
% |      -- 3
% v   --
% 4 ------- 2 (length = 2)
%     --
%        -- 1
% %
% V = [2 0;
%     2 1;
%     2 2;
%     0 1];
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
% V = [0 0;
%     1 0;
%     2 0];
%  
% E = [1 2;
%     2 3];
% 
% f = zeros(size(V,1),2);
% f(3,2) = -9.8;
% bf = 1;

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
% bf = 1;

%%%%%%%%%%%%%
% |        %
% v        %
% 1 ------ 2 (length = 2)
%          %
% V = [0 0;
%     2 0];
%  
% E = [1 2];
% 
% f = zeros(size(V,1),2);
% f(1,2) = -9.8;
% bf = 2;

% plot_edges(V,E);
 
sC = 1e2;
sT = 1e2;
sB = 1e2;

dim=2;
n = size(V,1);
m = size(E,1);

EV = V(E(:,2),:)-V(E(:,1),:); % edge vectors
% EV = normalizerow(EV); % B matrix

PV = V(E(:,2),:)-V(E(:,1),:); % perpendicular vectors
% PV = normalizerow(PV);
PV = [PV(:,2) -PV(:,1)]; % C matrix

[fr,fc] = find(f~=0);
fsum = sum(f(fr,:),1);

MP = (V(E(:,2),:)+V(E(:,1),:))/2; % midpoints for plotting
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

Anb = [repmat(-sC*diag(1./l),nf,1), -I; ...
    repmat(-sT*diag(1./l),nf,1),  I];
bnb = zeros(2*m*nf,1);

A = [repmat(-sC*diag(1./l),nf,1), -I, Z; ...
    repmat(-sT*diag(1./l),nf,1) , I, Z;...
    repmat(-sB*diag(1./l),nf,1), Z, -I;...
    repmat(-sB*diag(1./l),nf,1), Z, I];
b = zeros(4*m*nf,1);

% ig_ind = vec(ignored_edges(:)+m*(0:1))+m*(0:nf-1);
% Anb(ig_ind,:) = [];
% bnb(ig_ind,:) = [];
% 
% big_ig_ind = [ignored_edges;
%     ignored_edges+m;
%     ignored_edges+2*m;
%     ignored_edges+3*m];
% A(big_ig_ind(:),:) = [];
% b(big_ig_ind(:),:) = [];

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
% Constraints = [A * x <= b, B * x + C * x == ff(:)];
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
% 
% return


% l(ignored_edges) = 0;

%%%%%%% without bending %%%%%%%%
[an,~,flags,output] = linprog([l;zeros(m*nf,1)],Anb,bnb,...
  Bnb,ff(:), ...
  [zeros(m,1);-inf(m,1)], ...
  [inf*ones(m,1);inf(m,1)]);
arnb = an(1:m);
axnb = reshape(an(m+(1:m*nf)),m,nf);

%%%%%%% with bending %%%%%%%%
[anwb,~,flags,output] = linprog([l;zeros(2*m*nf,1)],A,b,...
  B+C,ff(:), ...
  [zeros(m,1);-inf(2*m,1)], ...
  [inf*ones(m,1);inf(2*m,1)]);
ar = anwb(1:m);
ax = reshape(anwb(m+(1:m*nf)),m,nf);
be = reshape(anwb(2*m+(1:m*nf)),m,nf);


tiledlayout(2,1) 
nexttile
title(sprintf('Without Bending, Vol: %g, #bars: %d',l'*arnb,sum(log10(arnb)>-7)),'FontSize',20);
hold on;
plot_groundstructure(V,E,arnb,axnb);
scatter(V(:,1),V(:,2),'.k','SizeData',500);
scatter(V(bf,1),V(bf,2),'.b','SizeData',1000);
quiver(V(:,1),V(:,2),f(:,1),f(:,2),0.5,'r','LineWidth',3);
% quiver(MV(:,1),MV(:,2),PV(:,1),PV(:,2),0.5,'g','LineWidth',3);
axis equal

nexttile
title(sprintf('With Bending, Vol: %g, #bars: %d',l'*ar,sum(log10(ar)>-7)),'FontSize',20);
hold on;
plot_groundstructure(V,E,ar,be);
scatter(V(:,1),V(:,2),'.k','SizeData',500);
scatter(V(bf,1),V(bf,2),'.b','SizeData',1000);
% quiver(MV(:,1),MV(:,2),PV(:,1),PV(:,2),0.5,'g','LineWidth',3);
quiver(V(:,1),V(:,2),f(:,1),f(:,2),0.5,'r','LineWidth',3);
axis equal

colormap(flipud(cbrewer('RdBu',256)))

% plot_edges(V,E,'--k','LineWidth',3);
% scatter(V(bf,1),V(bf,2),'.b','SizeData',3000);
% quiver(V(:,1),V(:,2),f(:,1),f(:,2),0.5,'r','LineWidth',3);
% quiver(MV(:,1),MV(:,2),PV(:,1),PV(:,2),0.5,'g','LineWidth',3);
% quiver(MV(:,1),MV(:,2),EV(:,1),EV(:,2),0.5,'c','LineWidth',3);
% quiver(V(1+(1:numel(aa)),1),V(1+(1:numel(aa)),2),fa(:,1),fa(:,2),0,'b','LineWidth',3);
hold off;

 
