clear
close all
clf
hold on

tests = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%
% % example 0
% V = [0 1 0;
%     2 1 0];
% f = zeros(size(V,1),3);
% f(2,2) = -9.8;
% tests{1}{1} = V;
% tests{1}{2} = f;

% example 0
V = [0 1 0;
    2 1 0;
    3 1 0];
f = zeros(size(V,1),3);
f(3,3) = -9.8;
tests{1}{1} = V;
tests{1}{2} = f;

% % example 1
% V = [0 1 0;
%     2 1 0;
%     0 1 1;
%     2 1 1];
% f = zeros(size(V,1),3);
% f(4,2) = -9.8;
% tests{1}{1} = V;
% tests{1}{2} = f;

% example 1 broken
V = [0 1 0;
    1 1 0;
    2 1 0;
    0 1 1;
    1 1 1;
    2 1 1];
f = zeros(size(V,1),3);
f(6,2) = -9.8;
tests{2}{1} = V;
tests{2}{2} = f;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% example 2
V = [0 1 0;
    0 1 1;
    0 1 2;
    2 1 0;
    2 1 1;
    2 1 2];
f = zeros(size(V,1),3);
f(5,2) = -9.8;
tests{3}{1} = V;
tests{3}{2} = f;

% example 2 broken
V = [0 1 0;
    0 1 1;
    0 1 2;
    1 1 0;
    1 1 1;
    1 1 2;
    2 1 0;
    2 1 1;
    2 1 2];
f = zeros(size(V,1),3);
f(8,2) = -9.8;
tests{4}{1} = V;
tests{4}{2} = f;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% example 3
V = [0 0 0;
    0 1 0;
    0 0 1;
    2 0 0;
    2 0 1;
    2 1 0;
    0 1 1;
    2 1 1];
f = zeros(size(V,1),3);
f(8,2) = -9.8;
tests{5}{1} = V;
tests{5}{2} = f;

% example 3 broken
V = [0 0 0;
    0 1 0;
    0 0 1;
    1 0 0;
    1 0 1;
    1 1 0;
    0 1 1;
    1 1 1;
    2 0 0;
    2 0 1;
    2 1 0;
    2 1 1];
f = zeros(size(V,1),3);
f(size(V,1),2) = -9.8;
tests{6}{1} = V;
tests{6}{2} = f;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% % example 4
% nx = 4;
% ny = 2;
% nz = 3;
% 
% [V,E] = regular_grid_network([nx,ny,nz],3);
% V = (V-[1 (ny+1)/2 (nz+1)/2])/(nx-1);
% f = zeros(size(V,1),3);
% f(16,2) = -9.8;

tiledlayout(3,2)
for i = 1
V = tests{i}{1};
f = tests{i}{2};

E = [repmat((1:size(V,1))',size(V,1),1) repelem((1:size(V,1))',size(V,1),1)];
E = sort(E,2);
E = unique(E,'rows');
E(any(diff(E,[],2)==0,2),:)=[]; % get rid of rows of edges from a vertex to itself
[~,E] = prune_edges(V,E);

bf = find(V(:,1) == 0);
 
% plot boundary conditions
plot_edges(V,E);
scatter3(V(bf,1),V(bf,2),V(bf,3),'.b','SizeData',1000);

% plot force
quiver3(V(:,1),V(:,2),V(:,3),f(:,1),f(:,2),f(:,3),0.5,'r','LineWidth',3);


axis equal
camup([0 1 0])
cameratoolbar('setmode','orbit')
camproj('perspective');

% return

sC = 1e2;
sT = 1e2;
sB = 1e3;

dim=3;
n = size(V,1);
m = size(E,1);

EV = V(E(:,2),:)-V(E(:,1),:); % edge vectors
% EV = normalizerow(EV); % B matrix

[fr,fc] = find(f~=0);
fsum = sum(f(fr,:),1); % sum of forces
nforce = normalizerow(fsum);

MV = cross(EV,repmat(nforce,size(E,1),1)); % moment vectors
CV = cross(MV,normalizerow(EV)); % C matrix
% CV = cross(MV,EV); % C matrix
% CV = normalizerow(CV);

MP = (V(E(:,2),:)+V(E(:,1),:))/2; % midpoints for plotting

l = edge_lengths(V,E); % objective

BT = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[EV;-EV],dim*n,m);
CT = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[CV;-CV],dim*n,m);
% CT = sparse(E(:)+n*(0:dim-1),repmat(1:m,dim,2)',[MV;-MV],dim*n,m);

% remove fixed vertices
bf = bf(:);
BT(bf+n*(0:dim-1),:) = [];
CT(bf+n*(0:dim-1),:) = [];
ff = f;
ff(bf,:,:) = [];

I = speye(m,m);
Z = sparse(m,m);
nf = size(f,3);

Anb = [repmat(-sC*diag(1./l),nf,1), -I;...
       repmat(-sT*diag(1./l),nf,1),  I];
bnb = zeros(2*m*nf,1);

A = [repmat(-sC*diag(1./l),nf,1), -I,  Z;...
     repmat(-sT*diag(1./l),nf,1),  I,  Z;...
     repmat(-sB*diag(1./l),nf,1),  Z, -I;...
     repmat(-sB*diag(1./l),nf,1),  Z,  I];
% A = [repmat(-sC*I,nf,1), -I, Z;...
%     repmat(-sT*I,nf,1) , I, Z;...
%     repmat(-sB*I,nf,1), Z, -I;...
%     repmat(-sB*I,nf,1), Z, I];
b = zeros(4*m*nf,1);

B = [sparse(size(BT,1),m) BT sparse(size(BT,1),m)];
Bnb = [sparse(size(BT,1),m) BT];
C = [sparse(size(CT,1),2*m) CT];

% %%%%%%%%%%%%%% optimization %%%%%%%%%%%%%%%
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

% %%%%%%% without bending %%%%%%%%
% [an,~,flags,output] = linprog([l;zeros(m*nf,1)],Anb,bnb,...
%   Bnb,ff(:), ...
%   [zeros(m,1);-inf(m,1)], ...
%   [inf*ones(m,1);inf(m,1)]);
% arnb = an(1:m);
% axnb = reshape(an(m+(1:m*nf)),m,nf);

%%%%%%% with bending %%%%%%%%
[anb,~,flags,output] = linprog([l;zeros(2*m*nf,1)],A,b,...
  B+C,ff(:), ...
  [zeros(m,1);-inf(2*m,1)], ...
  [inf*ones(m,1);inf(2*m,1)]);
ar = anb(1:m);
ax = reshape(anb(m+(1:m*nf)),m,nf);
be = reshape(anb(2*m+(1:m*nf)),m,nf);


% tiledlayout(2,1)
% nexttile
% title(sprintf('Without Bending, Vol: %g, #bars: %d',l'*arnb,sum(log10(arnb)>-7)),'FontSize',20);
% hold on;
% plot_groundstructure(V,E,arnb,axnb);
% plot_edges(V,E,'Color','#848484','LineWidth',0.001);
% % plot boundary conditions
% scatter3(V(bf,1),V(bf,2),V(bf,3),'.b','SizeData',1000);
% % plot force
% quiver3(V(:,1),V(:,2),V(:,3),f(:,1),f(:,2),f(:,3),0.5,'r','LineWidth',3);
% axis equal
% view(115.69,53.858)
% colormap(flipud(cbrewer('RdBu',256)))
% caxis([-1 1])
% camup([0 1 0])
% cameratoolbar('SetCoordSys','y')
% cameratoolbar('setmode','orbit')
% camproj('perspective');
% 
nexttile
title(sprintf('With Bending, Vol: %g, #bars: %d',l'*ar,sum(max(ar,0)>1e-7)),'FontSize',20);
hold on;
plot_groundstructure(V,E,ar,be);
plot_edges(V,E,'Color','#848484','LineWidth',0.001);

% plot boundary conditions
scatter3(V(bf,1),V(bf,2),V(bf,3),'.b','SizeData',1000);

% plot force
quiver3(V(:,1),V(:,2),V(:,3),f(:,1)*0.1,f(:,2)*0.1,f(:,3)*0.1,'r','LineWidth',3,'AutoScale','off')
    
% plot other vertices
scatter3(V(:,1),V(:,2),V(:,3),'.k','SizeData',300);

% % plot direction of rotation
% quiver3(MP(:,1),MP(:,2),MP(:,3),CV(:,1),CV(:,2),CV(:,3),0.5,'c','LineWidth',3);

view(115.69,53.858)
colormap(flipud(cbrewer('RdBu',256)))
caxis([-1 1])
camup([0 1 0])
cameratoolbar('SetCoordSys','y') 
cameratoolbar('setmode','orbit')
camproj('perspective');
ylim([0 1])
zlim([0 2])
end

% plot_edges(V,E,'--k','LineWidth',3);
% scatter(V(bf,1),V(bf,2),'.b','SizeData',3000);
% quiver(V(:,1),V(:,2),f(:,1),f(:,2),0.5,'r','LineWidth',3);
% quiver(MV(:,1),MV(:,2),PV(:,1),PV(:,2),0.5,'g','LineWidth',3);
% quiver(MV(:,1),MV(:,2),EV(:,1),EV(:,2),0.5,'c','LineWidth',3);
% quiver(V(1+(1:numel(aa)),1),V(1+(1:numel(aa)),2),fa(:,1),fa(:,2),0,'b','LineWidth',3);
% hold off;
% axis equal
 
