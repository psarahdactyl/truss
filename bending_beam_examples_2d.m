clear
close all
clf


%         |
%         v
% 2 - 4 - 6 
% |   |   | 
% 1 - 3 - 5
%
groundstructure1.V = [0 0;
     0 1;
     1 0;
     1 1;
     2 0;
     2 1];
 
groundstructure1.E = [1 2;
     1 3;
     2 3;
     1 4;
     3 6;
     4 5;
     2 4;
     3 4;
     4 6;
     3 5;
     5 6];

groundstructure1.f = zeros(size(groundstructure1.V,1),2);
groundstructure1.f(6,2) = -9.8;
groundstructure1.bf = [1 2];

%%%%%%%%%%%%%
% 3 --      |
%      -    v
% 2 ------- 4 (length = 2)
%     -
% 1 --
%
groundstructure2.V = [0 0;
    0 1;
    0 2;
    2 1];
 
groundstructure2.E = [1 4;
    2 4;
    3 4];

groundstructure2.f = zeros(size(groundstructure2.V,1),2);
groundstructure2.f(4,2) = -9.8;
groundstructure2.bf = [1 2 3];

%%%%%%%%%%%%%
%           |
%           v
% 1 -- 2 -- 3 (lengths = 1,1)
%
groundstructure3.V = [0 0;
    1 0;
    2 0];
 
groundstructure3.E = [1 2;
    2 3];

groundstructure3.f = zeros(size(groundstructure3.V,1),2);
groundstructure3.f(3,2) = -9.8;
groundstructure3.bf = [1];

%%%%%%%%%%%%%
%          |
%          v
% 1 ------ 2 (length = 2)
%
groundstructure4.V = [0 0;
    2 0];
 
groundstructure4.E = [1 2];

groundstructure4.f = zeros(size(groundstructure4.V,1),2);
groundstructure4.f(2,2) = -9.8;
groundstructure4.bf = [1];

%%%%%%%%%%%%%
%        |        %
%        v        %
% 1 ---- 2 ---- 3 %(lengths = 2,2)
%
groundstructure5.V = [0 0;
    2 0;
    4 0];
 
groundstructure5.E = [1 2;
    2 3];

groundstructure5.f = zeros(size(groundstructure5.V,1),2);
groundstructure5.f(2,2) = -9.8;
groundstructure5.bf = [1 3];

%%%%%%%%%%%%%
%           |           %
%           v           %
% 1 -- 2 -- 3 -- 4 -- 5 %(lengths = 1,1,1,1)
% % 
groundstructure6.V = [0 0;
    1 0;
    2 0;
    3 0;
    4 0];
 
groundstructure6.E = [1 2;
    2 3;
    3 4;
    4 5];

groundstructure6.f = zeros(size(groundstructure6.V,1),2);
groundstructure6.f(3,2) = -9.8;
groundstructure6.bf = [1 5];

structures = {groundstructure1,groundstructure2,...
                groundstructure3,groundstructure4,...
                groundstructure5,groundstructure6,};

sC = 1e2;
sT = 1e2;
sB = 1e3;

dim=2;

tiledlayout(3,2) 

for i=1:6 
    V = structures{i}.V;
    E = structures{i}.E;
    f = structures{i}.f;
    bf = structures{i}.bf;
    
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
        repmat(-sB*diag(1./l),nf,1), Z, I];
    b = zeros(4*m*nf,1);


    B = [sparse(size(BT,1),m) BT sparse(size(BT,1),m)];
    Bnb = [sparse(size(BT,1),m) BT];
    C = [sparse(size(CT,1),2*m) CT];

    %%%%%%%%%%%%%% optimization %%%%%%%%%%%%%%%

%     %%%%%%% without bending %%%%%%%%
%     [an,~,flags,output] = linprog([l;zeros(m*nf,1)],Anb,bnb,...
%       Bnb,ff(:), ...
%       [zeros(m,1);-inf(m,1)], ...
%       [inf*ones(m,1);inf(m,1)]);
%     arnb = an(1:m);
%     axnb = reshape(an(m+(1:m*nf)),m,nf);

    %%%%%%% with bending %%%%%%%%
    [anb,~,flags,output] = linprog([l;zeros(2*m*nf,1)],A,b,...
      B+C,ff(:), ...
      [zeros(m,1);-inf(2*m,1)], ...
      [inf*ones(m,1);inf(2*m,1)]);
    ar = anb(1:m);
    ax = reshape(anb(m+(1:m*nf)),m,nf);
    be = reshape(anb(2*m+(1:m*nf)),m,nf);

    nexttile
    title(sprintf('With Bending, Vol: %g, #bars: %d',l'*ar,sum(log10(ar)>-7)),'FontSize',20);
    hold on;
    plot_edges(V,E,'k','LineWidth',0.1)
    plot_groundstructure(V,E,ar,be);
    quiver(V(:,1),V(:,2),f(:,1)*0.1,f(:,2)*0.1,'r','LineWidth',3,'AutoScale','off')
    scatter(V(:,1),V(:,2),'.k','SizeData',300);
    scatter(V(bf,1),V(bf,2),'.b','SizeData',1000);
    ylim([-2 2])
%     quiver(MV(:,1),MV(:,2),PV(:,1),PV(:,2),0.5,'g','LineWidth',3);
    colormap(flipud(cbrewer('RdBu',256)))
    caxis([-1 1])

end
hold off;
 
