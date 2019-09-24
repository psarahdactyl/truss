[MV,MF] = bwmesh('hammer.png','Tol',2);
MV = MV*0.05;
th = 0.4;
R = [cos(th) -sin(th);sin(th) cos(th)];
MV = MV*R;

% attachment points
%[~,la] = min(MV(:,1)+1e26*(MV(:,2)>0.25*max(MV(:,2))));
%[~,ua] = min(MV(:,1)+1e26*(MV(:,2)<0.5*max(MV(:,2))));
%aa = [la;ua];
na = 2;
aa = unique(snap_points([min(MV(:,1))*ones(na,1) linspace(min(MV(:,2)),max(MV(:,2)),na)'],MV));
%aa = (1:size(MV,1))';
na = numel(aa);

bbd = normrow(max(MV)-min(MV));
W = linspace(min(MV(:,2)),max(MV(:,2)),20)';
W(:,2) = min(MV(:,1))-0.5*bbd;
W = fliplr(W);
[I,J] = find(ones(size(W,1),numel(aa)));
[C,vol] = centroid(MV,outline(MF));
V = [C;MV(aa,:);W];
E = [[ones(1,numel(aa));1+(1:numel(aa))]';I+1+numel(aa) J+1];

f = zeros(size(V));
f(1,2) = -9.8;


bf = 1+numel(aa)+1:size(V,1);
sC = 5e1;
sT = 5e1;
[a,n,l,BT] = groundstructure(V,E,f,bf,sC,sT,'IgnoredEdges',1:numel(aa));


naa = numel(aa);
av = 1 + (1:naa)';
ae = (1:naa)';
fa = reshape(BT(av + size(V,1)*[0 1],naa+1:end)*n(naa+1:end),[],2);
cross2 = @(A,B) A(:,1).*B(:,2)-A(:,2).*B(:,1);
torque = sum(cross2(fa,V(E(ae,2),:)-V(E(ae,1),:)))

fbf = reshape(BT(bf(:) + [0 1]*size(V,1),:)*n,[],2);

clf;
hold on;
tsurf(MF,MV,'FaceColor',[0.6 0.5 0.3],falpha(1,0));
%tsurf(E,V,'LineWidth',2);
plot_groundstructure(V,E,a,n);
plot_edges(V,E(1:numel(aa),:),'--k','LineWidth',3);
scatter(V(1,1),V(1,2),'.r','SizeData',bbd);
scatter(V(1+(1:numel(aa)),1),V(1+(1:numel(aa)),2),'.k','SizeData',bbd);
scatter(V(bf,1),V(bf,2),'.b','SizeData',bbd);
quiver(V(:,1),V(:,2),f(:,1),f(:,2),0,'r','LineWidth',3);
quiver(V(bf,1),V(bf,2),fbf(:,1),fbf(:,2),0,'.b','LineWidth',3);
quiver(V(1+(1:numel(aa)),1),V(1+(1:numel(aa)),2),fa(:,1),fa(:,2),0,'b','LineWidth',3);
hold off;
axis equal
title(sprintf('Vol: %g, #bars: %d',l'*a,sum(log10(a)>-4)),'FontSize',20);
