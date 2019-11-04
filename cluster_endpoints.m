% hold on
function [nV, nE, nNZ, nn, na] = cluster_endpoints(V,E,NZ,n,H,GV,sC,sT,cm,ideal_k)
    % separate rods in tension and rods in compression
    % and also separate them by object
    RT1 = V(E(NZ(sign(n(NZ))==-1),1),:);
    RT2 = V(E(NZ(sign(n(NZ))==-1),2),:);
    RC1 = V(E(NZ(sign(n(NZ))== 1),1),:);
    RC2 = V(E(NZ(sign(n(NZ))== 1),2),:);

    RT = [RT1 RT2];
    RC = [RC1 RC2];

    kt = ideal_k;
    kc = ideal_k;
    if ideal_k > size(RC,1)
        kc = size(RC,1);
    end
    if ideal_k > size(RT,1)
        kt = size(RT,1);
    end

    [~,cT] = kmeans(RT,kt);
    [~,cC] = kmeans(RC,kc);
%     cT = unique(cT,'rows');
%     cC = unique(cC,'rows');

    % sanity check here
%     cT = RT;
%     cC = RC;

    % split cT and cC into two sides
    C1 = [cT(:,1:3); cC(:,1:3)];
    C2 = [cT(:,4:6); cC(:,4:6)];
%     C1 = V(E(NZ(sign(n(NZ))==-1|sign(n(NZ))==1),1),:);
%     C2 = V(E(NZ(sign(n(NZ))==-1|sign(n(NZ))==1),2),:);
%     C1 = unique(C1,'rows');
%     C2 = unique(C2,'rows');
    
    scatter3(cT(:,1),cT(:,2),cT(:,3),'.g','SizeData',1000); % first part C1
    scatter3(cT(:,4),cT(:,5),cT(:,6),'.y','SizeData',1000); % first part C2
    scatter3(cC(:,1),cC(:,2),cC(:,3),'.r','SizeData',1000); % second part C1
    scatter3(cC(:,4),cC(:,5),cC(:,6),'.b','SizeData',1000); % second part C2

    nV = [cm;C1;C2];
    nE = [ones(size(C1,1),1) 1+(1:size(C1,1))';...
        1+[(1:size(cT,1))' ((1:size(cT,1))+size(cT,1)+size(cC,1))'];...
        1+size(cT,1)+[(1:size(cC,1))' ((1:size(cC,1))+size(cC,1)+size(cT,1))']];

%     [I,J] = find(ones(size(C1,1),size(C2,1)));
%     nE = [ones(size(C1,1),1) 1+(1:size(C1,1))';1+[I size(C1,1)+J]];

    hold on;
%     plot_edges(nV,nE,'Color',[0.4940 0.1840 0.5560],'LineWidth',3);

    f = zeros(size(nV));
    bl = snap_points(cm,nV);
    f(bl,2) = -9.8;
    bf = 1+size(C1,1)+(1:size(C2,1));

    [na,nn,l,h,BT] = groundstructure(nV,nE,H,GV,f,bf,sC,sT,'IgnoredEdges',1:size(C1,1));

    naa = size(C1,1);
    av = 1 + (1:naa)';
    ae = (1:naa)';
    dim = size(nV,2);
    fa = reshape(BT(av + size(nV,1)*[0:dim-1],naa+1:end)*nn(naa+1:end),[],dim);
    ntorque = normrow(sum(cross(fa,nV(nE(ae,2),:)-nV(nE(ae,1),:),2)))

    nNZ = find(max(na,0)>1e-4);
    num_rods_clustered = size(nNZ,1)
    num_compression_clustered = sum(sign(nn(nNZ))==1)
    num_tension_clustered = sum(sign(nn(nNZ))==-1)
    
    plot_edges(nV,nE(nNZ,:),'Color',[0.4940 0.1840 0.5560],'LineWidth',3);


%     [CV,CF,CJ,CI] = edge_cylinders(nV,nE(nNZ,:),...
%             'PolySize',10,'Thickness',sqrt(max(a(nNZ),0)/pi));

%     tsurf(CF,CV,falpha(1,0),'CData',sign(n(nNZ(CJ))),fsoft);

end
