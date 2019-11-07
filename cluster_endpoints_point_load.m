% hold on
function [nV, nE, nNZ, nn, na] = cluster_endpoints_point_load(V,E,NZ,n,H,GV,sC,sT,ideal_k,force_vertices)
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

%     [~,cT] = kmeans(RT,kt);
%     [~,cC] = kmeans(RC,kc);
    [~,cB] = kmeans([RT;RC],ideal_k);

    % split cT and cC into two sides
%     C1 = [cT(:,1:3); cC(:,1:3)];
%     C2 = [cT(:,4:6); cC(:,4:6)];
    C1 = cB(:,1:3);
    C2 = cB(:,4:6);
    
    scatter3(cB(:,1),cB(:,2),cB(:,3),'.g','SizeData',1000); % C1
    
    nV = [C1;C2];
    nE = [(1:size(cB,1))' ((1:size(cB,1))+size(cB,1))'];
    [SV,SVI,SVJ] = remove_duplicate_vertices(nV);
    SE = SVJ(nE);

    hold on;
    plot_edges(SV,SE,'Color',[0.4940 0.1840 0.5560],'LineWidth',3);

    f = zeros(size(SV));
    [tf, index] = ismember(SV,force_vertices,'rows');
    f(tf,2) = -9.8;
    [bv, ~] = ismember(SV,C2,'rows');
    bf = find(bv);

    [na,nn,l,h,BT] = groundstructure(SV,SE,H,GV,f,bf,sC,sT);

    nNZ = find(max(na,0)>1e-4);
    num_rods_clustered = size(nNZ,1)
    num_compression_clustered = sum(sign(nn(nNZ))==1)
    num_tension_clustered = sum(sign(nn(nNZ))==-1)
    
    plot_edges(nV,nE(nNZ,:),'Color',[0.4940 0.1840 0.5560],'LineWidth',3);


%     [CV,CF,CJ,CI] = edge_cylinders(nV,nE(nNZ,:),...
%             'PolySize',10,'Thickness',sqrt(max(a(nNZ),0)/pi));

%     tsurf(CF,CV,falpha(1,0),'CData',sign(n(nNZ(CJ))),fsoft);

end
