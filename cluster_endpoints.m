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

%     kt = ideal_k;
%     kc = ideal_k;
%     if ideal_k > size(RC,1)
%         kc = size(RC,1);
%     end
%     if ideal_k > size(RT,1)
%         kt = size(RT,1);
%     end

%     [~,cT] = kmeans(RT,kt);
%     [~,cC] = kmeans(RC,kc);

    rods = [RT;RC];
%     [rods, ridx] = sortrows(rods);
    [~,cB] = kmeans(rods,ideal_k);
    how_many_clustered_rods = size(cB)
%     cB = cB(ridx,:);

    % sanity check
%     cB = rods;

    % split cT and cC into two sides
%     C1 = [cT(:,1:3); cC(:,1:3)];
%     C2 = [cT(:,4:6); cC(:,4:6)];
    C1 = [cB(:,1:3)];
    C2 = [cB(:,4:6)];
    
%     scatter3(cB(:,1),cB(:,2),cB(:,3),'.g','SizeData',1000); % C1
    
    nV = [cm;C1;C2];
    nE = [ones(size(C1,1),1) 1+(1:size(C1,1))';...
        1+[(1:size(cB,1))' ((1:size(cB,1))+size(cB,1))']];
    [SV,SVI,SVJ] = remove_duplicate_vertices(nV);
    SE = SVJ(nE);
%     SE = unique(SE,'rows');

    hold on;
%     plot_edges(SV,SE,'Color',[0.4940 0.1840 0.5560],'LineWidth',3);

    f = zeros(size(SV));
    bl = snap_points(cm,SV);
    f(bl,2) = -9.8;
    [bv, ~] = ismember(SV,C2,'rows');
    bf = find(bv);
    ie = find(SE == 1);

    [na,nn,l,h,BT] = groundstructure(SV,SE,H,GV,f,bf,sC,sT,'IgnoredEdges',ie');

    naa = size(C1,1);
    av = 1 + (1:naa)';
    ae = (1:naa)';
    dim = size(SV,2);
%     fa = reshape(BT(av + size(SV,1)*[0:dim-1],naa+1:end)*nn(naa+1:end),[],dim);
%     ntorque = normrow(sum(cross(fa,SV(SE(ae,2),:)-SV(SE(ae,1),:),2)))

    nNZ = find(max(na,0)>1e-4);
    num_rods_clustered = size(nNZ,1)
    num_compression_clustered = sum(sign(nn(nNZ))==1)
    num_tension_clustered = sum(sign(nn(nNZ))==-1)
    
%     plot_edges(SV,SE(nNZ,:),'Color',[0.4940 0.1840 0.5560],'LineWidth',3);

end
