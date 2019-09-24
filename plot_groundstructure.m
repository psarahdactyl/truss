function plot_groundstructure(V,E,a,n)

EV = normalizerow(V(E(:,2),:)-V(E(:,1),:));
m = size(E,1);


dim = size(V,2);
switch dim
  case 2
  VV = [ ...
    V(E(:,1),:) + 0.5*a.*(EV*[0 1;-1 0]); ...
    V(E(:,1),:) - 0.5*a.*(EV*[0 1;-1 0]); ...
    V(E(:,2),:) - 0.5*a.*(EV*[0 1;-1 0]); ...
    V(E(:,2),:) + 0.5*a.*(EV*[0 1;-1 0])];
  QQ = (1:m)'+m*(0:3);
tsurf(QQ,VV,'Tets',0,falpha(1,0.05),'CData',repmat(n,4,1));
  case 3
  NZ = find(max(a,0)>1e-7);
  [CV,CF,CJ,CI] = edge_cylinders(V,E(NZ,:), ...
    'PolySize',10,'Thickness',sqrt(max(a(NZ),0)/pi));
  tsurf(CF,CV,falpha(1,0.05),'CData',n(NZ(CJ),1));
end
  

end

