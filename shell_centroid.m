function [cen,vol] = shell_centroid(V,F)
  M = diag(massmatrix(V,F));
  vol = sum(M);
  cen = M'*V/vol;
end
