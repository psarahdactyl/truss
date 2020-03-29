function [h segments] = edge_visibilities(V, E, GV, side, w, H, l)
% H : hiddenness on voxel grid
% V, E : vertices and edges
% GV : voxel grid bottom left corners
% l : edge lengths

h = zeros(size(l));

X = reshape(GV(:,1),side);
Y = reshape(GV(:,2),side);
Z = reshape(GV(:,3),side);

origins = V(E(:,1),:);
dests = V(E(:,2),:);
segments = ceil(l/(w(1)/2));

rep_origins = repelem(origins,max(segments),1);
rep_dests = repelem(dests,max(segments),1);

rep_slopes = rep_dests - rep_origins;

ts = zeros(size(segments,1),max(segments));

% yuck
for i=1:size(segments,1)
    p = linspace(0,1,segments(i));
    ts(i,1:numel(p)) = p;
end

% size(ts)

tss = reshape(ts',size(ts,1)*size(ts,2),1);

sts = rep_origins + rep_slopes.*tss;

scs = interp3(X,Y,Z,reshape(H,side),sts(:,1),sts(:,2),sts(:,3));
% size(H)
% sc = scatteredInterpolant(GV(:,1),GV(:,2),GV(:,3),-H');
% scs = sc(sts(:,1),sts(:,2),sts(:,3));

h = sum(reshape(1-scs,max(segments),size(E,1)));
% max(h)
% min(h)
% 
% h = h./(2*segments)';
% max(h)
% min(h)

end