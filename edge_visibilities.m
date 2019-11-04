function h = edge_visibilities(V, E, GV, H, l)
% H : hiddenness on voxel grid
% V, E : vertices and edges
% GV : voxel grid bottom left corners
% l : edge lengths

    h = zeros(size(l));
%     side = [19 16 20];
%     side = [20 20 15];
%     side = [20 20 20];
    side = [18 26 50];
    half_voxel_width = abs(GV(1,1) - GV(2,1))/2;
    voxel_width = abs(GV(1,1) - GV(2,1));
    for i=1:size(E,1)
        num_segments = ceil(l(i)/half_voxel_width)+1;
        origin = V(E(i,1),:);
        rep_origin = repmat(origin, num_segments, 1);
        dest = V(E(i,2),:);
        rep_end = repmat(dest, num_segments, 1);

        rep_slope = rep_end - rep_origin;
%         rep_norm_slope = rep_slope ./ vecnorm(rep_slope,2,2);

        t = linspace(0,1,num_segments);
        st = rep_origin + rep_slope.*(transpose(t));
        st = st-(repmat(min(GV),num_segments,1));
        %st = unique(st,'rows');
        st = round(st / voxel_width);
        ind = st(:,1) + st(:,2)*side(1) + st(:,3)*side(1)*side(2);
        ind = ind(ind > 0);
        h(i) = -sum(H(ind+1))/size(st,1);
    end

%     for i=1:size(E,1)
%         for t=1:floor(l(i)/half_voxel_width)
%             o = V(E(i,1),:);
%             st = o+t*(half_voxel_width)*(V(E(i,2),:)-o)
%             st = ceil(st-min(GV))
%             ind = st(1) + st(2)*side(1) + st(3)*side(1)*side(2);
%             h(i) = h(i) + H(ind);
%         end
%     end
    
end