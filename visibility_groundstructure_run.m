clf

filename = ...
     "~/Documents/siggraph2020/visibility-app/viewer/meshes/canonical/canonical_scene.txt";

% read scene and plot the objects
[objs,bb] = read_scene(filename);
% read grid
GV = readDMAT("grid.dmat");
Vs = readDMAT("grid_visibilities.dmat");
voxel_width = GV(1,1) - GV(2,1);
% GV = GV - abs(repmat(repelem(voxel_width/2,3),size(GV,1),1));


hold on;
% axis equal tight manual % this ensures that getframe() returns a consistent size
% camlight;
% camup([0 1 0]);
% camproj('perspective')
% view(26,15)

% gif_filename = 'k-means-rods.gif';
    
% plot background wall
wall = tsurf(objs{1}{2},objs{1}{1},...
    'FaceColor',0.5+0.5*blue,falpha(1,0),fsoft);

% get all vertices/faces in the scene
% AV = objs{1}{1};
% AF = objs{1}{2};
% for i=2:size(objs,2)
%    AV = [AV; objs{i}{1}];
%    AF = [AF; objs{i}{2}+max(max(objs{i-1}{2}))];
% end

rng(0);
writeSTL("./mesh0.stl",objs{1}{1},objs{1}{2});

% assumes "wall" shape is the first object, as does the c++ code
for i=2:size(objs,2)
    MC = repmat(0.8*[0.99 1 0.99],size(objs{i}{2},1),1); % colors
    % plot object
    object = tsurf(objs{i}{2},objs{i}{1},...
        'FaceVertexCData',MC,falpha(0.8,0),fsoft);
    
    % find centroid of object on which to apply force (gravity)
    cm = centroid(objs{i}{1},objs{i}{2});

%     DV = GV((SVs{1}==-1),:); % wall mesh
% 
%     RV = GV((SVs{i}==-1),:); % object mesh
%     RVI = snap_points(RV,objs{i}{1});
%     RVV = objs{i}{1}(RVI,:);
    
    % find points on other objects in scene
%     AV = GV(GV(:,3)<min(RV(:,3)),:);
%     AF = zeros(size(AV,1),3);
%     [~,T] = ray_mesh_intersect(...
%            repmat(RV,size(AV,1),1),...
%            repmat(AV,size(RV,1),1) - repelem(RV,size(AV,1),1),...
%            AV,AF);
%     % RV is the points on hammer
%     DV = GV(~isinf(T)&T~=0,:)+T(~isinf(T)&T~=0);%.*repmat(RV,size(AV,1),1)-repmat(RV,size(AV,1),1);

    % hidden-ness/visibility
    H = Vs{i};

    % this creates the ground structure, (V,E), including the center of mass
    % and connection points on both wall and hammer

    RVV = readDMAT("set1_"+(i-1)+".dmat");
    DV = readDMAT("set2_"+(i-1)+".dmat");
    [I,J] = find(ones(size(RVV,1),size(DV,1)));
    % find centroid of object on which to apply force (gravity)
%     V = [cm;RVV;DV];
%     E = [ones(size(RVV,1),1) 1+(1:size(RVV,1))';1+[I size(RVV,1)+J]];

%     V = mmread('../visibility-app/viewer/build/V_gsm1.mtx');
%     E = mmread('../visibility-app/viewer/build/E_gsm1.mtx');
%     V = full(V);
%     E = full(E);
%     E = E+1;
% 
%     % nodal forces, boundary vertices, stress limits
%     f = zeros(size(V));
%     bl = snap_points(cm,V);
%     f(bl,2) = -9.8;
% %     bf = 1+size(RVV,1)+(1:size(DV,1));
%     bf = mmread('../visibility-app/viewer/build/bf1.mtx');
%     bf = full(bf)+1;
%     sC = 5e3;
%     sT = 5e3;
% 
%     fe = mmread('../visibility-app/viewer/build/fe1.mtx');
%     fe = full(fe)+1;
    % get areas, axial forces, lengths, nodal equilibrium matrix
%     [a,n,l,h,BT] = groundstructure(V,E,H,GV,f,bf,sC,sT,'IgnoredEdges',fe);
%     m = size(E,1);
%     an = mmread('../visibility-app/viewer/build/solution1.mtx');
%     an = full(an);
%     a = an(1:m);
%     nf = 1;
%     n = reshape(an(m+(1:m*nf)),m,nf);

%     naa = size(RVV,1);
%     av = 1 + (1:naa)';
%     ae = (1:naa)';
%     dim = size(V,2);
%     fa = reshape(BT(av + size(V,1)*[0:dim-1],naa+1:end)*n(naa+1:end),[],dim);
%     torque = normrow(sum(cross(fa,V(E(ae,2),:)-V(E(ae,1),:),2)))

    V = readDMAT('V.dmat');
    E = readDMAT('E.dmat');
    a = readDMAT("areas.dmat");
    % find areas bigger than a certain threshold and place a cylinder there
    NZ = find(max(a,0)>1e-4);
    num_rods_before_cluster = size(NZ,1)
%     num_compression = sum(sign(n(NZ))==1)
%     num_tension = sum(sign(n(NZ))==-1)
    [CV,CF,CJ,CI] = edge_cylinders(V,E(NZ,:),...
        'PolySize',10,'Thickness',sqrt(max(a(NZ),0)/pi));

%     fbf = reshape(BT(bf(:) + [0:dim-1]*size(V,1),:)*n,[],dim);
%     plot_edges(V,E(NZ,:),'Color',[0.9 0.1840 0.5560],'LineWidth',3);

    % plot the cylinders
%     color = sign(n(NZ(CJ)));

    cylinders = tsurf(CF,CV);%,falpha(1,0),'CData',color,fsoft);
%     hold off;
%     caxis([-1 1])
%     colormap(flipud(cbrewer('RdBu',256)))
%     frame = getframe; 
%     im = frame2im(frame); 
%     [imind,colorm] = rgb2ind(im,256);
%     imwrite(imind,colorm,gif_filename,'gif', 'Loopcount',inf); 
    
%     writeSTL("./mesh"+(i-1)+".stl",objs{i}{1},objs{i}{2});
    
%     % plot clustered cylinders
%     % ideal k (last argument) is the "ideal" (by user standards)
%     % k is how many TOTAL rods you want to end up with
%     num_rods_collection = ones(size(NZ,1),1);
%     num_rods_collection(1) = size(NZ,1);
%     for fr = 1:size(NZ,1)-1
%         hold on
%         ideal_k = size(NZ,1)
%         [SV, SE, nNZ, nn, na] = cluster_endpoints(V,E,NZ,n,H,GV,sC,sT,cm,ideal_k);
%         [nCV,nCF,nCJ,nCI] = edge_cylinders(SV,SE(nNZ,:),...
%             'PolySize',10,'Thickness',sqrt(max(na(nNZ),0)/pi));
%         
%         num_rods_collection(fr+1) = size(nNZ,1);
%         num_rods_after_cluster = size(nNZ,1)
% 
% 
%         % plot the new cylinders
%         if (size(nNZ,1) ~= size(NZ,1) && size(nNZ,1)~=0)
%             delete(cylinders)
%             color = sign(nn(nNZ(nCJ)));
%             cylinders = tsurf(nCF,nCV,falpha(1,0),'CData',color,fsoft);
%             writeSTL("./cylinders"+(i-1)+".stl",CV,CF);
% 
% %             plot COM
%             scatter3(cm(:,1),cm(:,2),cm(:,3),'.r','SizeData',1000);   
% 
%             % Capture the plot as an image 
%             frame = getframe; 
%             im = frame2im(frame); 
%             [imind,colorm] = rgb2ind(im,256); 
%             % Write to the GIF File 
%             imwrite(imind,colorm,gif_filename,'gif','WriteMode','append'); 
%         end
%     
%     end
    
%     % plot visibilities
%     GV_h = GV(H~=0,:);
%     H_h = H(H~=0);
%     scatter3(GV_h(:,1),GV_h(:,2),GV_h(:,3),...
%         '.','CData', abs(H_h),'SizeData',300);

%     quiver3( ...
%       V(:,1),V(:,2),V(:,3), ...
%       f(:,1),f(:,2),f(:,3),'r','LineWidth',2);
  
      % plot isosurface
      % TODO: make this nicer than hardcoding
      % aka turn GV into a meshgrid somehow
%     [X,Y,Z] = meshgrid(-7:.73:7,-7:.73:7,-7:.73:7);
%     H_h = reshape(H,20,20,15);
%     H_r = fliplr(H_h);
%     H_r = imrotate3(H_r,90,[0 0 1]);

%     [faces,verts,colors] = isosurface(X,Y,Z,H_r,-.5,abs(H_r));
%     
%     patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
%     'FaceColor','interp','EdgeColor','none')
%     alpha(0.3);

end

hold off;
axis equal;
view(26,15);
camlight;
camup([0 1 0]);
caxis(max(abs(caxis))*[-1 1])
cameratoolbar('SetCoordSys','y')
cameratoolbar('setmode','orbit')
camproj('perspective')
caxis([-1 1])
colormap(flipud(cbrewer('RdBu',256)))
colorbar

