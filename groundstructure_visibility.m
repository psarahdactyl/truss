% what i need passed in from my C++ tool:

% the visibility values -- can do this by writing out dmat of GV (voxel grid
% center positions) and Ss (the visibility values for each voxel)

% viable attachment points on back of objects -- can do this by finding
% viable voxels for each object and writing out dmat of indices into GV for
% each object

% things i can do in matlab:
% V,E for every object in scene and its position in the scene -- can do
% this by reading in the stls

% forces on each object (gravity on top center point of object)
% create all edges
% ground structure
% visualize

% integrate visibilities along edges
% use visibility score in LP formulation

clf

filename = ...
    "~/Documents/visibility-app/viewer/data/dogfight/dogfight_scene.txt";

%     "~/Documents/visibility-app/viewer/data/planets/planet_scene.txt";
%      "~/Documents/visibility-app/viewer/data/dale-test-scene/dale-scene.txt";
%     "~/Documents/visibility-app/viewer/data/canonical/canonical_scene.txt";




% read scene and plot the objects
[objs,bb] = read_scene(filename);
% read grid
GV = readDMAT("grid.dmat");
voxel_width = GV(1,1) - GV(2,1);
GV = GV - abs(repmat(repelem(voxel_width/2,3),size(GV,1),1));
% read visibility values of GV for each object
Vs = read_visibilities(filename);
% read surface voxels of each object
SVs = read_surface_voxels(filename);

hold on
axis equal
% plot background wall
tsurf(objs{1}{2},objs{1}{1},...
    'FaceColor',0.5+0.5*blue,falpha(1,0),fsoft)

% get all vertices/faces in the scene
% AV = objs{1}{1};
% AF = objs{1}{2};
% for i=2:size(objs,2)
%    AV = [AV; objs{i}{1}];
%    AF = [AF; objs{i}{2}+max(max(objs{i-1}{2}))];
% end

% assumes "wall" shape is the first object, as does the c++ code
for i=2:size(objs,2)
    MC = repmat(0.8*[0.99 1 0.99],size(objs{i}{2},1),1); % colors
    % plot object
    tsurf(objs{i}{2},objs{i}{1},...
        'FaceVertexCData',MC,falpha(0.8,0),fsoft);
    % find centroid of object on which to apply force (gravity)
    cm = centroid(objs{i}{1},objs{i}{2});

    DV = GV((SVs{1}==-1),:); % wall mesh

    RV = GV((SVs{i}==-1),:); % object mesh
    RVI = snap_points(RV,objs{i}{1});
    RVV = objs{i}{1}(RVI,:);
    
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
    [I,J] = find(ones(size(RV,1),size(DV,1)));

    % find centroid of object on which to apply force (gravity)
    V = [cm;RVV;DV];
    E = [ones(size(RVV,1),1) 1+(1:size(RVV,1))';1+[I size(RVV,1)+J]];

    % nodal forces, boundary vertices, stress limits
    f = zeros(size(V));
    bl = snap_points(cm,V);
    f(bl,2) = -9.8;
    bf = 1+size(RVV,1)+(1:size(DV,1));
    sC = 5e2;
    sT = 5e2;

    % get areas, axial forces, lengths, nodal equilibrium matrix
    [a,n,l,h,BT] = groundstructure(V,E,H,GV,f,bf,sC,sT,'IgnoredEdges',1:size(RVV,1));

    naa = size(RVV,1);
    av = 1 + (1:naa)';
    ae = (1:naa)';
    dim = size(V,2);
    fa = reshape(BT(av + size(V,1)*[0:dim-1],naa+1:end)*n(naa+1:end),[],dim);
    torque = normrow(sum(cross(fa,V(E(ae,2),:)-V(E(ae,1),:),2)))

    % find areas bigger than a certain threshold and place a cylinder there
    NZ = find(max(a,0)>1e-4);
    num_rods_before_cluster = size(NZ,1)
    num_compression = sum(sign(n(NZ))==1)
    num_tension = sum(sign(n(NZ))==-1)
    [CV,CF,CJ,CI] = edge_cylinders(V,E(NZ,:),...
        'PolySize',10,'Thickness',sqrt(max(a(NZ),0)/pi));

    fbf = reshape(BT(bf(:) + [0:dim-1]*size(V,1),:)*n,[],dim);

    % plot the cylinders
%     color = sign(n(NZ(CJ)));
%     tsurf(CF,CV,falpha(1,0),'CData',color,fsoft);
    
    % plot clustered cylinders
    % ideal k (last argument) is the "ideal" (by user standards)
    % number of RODS for tension
    % and the ideal number of RODS for COMPRESSION
    % i.e. if k = 3, you are asking for 3 compression rods and 3 tension
    % rods
    % nevermind, k is now how many TOTAL rods you want to end up with
    [nV, nE, nNZ, nn, na] = cluster_endpoints(V,E,NZ,n,H,GV,sC,sT,cm,7);
    [nCV,nCF,nCJ,nCI] = edge_cylinders(nV,nE(nNZ,:),...
        'PolySize',10,'Thickness',sqrt(max(na(nNZ),0)/pi));
    
    % plot the new cylinders
    color = sign(nn(nNZ(nCJ)));
    tsurf(nCF,nCV,falpha(1,0),'CData',color,fsoft);
    
    % plot COM
    scatter3(cm(:,1),cm(:,2),cm(:,3),'.r','SizeData',1000);
    
    % plot visibilities
%     GV_h = GV(H~=0,:);
%     H_h = H(H~=0);
%     scatter3(GV_h(:,1),GV_h(:,2),GV_h(:,3),...
%         '.','CData', abs(H_h),'SizeData',300);

    quiver3( ...
      V(:,1),V(:,2),V(:,3), ...
      f(:,1),f(:,2),f(:,3),'r','LineWidth',2);
  
      % plot isosurface
      % TODO: make this nicer than hardcoding
      % aka turn GV into a meshgrid somehow
%     [X,Y,Z] = meshgrid(-7:.73:7,-7:.73:7,-7:.73:7);
%     H_h = reshape(H,20,20,20);
%     H_r = fliplr(H_h);
%     H_r = imrotate3(H_r,90,[0 0 1]);
% 
%     [faces,verts,colors] = isosurface(X,Y,Z,H_r,-.5,abs(H_r));
%     
%     patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
%     'FaceColor','interp','EdgeColor','none')
%     alpha(0.3);
  

% scatter3(RVV(:,1),RVV(:,2),RVV(:,3),'.b','SizeData',1000);

% quiver3( ...
%   V(bf,1),V(bf,2),V(bf,3), ...
%   fbf(:,1),fbf(:,2),fbf(:,3), ...
%   0,'.b','LineWidth',3);
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
% colormap([cbrewer('RdBu',256);cbrewer('PiGn',256)])
% colormap([parula(256);jet(256)])
% colormap(flipud(cbrewer('RdBu',256)))
colormap(flipud(cbrewer('PiYG',256)))
colorbar

