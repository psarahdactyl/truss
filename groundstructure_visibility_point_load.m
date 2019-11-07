% this file will be using point loads instead of rigid body loads
clf
filename = ...
    "~/Documents/visibility-app/viewer/data/canonical/canonical_scene.txt";


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

% assumes "wall" shape is the first object, as does the c++ code
for i=2:size(objs,2)
    MC = repmat(0.8*[0.99 1 0.99],size(objs{i}{2},1),1); % colors
    % plot object
    tsurf(objs{i}{2},objs{i}{1},...
        'FaceVertexCData',MC,falpha(0.8,0),fsoft);

    DV = GV((SVs{1}==-1),:); % wall mesh

    RV = GV((SVs{i}==-1),:); % object mesh
    RVI = snap_points(RV,objs{i}{1});
    RVV = objs{i}{1}(RVI,:);

    % hidden-ness/visibility
    H = Vs{i};

    % this creates the ground structure, (V,E), including the center of mass
    % and connection points on both wall and hammer
    [I,J] = find(ones(size(RV,1),size(DV,1)));

    % find centroid of object on which to apply force (gravity)
    V = [RVV;DV];
    E = [I size(RVV,1)+J];

    % nodal forces, boundary vertices, stress limits
    % place nodal force on first V
    f = zeros(size(V));
    f(1:3,2) = -9.8;
    V_forces = V(1:3,:);
    bf = size(RVV,1)+(1:size(DV,1));
    sC = 5e2;
    sT = 5e2;

    % get areas, axial forces, lengths, nodal equilibrium matrix
    [a,n,l,h,BT] = groundstructure(V,E,H,GV,f,bf,sC,sT);

    % find areas bigger than a certain threshold and place a cylinder there
    NZ = find(max(a,0)>1e-4);
    num_rods_before_cluster = size(NZ,1)
    num_compression = sum(sign(n(NZ))==1)
    num_tension = sum(sign(n(NZ))==-1)
    [CV,CF,CJ,CI] = edge_cylinders(V,E(NZ,:),...
        'PolySize',10,'Thickness',sqrt(max(a(NZ),0)/pi));
    
    % plot clustered cylinders
    % ideal k (last argument) is the "ideal" (by user standards)
    % k is how many TOTAL rods you want to end up with
    [nV, nE, nNZ, nn, na] = cluster_endpoints_point_load(V,E,NZ,n,H,GV,sC,sT,size(NZ,1),V_forces);
    [nCV,nCF,nCJ,nCI] = edge_cylinders(nV,nE(nNZ,:),...
        'PolySize',10,'Thickness',sqrt(max(na(nNZ),0)/pi));
    
    % plot the new cylinders
    color = sign(nn(nNZ(nCJ)));
    tsurf(nCF,nCV,falpha(1,0),'CData',color,fsoft);
    
    % plot visibilities
%     GV_h = GV(H~=0,:);
%     H_h = H(H~=0);
%     scatter3(GV_h(:,1),GV_h(:,2),GV_h(:,3),...
%         '.','CData', abs(H_h),'SizeData',300);

    quiver3( ...
      V(:,1),V(:,2),V(:,3), ...
      f(:,1),f(:,2),f(:,3),'r','LineWidth',2);
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

