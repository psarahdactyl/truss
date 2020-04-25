function Zint = groundstructure_visibility(P,VV,FF,XX,XE,varargin)
  % Zint = groundstructure_visibility(P,VV,FF,XX,XE,varargin)
  %
  % Inputs:
  %   P  #P by 3 list of view (samples)
  %   VV  #VV by 3 list of mesh vertex positions
  %   FF  #FF by 3 list of triangle indices into VV
  %   XX  #XX by 3 list of ground structure positions
  %   XE  #XE by 2 list of indices into XX
  %   Optional:
  %     'SampleSize' followed by length of subsample on groundstructure edge
  %       {0.1m}
  %     'Weight'  followed by #P list of weights (should sum to 1) {1/#P}
  % Outputs:
  %   Zint  #XE list of integrated visibility values
  %

  function [Zs] = many_ray_mesh_intersect(P,W,SBC,VV,FF,th,SI)
    % Attempt to keep things in memory
    maxsbc = 500000;
    if size(SBC,1)>maxsbc
      Zs = [
        many_ray_mesh_intersect(P,W,SBC(1:maxsbc,:),VV,FF,th,SI(1:maxsbc)); ...
        many_ray_mesh_intersect(P,W,SBC(maxsbc+1:end,:),VV,FF,th,SI(maxsbc+1:end)); ...
        ];
      return;
    end
    RO = repmat(permute(P,[3 1 2]),[size(SBC,1) 1 1]);
    RD =        permute(SBC,[1 3 2])-permute(P,[3 1 2]);
    [H,T] = ray_mesh_intersect(reshape(RO,[],3),reshape(RD,[],3),VV,FF);
    H = H & T<0.999999;
    Zs = 1*~H;
    Zs = reshape(Zs,size(SBC,1),size(W,1));
    Zs=(Zs.*th(SI,:))*W/sum(W);
  end

  W = [];
  r = 0.1;
  vis_method = 'exact-sampling';
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Method','SampleSize','Weight',}, ...
    {'vis_method','r','W'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end
  if isempty(W)
    W = repmat(1/size(P,1),size(P,1),1);
  end

  tic;
  [SX,SE,SI] = upsample(XX,XE, ...
    'Iterations',5,'OnlySelected',@(V,E) normrow(V(E(:,2),:)-V(E(:,1),:))>r);
  SBC = barycenter(SX,SE);
  
  %Zint = zeros(size(XE,1),1);
  %for p = 1:size(P,1)
  %  Pp = P(p,:);
  %  switch vis_method
  %  case 'approx-sampling'
  %    side_x = 10*(max(VV(:,1))-min(VV(:,1)))/r;
  %    side_x
  %    [GX,side] = voxel_grid(VV,side_x,'Pad',1);
  %    % precompute Z on grid
  %    [H,T] = ray_mesh_intersect(repmat(Pp,size(GX,1),1),GX-Pp,VV,FF);
  %    H = H & T<0.999999;
  %    GZ = 1*~H;
  %    GX = reshape(GX,[side([2 3 1]) 3]);
  %    Zs = interp3( ...
  %      GX(:,:,:,1),GX(:,:,:,2),GX(:,:,:,3),reshape(GZ,size(GX(:,:,:,1))), ...
  %      SBC(:,1),SBC(:,2),SBC(:,3));
  %  case 'exact-sampling'
  %    [H,T] = ray_mesh_intersect(repmat(Pp,size(SBC,1),1),SBC-Pp,VV,FF);
  %    H = H & T<0.999999;
  %    Zs = 1*~H;
  %    Zs = W(p)/sum(W)*Zs;
  %  end
  %  % integrated viZibility  
  %  Zavg = sparse(SI,1,Zs,size(XE,1),1)./sparse(SI,1,1,size(XE,1),1);
  %  Xl = normrow(XX(XE(:,2),:)-XX(XE(:,1),:));
  %  Zintp = Xl.*Zavg;
  %  Zint = Zint+Zintp;
  %end


  D2 = permute(XX(XE(:,2),:),[1 3 2])-permute(P,[3 1 2]);
  D1 = permute(XX(XE(:,1),:),[1 3 2])-permute(P,[3 1 2]);
  th = acos(sum(D1.*D2,3)./(sqrt(sum(D1.^2,3)).*sqrt(sum(D2.^2,3))));

  switch vis_method
  %% Commented out because it's not accounting for th
  %case 'approx-sampling'
  %  side_x = 10*(max(VV(:,1))-min(VV(:,1)))/r;
  %  % make tighter around SBC instead of VV
  %  [GX,side] = voxel_grid(VV,side_x,'Pad',1);
  %  % precompute Z on grid
  %  %GZ = many_ray_mesh_intersect(P,W,GX,VV,FF);
  %  GZ = 0;
  %  tic;
  %  for p = 1:size(P,1)
  %    GZp = many_ray_mesh_intersect(P(p,:),1,GX,VV,FF);
  %    GZ = GZ + GZp * W(p)/sum(W);
  %  end
  %  toc
  %  GX = reshape(GX,[side([2 3 1]) 3]);
  %  tic;
  %  Zs = interp3( ...
  %    GX(:,:,:,1),GX(:,:,:,2),GX(:,:,:,3),reshape(GZ,size(GX(:,:,:,1))), ...
  %    SBC(:,1),SBC(:,2),SBC(:,3));
  %  toc
  %  Zs = W(1)*Zs;
  %  Zs = Zs/sum(W);
  case 'exact-sampling'
    Zs = many_ray_mesh_intersect(P,W,SBC,VV,FF,th,SI);
  end
  % integrated viZibility  
  Zavg = accumarray(SI,Zs,[size(XE,1) 1])./accumarray(SI,1,[size(XE,1) 1]);
  %Xl = normrow(XX(XE(:,2),:)-XX(XE(:,1),:));
  Zint = Zavg;
  fprintf('visibility: %g secs\n',toc);
end
