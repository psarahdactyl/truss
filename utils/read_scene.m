function [new_objs, bb, cams] = read_scene(filename)
    % open scene file
    fid = fopen(filename,'r');
    txt = textscan(fid,'%s','delimiter',"\n");
    fclose(fid);
    
    file_arr = split(filename,"/");
    just_path = file_arr(1:end-1);
    path = join(just_path,"/");
    
    % read STLs
    objs = {}; %cell(size(txt{1},1),3);
    mm = []; %zeros(size(txt{1},1),3);
    cams = []; %zeros(size(txt{1},1),3);
    oi=1;
    ci=1;
    for i=1:size(txt{1},1)
        a = split(txt{1}{i});
        if(size(a,1)==2) 
          p =path+"/"+a{1};
          [V,F] = readSTL(p);
  %         [V,F] = readSTL(a{1});
          [SV,~,SVJ] = remove_duplicate_vertices(V,1e-7);
          SF = SVJ(F);
          mm(oi,:) = min(SV);
          mm(oi+size(txt{1},1),:) = max(SV);
          objs{oi,1} = SV;
          objs{oi,2} = SF;
          objs{oi,3} = a{2};
          oi = oi+1;
        else
          cam=str2double(a)
          cams(ci,:) = cam;
          ci = ci+1;
        end
    end
    
    [bb,~] = bounding_box(mm);
    
    % find fixed object (will have a 1 in the third col)
    fi=find(str2double(objs(:,3))==1)
    % reorder so first object is fixed
    new_objs = [objs(fi,:); objs(1:fi-1,:); objs(fi:end-1,:)]
  
    
    
% %     translate STLs
%     for i=1:size(objs,2)
%         objs{i}{1} = objs{i}{1}...
%             - repmat(((max(bb)+min(bb))/2),size(objs{i}{1},1),1);
%     end

end