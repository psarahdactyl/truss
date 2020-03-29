function [objs, bb] = read_scene(filename)
    % open scene file
    fid = fopen(filename,'r');
    txt = textscan(fid,'%s','delimiter',"\n");
    fclose(fid);
    
    file_arr = split(filename,"/");
    just_path = file_arr(1:end-1);
    path = join(just_path,"/");
    
    % read STLs
    objs = {};
    mm = zeros(size(txt{1},1),3);
    for i=1:size(txt{1},1)
        p =path+"/"+txt{1}{i};
        p
        [V,F] = readSTL(p);
        [SV,SVI,SVJ] = remove_duplicate_vertices(V,1e-7);
        SF = SVJ(F);
        mm(i,:) = min(SV);
        mm(i+size(txt{1},1),:) = max(SV);
        objs{i} = {SV,SF};
        size(SV)
    end
    
    % translate STLs
    [bb,~] = bounding_box(mm);
    
%     for i=1:size(objs,2)
%         objs{i}{1} = objs{i}{1}...
%             - repmat(((max(bb)+min(bb))/2),size(objs{i}{1},1),1);
%     end

end