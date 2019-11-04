function SV = read_surface_voxels(filename)
    fid = fopen(filename,'r');
    txt = textscan(fid,'%s','delimiter',"\n");
    fclose(fid);
    
    SV = {};
    for i=1:size(txt{1},1)
        W = readDMAT("object_voxels_"+(i-1)+".dmat");
        SV{i} = W;
    end
end