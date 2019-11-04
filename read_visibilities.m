function visibilities = read_visibilities(filename)
    fid = fopen(filename,'r');
    txt = textscan(fid,'%s','delimiter',"\n");
    fclose(fid);
    
    visibilities = {};
    for i=1:size(txt{1},1)
        W = readDMAT("vis_object_"+(i-1)+".dmat");
        visibilities{i} = W;
    end
end