% input:
% objs - a cell array of the form 
%    objs{mesh_number}{1} are the vertices
%    objs{mesh_number}{2} are the faces
% output:
% AV, AF - the vertices and faces
function [AV,AF] = list_to_mesh(objs)

AV=[];
AF=[];

for m = 1:length(objs)
    Vm = objs{m}{1};
    Fm = objs{m}{2};
    AF = [AF; Fm+size(AV,1)];
    AV = [AV; Vm];
end

end