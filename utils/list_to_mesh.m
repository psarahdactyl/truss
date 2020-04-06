% input:
% objs - a cell array of the form 
%    objs{mesh_number}{1} are the vertices
%    objs{mesh_number}{2} are the faces
% output:
% AV, AF - the vertices and faces
function [AV,AF,ACV,ACF,coms] = list_to_mesh(objs)

AV=[];
AF=[];
ACV=[];
ACF=[];
coms = zeros(size(objs,1),3);

for m = 1:length(objs)
  Vm = objs{m,1};
  Fm = objs{m,2};
  AF = [AF; Fm+size(AV,1)];
  AV = [AV; Vm];
  ACV = [ACV; repmat(m,size(Vm,1),1)];
  ACF = [ACF; repmat(m,size(Fm,1),1)];
  % find centroid of object on which to apply force (gravity)
  com = centroid(Vm,Fm);
  coms(m,:) = com;
end

end