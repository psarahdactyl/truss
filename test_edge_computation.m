% a small example to test vectorized(-ish) edge visibility computation

v = [0,0,0;1,3,1;4,1,2;3,1,1];
e = [1,2;3,4];

w = 0.1;

os = v(e(:,1),:);
ds = v(e(:,2),:);

l = vecnorm(ds-os,2,2);
ss = ceil(l / w);

oss = repelem(os,max(ss),1);
dss = repelem(ds,max(ss),1);

sl = dss-oss;

ts = zeros(size(ss,1),max(ss));

% yuck
for i=1:size(ss,1)
    p = linspace(0,1,ss(i));
    ts(i,:) = [p zeros(1,max(ss)-size(p,2))];
end

size(ts)

tss = reshape(ts',size(ts,1)*size(ts,2),1);

sts = oss + sl.*tss;
size(sts)