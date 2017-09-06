function out=sumdims(r,dims)

out=r;

for i=1:length(dims)
    out=sum(out,dims(i));
end

out = squeeze(out);