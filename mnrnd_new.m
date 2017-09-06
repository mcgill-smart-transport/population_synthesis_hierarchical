function z = mnrnd_new(prob_mat)

[number,dim] = size(prob_mat);
cump = [zeros(number,1),cumsum(prob_mat,2)];
uf = rand(number,1);
z = zeros(number,dim);
for f = 1:dim
    ind = uf > cump(:,f) & uf <= cump(:,f+1);
    z(ind,f) = 1;
end

