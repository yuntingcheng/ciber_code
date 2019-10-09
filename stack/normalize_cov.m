function cov_rho = normalize_cov(cov)
% given a cov matrix, return the correlation coeff. (normalized cov), 
% i.e. cov_rho(i,j) = cov(i,j)/sqrt(cov(i,i)cov(j,j))

N = size(cov,1);
cov_rho = cov;
for i=1:N
    for j=1:N
        if cov(i,i)==0 | cov(j,j)==0
            cov_rho(i,j) = cov(i,j);
        else
            cov_rho(i,j) = cov(i,j)./sqrt(cov(i,i)*cov(j,j));
        end
    end
end

return