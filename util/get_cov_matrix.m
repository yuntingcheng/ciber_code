function cov = get_cov_matrix(data)
%%% given a data with size (Nsample, Nbin) return Nbin x Nbin Cov

Nbin = size(data,2);

if Nbin>1
    cov = zeros(Nbin);
    for i=1:Nbin
        for j=1:Nbin
            datai = data(:,i);
            dataj = data(:,j);
            cov(i,j) = mean(datai.*dataj) - mean(datai)*mean(dataj);
        end
    end
else
    cov = mean(data.^2) - mean(data)^2;
end


return