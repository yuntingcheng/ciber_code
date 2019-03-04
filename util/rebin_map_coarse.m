function rbmap=rebin_map_coarse(map,Nsub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bin the map into coarse pixel with mean 
%of NsubxNsub original pixel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nfine=size(map,1);
Ncoarse=Nfine./Nsub;

rbmap=zeros(Ncoarse);
for i=1:Ncoarse
    for j=1:Ncoarse
        rbmap(i,j)=mean(mean(map(i*Nsub-Nsub+1:i*Nsub,j*Nsub-Nsub+1:j*Nsub)));
    end
end

return