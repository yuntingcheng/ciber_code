function bigpsf=get_coarse_map(psfmap,bigsize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make pixscale =0.7 psfmap become pixscale=7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    N=size(psfmap,1);
    finepsf=zeros(2*N);
    for i=1:2*N
        for j=1:2*N
            finepsf(i,j)=psfmap(ceil(i/2),ceil(j/2))/4;
        end
    end

    Nsubpix=2*(floor((N-bigsize)/bigsize/2)*bigsize*2+bigsize);
    Ncut=(2*N-Nsubpix)/2;
    finepsf=finepsf(Ncut+1:Ncut+Nsubpix,Ncut+1:Ncut+Nsubpix);

    bigpsf=zeros(Nsubpix/bigsize/2);
    for i=1:size(bigpsf,1)
        for j=1:size(bigpsf,2)
            submap=finepsf((i-1)*2*bigsize+1:i*2*bigsize,...
                           (j-1)*2*bigsize+1:j*2*bigsize);
            bigpsf(i,j)=sum(submap(:));
        end
    end
return