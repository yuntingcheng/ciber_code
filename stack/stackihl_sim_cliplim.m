function [cliplim_arrin,binedgesin,cliplim_arrout,binedgesout]=...
    stackihl_sim_cliplim(flight,inst,ifield,type,m_min,m_max,...
    cbmap,psmap,mask_inst,strmask,strnum)

nbins = 6;
sig = 10;
rmin = 10;

dx = 50;
unmask = 1;
Nsrc = 20;
verbose = 0;

[stampercb,stamperps,hitmap]=...
    stackihl_sim0(flight,inst,ifield,type,m_min,m_max,dx,...
    cbmap,psmap,mask_inst,strmask,strnum,unmask,Nsrc,verbose,0,0,[]);
stackcb=stampercb./hitmap;
stackps=stamperps./hitmap;
rad = make_radius_map(stackcb,dx+1,dx+1);
binedges=logspace(log10(rmin),log10(max(rad(:))),nbins+1);

maskcb = ones(size(stackcb));
maskcb(find(stackcb==0 | stackcb~=stackcb)) = 0;
maskps = ones(size(stackps));
maskps(find(stackps==0 | stackps~=stackps)) = 0;

cliplim_arr = zeros([4,nbins]);
for ibin=1:nbins
    spcb=find(rad >= binedges(ibin) ...
        & rad < binedges(ibin+1) & maskcb);
    spps=find(rad >= binedges(ibin) ...
        & rad < binedges(ibin+1) & maskps);
    bcb = stackcb(spcb);
    bps = stackps(spps);
    clipmaxcb = nanmedian(bcb(:))+sig*std(bcb(:));
    clipmincb = nanmedian(bcb(:))-sig*std(bcb(:));        
    clipmaxps = nanmedian(bps(:))+sig*std(bps(:));
    clipminps = nanmedian(bps(:))-sig*std(bps(:));
    cliplim_arr(1,ibin) = clipmaxcb;
    cliplim_arr(2,ibin) = clipmincb;
    cliplim_arr(3,ibin) = clipmaxps;
    cliplim_arr(4,ibin) = clipminps;
end
cliplim_arrout  = cliplim_arr;
binedgesout = binedges;

sig = 1;
binedges=logspace(log10(1),log10(rmin),nbins+1);
cliplim_arr = zeros([4,nbins]);
for ibin=1:nbins
    spcb=find(rad >= binedges(ibin) ...
        & rad < binedges(ibin+1) & maskcb);
    spps=find(rad >= binedges(ibin) ...
        & rad < binedges(ibin+1) & maskps);
    bcb = stackcb(spcb);
    bps = stackps(spps);
    clipmaxcb = nanmedian(bcb(:))+sig*std(bcb(:));
    clipmincb = nanmedian(bcb(:))-sig*std(bcb(:));        
    clipmaxps = nanmedian(bps(:))+sig*std(bps(:));
    clipminps = nanmedian(bps(:))-sig*std(bps(:));
    cliplim_arr(1,ibin) = clipmaxcb;
    cliplim_arr(2,ibin) = clipmincb;
    cliplim_arr(3,ibin) = clipmaxps;
    cliplim_arr(4,ibin) = clipminps;
end
cliplim_arrin = cliplim_arr;
binedgesin = binedges;

return