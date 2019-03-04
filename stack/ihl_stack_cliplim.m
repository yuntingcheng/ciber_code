function cliplim_arr = ihl_stack_cliplim(flight,inst,ifield,type,...
    m_max,m_min,dx,nbins,sig,rmin)

mypaths=get_paths(flight);
savedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));
dt=get_dark_times(flight,inst,ifield);

if type==-1
    datacb=fitsread(strcat(savedir,'ciber_ps/',dt.name,...
        '_stampers',num2str(m_min),'_',num2str(m_max),'_1.fits'));
    dataps=fitsread(strcat(savedir,'panstarrs_ps/',dt.name,...
        '_stampers',num2str(m_min),'_',num2str(m_max),'_1.fits'));
    mdata=load(strcat(savedir,'ciber_ps/',dt.name,...
        '_hitmaps',num2str(m_min),'_',num2str(m_max),'_1.mat'));
    mdata = mdata.hitmap;
    stackcb=datacb./mdata;
    stackps=dataps./mdata;
elseif type==1
    datacb=fitsread(strcat(savedir,'ciber_ps/',dt.name,...
        '_stamperg',num2str(m_min),'_',num2str(m_max),'_1.fits'));
    dataps=fitsread(strcat(savedir,'panstarrs_ps/',dt.name,...
        '_stamperg',num2str(m_min),'_',num2str(m_max),'_1.fits'));
    mdata=load(strcat(savedir,'ciber_ps/',dt.name,...
        '_hitmapg',num2str(m_min),'_',num2str(m_max),'_1.mat'));
    mdata = mdata.hitmap;
    stackcb=datacb./mdata;
    stackps=dataps./mdata;
end

rad = make_radius_map(stackcb,dx+1,dx+1);
binedge=logspace(log10(rmin),log10(max(rad(:))),nbins+1);
if rmin==0 ; binedge(1)=0 ; end

maskcb = ones(size(stackcb));
maskcb(find(stackcb==0 | stackcb~=stackcb)) = 0;
maskps = ones(size(stackps));
maskps(find(stackps==0 | stackps~=stackps)) = 0;

cliplim_arr = zeros([4,nbins]);
for ibin=1:nbins
    spcb=find(rad >= binedge(ibin) ...
        & rad < binedge(ibin+1) & maskcb);
    spps=find(rad >= binedge(ibin) ...
        & rad < binedge(ibin+1) & maskps);
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

return