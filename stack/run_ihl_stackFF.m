function run_ihl_stackFF(flight,inst,ifield)
mypaths=get_paths(flight);
savedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

dx = 1200;
verbose = 1;
dt=get_dark_times(flight,inst,ifield);
psmap = stackmapdat(ifield).psmap;
psmapFF = stackmapdat(ifield).psmap_FFerr;
mask_inst = stackmapdat(ifield).mask_inst_clip;
strmask = stackmapdat(ifield).strmask;
strnum = stackmapdat(ifield).strnum;
m_min_arr = stackmapdat(ifield).m_min_arr;
m_max_arr = stackmapdat(ifield).m_max_arr;

nbins = 25;
sig = 5;
rmin = 20;

for im=10:12%1:numel(m_min_arr)

    m_min = m_min_arr(im);
    m_max = m_max_arr(im);

    type = 1;
    cliplim_arr = ihl_stack_cliplim(flight,inst,ifield,type,...
    m_max,m_min,dx,nbins,sig,rmin);
    [stamperps,stamperpsFF,hitmap,count_totg,count_stackg]=...
        stackihl_psFF(flight,inst,ifield,type,m_min,m_max,dx,...
        psmap,psmapFF,mask_inst,strmask,strnum,cliplim_arr,verbose); 
    
    fits_write(strcat(savedir,'FFerr_test/',dt.name,'_stamperg',...
        num2str(m_min),'_',num2str(m_max),'.fits'),stamperps);
    save(strcat(savedir,'FFerr_test/',dt.name,'_hitmapg',...
        num2str(m_min),'_',num2str(m_max),'.mat'),'hitmap');
    fits_write(strcat(savedir,'FFerr_test/',dt.name,'_stamperg',...
        num2str(m_min),'_',num2str(m_max),'FF.fits'),stamperpsFF);

    type = -1;
    cliplim_arr = ihl_stack_cliplim(flight,inst,ifield,type,...
    m_max,m_min,dx,nbins,sig,rmin);
    [stamperps,stamperpsFF,hitmap,count_totg,count_stackg]=...
        stackihl_psFF(flight,inst,ifield,type,m_min,m_max,dx,...
        psmap,psmapFF,mask_inst,strmask,strnum,cliplim_arr,verbose); 
    
    fits_write(strcat(savedir,'FFerr_test/',dt.name,'_stampers',...
        num2str(m_min),'_',num2str(m_max),'.fits'),stamperps);
    save(strcat(savedir,'FFerr_test/',dt.name,'_hitmaps',...
        num2str(m_min),'_',num2str(m_max),'.mat'),'hitmap');
    fits_write(strcat(savedir,'FFerr_test/',dt.name,'_stampers',...
        num2str(m_min),'_',num2str(m_max),'FF.fits'),stamperpsFF);

end

return