function run_psf_stack(flight,inst,ifield)
mypaths=get_paths(flight);
savedir=(strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf/TM',...
    num2str(inst),'/'));

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

dx = 1200;
verbose = 1;
dt=get_dark_times(flight,inst,ifield);
cbmap = stackmapdat(ifield).cbmap;
mask_inst_clip = stackmapdat(ifield).mask_inst_clip;
[tmmask,tmnum] = make_mask_2m(flight,inst,'y',ifield,0,18);
strmask = tmmask;
strnum = tmnum;

m_min = 0;
m_max = 14;

[stamper,hitmap]=stackpsf_2m(flight,inst,ifield,m_min,m_max,dx,...
    cbmap,mask_inst_clip,strmask,strnum,verbose);

fits_write(strcat(savedir,dt.name,'_stamper.fits'),stamper);
fits_write(strcat(savedir,dt.name,'_hitmap.fits'),hitmap);

return