function make_quadmap_for_astrometry(flight,inst)
mypaths=get_paths(flight);
savedir=(strcat(mypaths.ciberdir,'doc/20170617_Stacking/maps/astroutputs/inst',...
    num2str(inst),'/'));

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

for ifield = 5%4:8
dt=get_dark_times(flight,inst,ifield);
cbmap = stackmapdat(ifield).cbmap;
mask = stackmapdat(ifield).mask_inst_clip;
map = cbmap.*mask;

for quad=['A','B','C','D']
    if quad == 'A'
        qmap = map(1:512,1:512);
    elseif quad == 'B'
        qmap = map(513:1024,1:512);
    elseif quad == 'C'
        qmap = map(1:512,513:1024);
    elseif quad == 'D'
        qmap = map(513:1024,513:1024);
    end
    fits_write(strcat(savedir,dt.name,'_',quad,'_astr.fits'),qmap);
end

end

return
