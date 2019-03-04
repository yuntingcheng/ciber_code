function run_ihl_stack_mock_star(flight,inst)
mypaths=get_paths(flight);
savedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

dx = 1200;
verbose = 1;

for ifield=4:8
dt=get_dark_times(flight,inst,ifield);
cbmap = stackmapdat(ifield).cbmap;
psmap = stackmapdat(ifield).psmap;
mask_inst_clip = stackmapdat(ifield).mask_inst_clip;
strmask = stackmapdat(ifield).strmask;
m_min_arr = stackmapdat(ifield).m_min_arr(10:12);
m_max_arr = stackmapdat(ifield).m_max_arr(10:12);
counts_arr = stackmapdat(ifield).count_stacks_arr(10:12);
countg_arr = stackmapdat(ifield).count_stackg_arr(10:12);
interp = 1;
for im=find(countg_arr>counts_arr)
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    counts = counts_arr(im);
    countg = countg_arr(im);
    Nsrc = countg - counts;
    
    disp(sprintf('%s,stack %d src between %d < m < %d '...
        ,dt.name,Nsrc,m_min,m_max));

    [stampercb,stamperps,hitmap]=...
    stackihl_mock_star(flight,inst,ifield,Nsrc,m_min,m_max,interp,dx,...
    cbmap,psmap,mask_inst_clip,strmask,verbose);

    fits_write(strcat(savedir,'ciber_ps/',dt.name,'_stampers',...
        num2str(m_min),'_',num2str(m_max),'mock.fits'),stampercb);
    save(strcat(savedir,'ciber_ps/',dt.name,'_hitmaps',...
        num2str(m_min),'_',num2str(m_max),'mock.mat'),'hitmap');
    fits_write(strcat(savedir,'panstarrs_ps/',dt.name,'_stampers',...
        num2str(m_min),'_',num2str(m_max),'mock.fits'),stamperps);
end
end

return