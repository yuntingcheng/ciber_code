function run_ihl_stack_sides(flight,inst,ifield)
mypaths=get_paths(flight);
savedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/sides/'));
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat_sides',loaddir),'stackmapdat');

dx = 1200;
verbose = 1;
dt=get_dark_times(flight,inst,ifield);
cbmap = stackmapdat(ifield).cbmap;
psmap = stackmapdat(ifield).psmap;
mask_inst_clip = ones(size(cbmap));
strmask = stackmapdat(ifield).strmask;
strnum = stackmapdat(ifield).strnum;
m_min_arr = stackmapdat(ifield).m_min_arr;
m_max_arr = stackmapdat(ifield).m_max_arr;
for im=1:numel(m_min_arr)
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
 
    [stampercb,stamperps,maskstamper,count_tot,count_stack]=...
    stackihl_sides(flight,inst,ifield,m_min,m_max,dx,...
    cbmap,psmap,mask_inst_clip,strmask,strnum,verbose);
     
    fits_write(strcat(savedir,dt.name,'_stampercb',...
        num2str(m_min),'_',num2str(m_max),'.fits'),stampercb);
    fits_write(strcat(savedir,dt.name,'_hitmap',...
        num2str(m_min),'_',num2str(m_max),'.fits'),maskstamper);
    fits_write(strcat(savedir,dt.name,'_stamperps',...
        num2str(m_min),'_',num2str(m_max),'.fits'),stamperps);
end

return