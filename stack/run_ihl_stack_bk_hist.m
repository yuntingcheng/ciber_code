function run_ihl_stack_bk_hist(flight,inst,run,spire)
mypaths=get_paths(flight);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

savedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

dx = 1200;
for ifield=8:-1:4
    dt=get_dark_times(flight,inst,ifield);
    load(sprintf('%shistdat_%s',loaddir,dt.name),'histdat');
    cbmap = stackmapdat(ifield).cbmap;
    psmap = stackmapdat(ifield).psmap;
    mask_inst_clip = stackmapdat(ifield).mask_inst_clip;
    strmask = stackmapdat(ifield).strmask;
    
    counts_arr = [];
    countg_arr = [];
    for im = 10:12
        counts_arr = [counts_arr, histdat(im).counts];
        countg_arr = [countg_arr, histdat(im).countg];
    end
    
    N_arr = unique([counts_arr, countg_arr]);
    N_arr = N_arr(find(N_arr>0));
%     N_arr = [N_arr, max(N_arr)*3];
    Ibinedges_cb = histdat(im).Ibinedges_cb;
    Ibinedges_ps = histdat(im).Ibinedges_ps;
    for iter=(0:4)*10 + run
        [histbkdat]=stackihl_ps_randomN_hist(N_arr,dx,cbmap,psmap,...
            Ibinedges_cb,Ibinedges_ps,mask_inst_clip.*strmask,spire);
        
        save(strcat(savedir,'bk_ps/',dt.name,'_histbkdat','_',num2str(iter)),...
            'histbkdat');
        
    end
end

return