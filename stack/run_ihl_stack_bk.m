function run_ihl_stack_bk(flight,inst,run)
mypaths=get_paths(flight);
savedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

dx = 1200;
for ifield=8:-1:4
    dt=get_dark_times(flight,inst,ifield);
    load(sprintf('%s/stackcountdat_%s',loaddir,dt.name),'stackcountdat');
    cbmap = stackmapdat(ifield).cbmap;
    psmap = stackmapdat(ifield).psmap;
    mask_inst_clip = stackmapdat(ifield).mask_inst_clip;
    strmask = stackmapdat(ifield).strmask;
    
    counts_arr = [];
    countg_arr = [];
    for im = 10:12        
        counts_arr = [counts_arr, stackcountdat(im).counts];
        countg_arr = [countg_arr, stackcountdat(im).countg];
    end
    N_arr = unique([counts_arr, countg_arr]);
    N_arr = N_arr(find(N_arr>0));
    N_arr = [N_arr, max(N_arr)*3];
    savename = strcat(savedir,'bk_ps/',dt.name);
    for iter=(0:4)*10 + run
        [~,~,~]=...
        stackihl_ps_randomN(N_arr,dx,cbmap,psmap,...
        mask_inst_clip.*strmask,savename,iter);
    end
end

return