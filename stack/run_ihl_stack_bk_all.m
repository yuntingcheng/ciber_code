function run_ihl_stack_bk_all(flight,inst,ifield_arr,run)
mypaths=get_paths(flight);
savedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

dx = 1200;
for ifield=ifield_arr
dt=get_dark_times(flight,inst,ifield);
cbmap = stackmapdat(ifield).cbmap;
psmap = stackmapdat(ifield).psmap;
mask_inst_clip = stackmapdat(ifield).mask_inst_clip;
strmask = stackmapdat(ifield).strmask;
counts_arr = stackmapdat(ifield).count_stacks_arr;
countg_arr = stackmapdat(ifield).count_stackg_arr;

N_arr = [sum(countg_arr(10:12)),sum(counts_arr(10:12))];


savename = strcat(savedir,'bk_ps/',dt.name);
for iter=(1:10)+(run-1)*10
     [stampercb,stamperps,maskstamper]=...
     stackihl_ps_randomN(N_arr,dx,cbmap,psmap,...
     mask_inst_clip.*strmask,savename,iter);
end

end


return