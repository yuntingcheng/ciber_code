function run_ihl_stack_hist_subreg(flight,inst,ifield,spire)
mypaths=get_paths(flight);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

dx = 1200;
verbose = 1;
dt=get_dark_times(flight,inst,ifield);
cbmap = stackmapdat(ifield).cbmap;
psmap = stackmapdat(ifield).psmap;
mask_inst = stackmapdat(ifield).mask_inst_clip;
strmask = stackmapdat(ifield).strmask;
strnum = stackmapdat(ifield).strnum;
m_min_arr = stackmapdat(ifield).m_min_arr;
m_max_arr = stackmapdat(ifield).m_max_arr;

if inst == 1
    stackband = 'I';
else
    stackband = 'H';
end

for im= 10:12
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    
    type = -1;
    [histcb, histps,Ibinedges_cb,Ibinedges_ps,stackcount]=...
    stackihl_ps0_hist_subreg(flight,inst,ifield,type,m_min,m_max,dx,cbmap,psmap,...
    mask_inst,strmask,strnum,1,0,verbose,10,stackband,[],spire);
    histdatall(im).histcbs = histcb;
    histdatall(im).histpss = histps;
    histdatall(im).Ibinedges_cb = Ibinedges_cb;
    histdatall(im).Ibinedges_ps = Ibinedges_ps;
    histdatall(im).counts_arr = stackcount;
    
    type = 1;
    [histcb, histps,Ibinedges_cb,Ibinedges_ps,stackcount]=...
    stackihl_ps0_hist_subreg(flight,inst,ifield,type,m_min,m_max,dx,cbmap,psmap,...
    mask_inst,strmask,strnum,1,0,verbose,10,stackband,[],spire);
    histdatall(im).histcbg = histcb;
    histdatall(im).histpsg = histps;
    histdatall(im).countg_arr = stackcount;

end

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
for isub=1:8
    for im=10:12
        histdat(im).Ibinedges_cb = histdatall(im).Ibinedges_cb;
        histdat(im).Ibinedges_ps = histdatall(im).Ibinedges_ps;
        histdat(im).counts=histdatall(im).counts_arr(isub);
        histdat(im).countg=histdatall(im).countg_arr(isub);
        dat = histdatall(im).histcbs(isub,:,:,:);s = size(dat);
        histdat(im).histcbs = reshape(dat,s(2:end));
        dat = histdatall(im).histcbg(isub,:,:,:);s = size(dat);
        histdat(im).histcbg = reshape(dat,s(2:end));
        dat = histdatall(im).histpss(isub,:,:,:);s = size(dat);
        histdat(im).histpss = reshape(dat,s(2:end));
        dat = histdatall(im).histpsg(isub,:,:,:);s = size(dat);
        histdat(im).histpsg = reshape(dat,s(2:end));
        save(sprintf('%s/histdat_%s_subreg%d',savedir,dt.name,isub),'histdat');
    end
end
return