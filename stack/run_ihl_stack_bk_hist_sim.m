function run_ihl_stack_bk_hist_sim(flight,inst,field,hsc_idx,...
    run,iset,f_ihl,rvir,spire)

ifield=8;
mypaths=get_paths(flight);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');

if strcmp(field,'SIDES')
    if rvir==1
        load(sprintf('%sstackmapdatsim%d',loaddir,f_ihl*100),'stackmapdatsim');
    else
        load(sprintf('%sstackmapdatsim%d_rv%d',loaddir,f_ihl*100,rvir),...
            'stackmapdatsim');
    end
elseif strcmp(field,'HSC')
    name = HSC_field_name(hsc_idx);
    if rvir==1
        load(sprintf('%sstackmapdathsc_%s%d',loaddir,name,f_ihl*100),...
            'stackmapdatsim');
    else
        load(sprintf('%sstackmapdathsc_%s%d_rv%d',loaddir,name,f_ihl*100,rvir),...
            'stackmapdatsim');
    end
end


savedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

if iset==0
    cbmap = stackmapdatsim(ifield).all.cbmap;
    psmap = stackmapdatsim(ifield).all.psmap;
    mask_inst = stackmapdatsim(ifield).all.mask_inst_clip;
    strmask = stackmapdatsim(ifield).all.strmask;
elseif iset==1
    cbmap = stackmapdatsim(ifield).sub.cbmap;
    psmap = stackmapdatsim(ifield).sub.psmap;
    mask_inst = stackmapdatsim(ifield).sub.mask_inst_clip;
    strmask = stackmapdatsim(ifield).sub.strmask;
end

dx = 1200;
dt=get_dark_times(flight,inst,ifield);

if strcmp(field,'SIDES')
    if rvir==1
    load(sprintf('%s/histdatsim_%d',loaddir,f_ihl*100),'histdat');
    else
    load(sprintf('%s/histdatsim_%d_rv%d',loaddir,f_ihl*100,rvir),...
        'histdat');
    end
elseif strcmp(field,'HSC')
    if rvir==1
    load(sprintf('%s/histdathsc_%s%d',loaddir,name,f_ihl*100),'histdat');
    else
    load(sprintf('%s/histdathsc_%s%d_rv%d',loaddir,name,f_ihl*100,rvir),...
        'histdat');
    end
end

counts_arr = [];
countg_arr = [];
for im = 1:3
    counts_arr = [counts_arr, histdat(im).counts];
    countg_arr = [countg_arr, histdat(im).countg];
end

N_arr = unique([counts_arr, countg_arr]);
N_arr = N_arr(find(N_arr>0));
Ibinedges_cb = histdat(im).Ibinedges_cb;
Ibinedges_ps = histdat(im).Ibinedges_ps;

if strcmp(field,'SIDES')
    Npix = 720;
elseif strcmp(field,'HSC')
    Npix = 642;
end

for iter=(0:4)*10 + run
    [histbkdat]=stackihl_ps_randomN_hist(N_arr,dx,cbmap,psmap,...
        Ibinedges_cb,Ibinedges_ps,mask_inst.*strmask,Npix,spire);
    
    if strcmp(field,'SIDES')
        if rvir==1
        save(strcat(savedir,'bk_ps/','histbkdatsim',...
            num2str(f_ihl*100),'_',num2str(iter)),'histbkdat');    
        else
        save(strcat(savedir,'bk_ps/','histbkdatsim',...
            num2str(f_ihl*100),'_rv',num2str(rvir),'_',num2str(iter)),'histbkdat');
        end
    elseif strcmp(field,'HSC')
        if rvir==1
        save(strcat(savedir,'bk_ps/','histbkdathsc_',name,...
            num2str(f_ihl*100),'_',num2str(iter)),'histbkdat');    
        else
        save(strcat(savedir,'bk_ps/','histbkdathsc_',name,...
            num2str(f_ihl*100),'_rv',num2str(rvir),'_',num2str(iter)),'histbkdat');
        end
    end
    
end
return