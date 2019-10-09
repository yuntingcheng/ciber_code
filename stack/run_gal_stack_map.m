function run_gal_stack_map(flight,inst,ifield,varargin)
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('ifield',@isnumeric);
  p.addOptional('masklim',false,@islogical);
  p.addOptional('sample_type','jack_random',@ischar);
  p.addOptional('subpix',true,@isnumeric);
  
  p.parse(flight,inst,ifield,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  ifield   = p.Results.ifield;
  masklim = p.Results.masklim;
  sample_type=p.Results.sample_type;
  subpix   = p.Results.subpix;
  clear p varargin;
%%
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);
load(sprintf('%s/TM%d/stackmapdat',mypaths.alldat,1),'stackmapdat');
stackmapdat1 = stackmapdat;
load(sprintf('%s/TM%d/stackmapdat',mypaths.alldat,2),'stackmapdat');
stackmapdat2 = stackmapdat;

if inst==1
    stackmapdat = stackmapdat1;
else
    stackmapdat = stackmapdat2;
end

dx = 1200;
verbose = false;
cbmap = stackmapdat(ifield).cbmap;
psmap = stackmapdat(ifield).psmap;
m_min_arr = 16:19;
m_max_arr = 17:20;
Njack = 50;
Nbg = 100;

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
%%
for im= 1:numel(m_min_arr)
%%
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    mask_inst = zeros([2,1024,1024]);
    if masklim
        mask_inst(1,:,:) = stackmapdat1(ifield).m_max(m_max).mask_inst_clip;
        mask_inst(2,:,:) = stackmapdat2(ifield).m_max(m_max).mask_inst_clip;
        
        strmask = stackmapdat(ifield).m_max(m_max).strmask;
        strnum = stackmapdat(ifield).m_max(m_max).strnum;
    else
        mask_inst(1,:,:) = stackmapdat1(ifield).mask_inst_clip;
        mask_inst(2,:,:) = stackmapdat2(ifield).mask_inst_clip;
        strmask = stackmapdat(ifield).strmask;
        strnum = stackmapdat(ifield).strnum;    
    end
    
    stackdat.m_min = m_min;
    stackdat.m_max = m_max;
    
    srcdat = ps_src_select(flight,inst,ifield,m_min,m_max,mask_inst,...
    'Nsub', Njack, 'sample_type',sample_type);
  
    [clipmaxs, clipmins, r_arr]=...
    stackihl_ps0_cliplim(flight,inst,ifield,m_min,m_max,cbmap,psmap,...
    mask_inst,strnum,1000,verbose,[],nan,false);
    stackdat.r_arr = r_arr;
    sp100 = find(r_arr>100);
    mask_inst = squeeze(mask_inst(inst,:,:));
    %%
    for isub=1:Njack
        [~,~,~,profcbg,profpsg,profhitg] = ...
            stackihl_ps0_hist_map(flight,inst,ifield,dx,cbmap,psmap,mask_inst,...
            strmask,strnum,1,verbose,nan,clipmaxs,clipmins,...
            srcdat.sub(isub).xg_arr,srcdat.sub(isub).yg_arr,...
            srcdat.sub(isub).mg_arr,subpix);
        
        fprintf('stack %s, %d<m<%d, %d srcs, isub %d\n',...
            dt.name,m_min,m_max, srcdat.sub(isub).Ng,isub);

        stackdat.sub(isub).counts = srcdat.sub(isub).Ns;
        stackdat.sub(isub).countg = srcdat.sub(isub).Ng;
        profcbg(profhitg==0) = 0;
        profpsg(profhitg==0) = 0;

        stackdat.sub(isub).profcbg = profcbg;
        stackdat.sub(isub).profpsg = profpsg;
        stackdat.sub(isub).profhitg = profhitg;
        
        [profcb,profps,~]=stackihl_ps0_hist_map_bk...
            (dx,cbmap,psmap,mask_inst,strmask,...
            [srcdat.sub(isub).Ng-1,srcdat.sub(isub).Ng],verbose,false);
        
        % interpolate missing data
        profcb = profcb(2,:);
        sp = find(profcb==profcb);
        profcb = spline(r_arr(sp),profcb(sp),r_arr);
        
        profps = profps(2,:);
        sp = find(profps==profps);
        profps = spline(r_arr(sp),profps(sp),r_arr);

        stackdat.bgsub(isub).profcbg = profcb;
        stackdat.bgsub(isub).profpsg = profps;

        fprintf('stack %s BG, %d<m<%d, %d gals, isub %d\n',...
            dt.name,m_min,m_max,srcdat.sub(isub).Ng,isub);
    end
    %% profile combining all subset
    profcbg = zeros(size(r_arr));
    profpsg = zeros(size(r_arr));
    profhitg = zeros(size(r_arr));
    profcbgbg = zeros(size(r_arr));
    profpsgbg = zeros(size(r_arr));
    counts = 0;
    countg = 0;
    for isub=1:Njack
        profcbg = profcbg + stackdat.sub(isub).profcbg.*stackdat.sub(isub).profhitg;
        profpsg = profpsg + stackdat.sub(isub).profpsg.*stackdat.sub(isub).profhitg;
        profcbgbg = profcbgbg + stackdat.bgsub(isub).profcbg;
        profpsgbg = profpsgbg + stackdat.bgsub(isub).profpsg;
        profhitg = profhitg + stackdat.sub(isub).profhitg;
        counts = counts + stackdat.sub(isub).counts;
        countg = countg + stackdat.sub(isub).countg;
    end
    stackdat.all.profcbg = profcbg./profhitg;
    stackdat.all.profpsg = profpsg./profhitg;
    stackdat.all.profcbgbg = profcbgbg./Njack;
    stackdat.all.profpsgbg = profpsgbg./Njack;
    stackdat.all.profhitg = profhitg;
    stackdat.all.counts = counts;
    stackdat.all.countg = countg;
    
    [rsub_arr,r100] = profile_radial_binning...
        (r_arr,stackdat.all.profhitg,sp100);
    stackdat.rsub_arr = rsub_arr;
    stackdat.r100 = r100;
    
    [prof15,prof100] = profile_radial_binning...
        (stackdat.all.profcbg,stackdat.all.profhitg,sp100);
    stackdat.all.profcbgsub = prof15;
    stackdat.all.profcbg100 = prof100;
    
    [prof15,prof100] = profile_radial_binning...
        (stackdat.all.profpsg,stackdat.all.profhitg,sp100);
    stackdat.all.profpsgsub = prof15;
    stackdat.all.profpsg100 = prof100;

    [prof15,prof100] = profile_radial_binning...
        (stackdat.all.profcbgbg,stackdat.all.profhitg,sp100);
    stackdat.all.profcbgbgsub = prof15;
    stackdat.all.profcbgbg100 = prof100;
    
    [prof15,prof100] = profile_radial_binning...
        (stackdat.all.profpsgbg,stackdat.all.profhitg,sp100);
    stackdat.all.profpsgbgsub = prof15;
    stackdat.all.profpsgbg100 = prof100;
    %% profile of jackknife samples (leave one out)
    for isub=1:Njack
        jackcbg = stackdat.all.profcbg.*stackdat.all.profhitg...
            - stackdat.sub(isub).profcbg.*stackdat.sub(isub).profhitg;
        jackpsg = stackdat.all.profpsg.*stackdat.all.profhitg...
            - stackdat.sub(isub).profpsg.*stackdat.sub(isub).profhitg;
        jackhitg = stackdat.all.profhitg - stackdat.sub(isub).profhitg;
        jackcbgbg = Njack.*stackdat.all.profcbgbg - stackdat.bgsub(isub).profcbg;
        jackpsgbg = Njack.*stackdat.all.profpsgbg - stackdat.bgsub(isub).profpsg;
        stackdat.jack(isub).profcbg = jackcbg./jackhitg;
        stackdat.jack(isub).profpsg = jackpsg./jackhitg; 
        stackdat.jack(isub).profhitg = jackhitg; 
        stackdat.bgjack(isub).profcbg = jackcbgbg./(Njack-1);
        stackdat.bgjack(isub).profpsg = jackpsgbg./(Njack-1); 

        [prof15,prof100] = profile_radial_binning...
            (stackdat.jack(isub).profcbg,stackdat.jack(isub).profhitg,sp100);
        stackdat.jack(isub).profcbgsub = prof15;
        stackdat.jack(isub).profcbg100 = prof100;
        
        [prof15,prof100] = profile_radial_binning...
            (stackdat.jack(isub).profpsg,stackdat.jack(isub).profhitg,sp100);
        stackdat.jack(isub).profpsgsub = prof15;
        stackdat.jack(isub).profpsg100 = prof100;        

        [prof15,prof100] = profile_radial_binning...
            (stackdat.bgjack(isub).profcbg,stackdat.jack(isub).profhitg,sp100);
        stackdat.bgjack(isub).profcbgsub = prof15;
        stackdat.bgjack(isub).profcbg100 = prof100;
        
        [prof15,prof100] = profile_radial_binning...
            (stackdat.bgjack(isub).profpsg,stackdat.jack(isub).profhitg,sp100);
        stackdat.bgjack(isub).profpsgsub = prof15;
        stackdat.bgjack(isub).profpsg100 = prof100;        
        
    end   
    %% stack Nbg of BG with Ng random positions
    r_arr = stackdat.r_arr;
    rsub_arr = stackdat.rsub_arr;
    profcbg_arr = zeros([Nbg,numel(r_arr)]);
    profpsg_arr = zeros([Nbg,numel(r_arr)]);
    profcbgsub_arr = zeros([Nbg,numel(rsub_arr)]);
    profpsgsub_arr = zeros([Nbg,numel(rsub_arr)]);
    profcbg100_arr = zeros([Nbg,1]);
    profpsg100_arr = zeros([Nbg,1]);
    for isim=1:Nbg
        [profcb,profps,hitmap]=stackihl_ps0_hist_map_bk...
            (dx,cbmap,psmap,mask_inst,strmask,[srcdat.Ng-1,srcdat.Ng],...
            verbose,false);
        
        % interpolate missing data
        profcb = profcb(2,:);
        sp = find(profcb==profcb);
        profcb = spline(r_arr(sp),profcb(sp),r_arr);
        profcbg_arr(isim,:) = profcb;
        
        profps = profps(2,:);
        sp = find(profps==profps);
        profps = spline(r_arr(sp),profps(sp),r_arr);
        profpsg_arr(isim,:) = profps;
        
        [prof15,prof100] = profile_radial_binning(profcb,hitmap,sp100);
        profcbgsub_arr(isim,:) = prof15;
        profcbg100_arr(isim) = prof100;

        [prof15,prof100] = profile_radial_binning(profps,hitmap,sp100);
        profpsgsub_arr(isim,:) = prof15;
        profpsg100_arr(isim) = prof100;
        
        fprintf('stack %s, %d<m<%d, %d gals, isim %d\n',...
            dt.name,m_min,m_max,srcdat.Ng,isim);
    end
    stackdat.bg.profcbg = mean(profcbg_arr);
    stackdat.bg.profpsg = mean(profpsg_arr);
    stackdat.bg.profcbgsub = mean(profcbgsub_arr);
    stackdat.bg.profpsgsub = mean(profpsgsub_arr);
    stackdat.bg.profcbg100 = mean(profcbg100_arr);
    stackdat.bg.profpsg100 = mean(profpsg100_arr);
    
    stackdat.bgcov.covcb = get_cov_matrix(profcbg_arr);
    stackdat.bgcov.covps = get_cov_matrix(profpsg_arr);
    stackdat.bgcov.covcbsub = get_cov_matrix(profcbgsub_arr);
    stackdat.bgcov.covpssub = get_cov_matrix(profpsgsub_arr);
    stackdat.bgcov.covcb100 = get_cov_matrix(profcbg100_arr);
    stackdat.bgcov.covps100 = get_cov_matrix(profpsg100_arr);
    %% save data  
    stackdatall(im).stackdat = stackdat;
    if masklim
        save(sprintf('%s/stackdat_%s_masklim',...
            savedir,dt.name),'stackdatall');        
    else
        save(sprintf('%s/stackdat_%s',...
            savedir,dt.name),'stackdatall');
    end
    clear stackdat
end

return
