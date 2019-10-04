function run_ihl_stack_map(flight,inst,ifield,varargin)
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

% load(sprintf('%s/TM%d/psfdat',mypaths.alldat,inst),'psfdatallfields');
load(sprintf('%s/TM%d/psfdat_%s',mypaths.alldat,inst,dt.namd),'psfdatall');
dx = 1200;
verbose = false;
cbmap = stackmapdat(ifield).cbmap;
psmap = stackmapdat(ifield).psmap;
m_min_arr = 16:19;
m_max_arr = 17:20;
Nbg = 100;
Njack = 100;

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');

%%%%%%%%%%%% just for testing !!!!!!!!!!
% if masklim
%     load(sprintf('%s/stackdat_%s_masklim',...
%         savedir,dt.name),'stackdatall');        
% else
%     load(sprintf('%s/stackdat_%s',...
%         savedir,dt.name),'stackdatall');
% end
% stackdatallold = stackdatall;

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
%     r_arr = stackdatallold(im).stackdat.r_arr;%%%%%%%
    stackdat.r_arr = r_arr;
    sp100 = find(r_arr>100);
    mask_inst = squeeze(mask_inst(inst,:,:));

%     psf = psfdatallfields(ifield).psf(im).psf;
%     psfps = psfdatallfields(ifield).psf(im).psfps;
%     psf_var = psfdatallfields(ifield).psf(im).psf_err.^2;
%     psfps_var = psfdatallfields(ifield).psf(im).psfps_err.^2;
    %%
    for isub=1:Njack
        [~,~,~,profcbg,profpsg,profhitg] = ...
            stackihl_ps0_hist_map(flight,inst,ifield,dx,cbmap,psmap,mask_inst,...
            strmask,strnum,1,verbose,nan,clipmaxs,clipmins,...
            srcdat.sub(isub).xg_arr,srcdat.sub(isub).yg_arr,...
            srcdat.sub(isub).mg_arr,subpix);
%         profcbg = stackdatallold(im).stackdat.sub(isub).profcbg;%%%%%%
%         profpsg = stackdatallold(im).stackdat.sub(isub).profpsg;%%%%%%
%         profhitg = stackdatallold(im).stackdat.sub(isub).profhitg;%%%%%%
        
        fprintf('stack %s, %d<m<%d, %d srcs, isub %d\n',...
            dt.name,m_min,m_max, srcdat.sub(isub).Ng,isub);

        stackdat.sub(isub).counts = srcdat.sub(isub).Ns;
        stackdat.sub(isub).countg = srcdat.sub(isub).Ng;
        profcbg(profhitg==0) = 0;
        profpsg(profhitg==0) = 0;

        stackdat.sub(isub).profcbg = profcbg;
        stackdat.sub(isub).profpsg = profpsg;
        stackdat.sub(isub).profhitg = profhitg;
    end
    %% profile combining all subset
    profcbg = zeros(size(stackdat.r_arr));
    profpsg = zeros(size(stackdat.r_arr));
    profhitg = zeros(size(stackdat.r_arr));
    counts = 0;
    countg = 0;
    for isub=1:Njack
        profcbg = profcbg + stackdat.sub(isub).profcbg.*stackdat.sub(isub).profhitg;
        profpsg = profpsg + stackdat.sub(isub).profpsg.*stackdat.sub(isub).profhitg;
        profhitg = profhitg + stackdat.sub(isub).profhitg;
        counts = counts + stackdat.sub(isub).counts;
        countg = countg + stackdat.sub(isub).countg;
    end
    stackdat.all.profcbg = profcbg./profhitg;
    stackdat.all.profpsg = profpsg./profhitg;
    stackdat.all.profhitg = profhitg;
    stackdat.all.counts = counts;
    stackdat.all.countg = countg;
    
    [prof15,prof100] = profile_radial_binning...
        (stackdat.r_arr,stackdat.all.profhitg,sp100);
    stackdat.rsub_arr = prof15;
    stackdat.r100 = prof100;
    
    [prof15,prof100] = profile_radial_binning...
        (stackdat.all.profcbg,stackdat.all.profhitg,sp100);
    stackdat.all.profcbgsub = prof15;
    stackdat.all.profcbg100 = prof100;
    
    [prof15,prof100] = profile_radial_binning...
        (stackdat.all.profpsg,stackdat.all.profhitg,sp100);
    stackdat.all.profpsgsub = prof15;
    stackdat.all.profpsg100 = prof100;
    
    %% profile of jackknife samples (leave one out)
    for isub=1:Njack
        jackcbg = stackdat.all.profcbg.*stackdat.all.profhitg...
            - stackdat.sub(isub).profcbg.*stackdat.sub(isub).profhitg;
        jackpsg = stackdat.all.profpsg.*stackdat.all.profhitg...
            - stackdat.sub(isub).profpsg.*stackdat.sub(isub).profhitg;
        jackhitg = stackdat.all.profhitg - stackdat.sub(isub).profhitg;
        stackdat.jack(isub).profcbg = jackcbg./jackhitg;
        stackdat.jack(isub).profpsg = jackpsg./jackhitg; 
        stackdat.jack(isub).profhitg = jackhitg; 

        [prof15,prof100] = profile_radial_binning...
            (stackdat.jack(isub).profcbg,stackdat.jack(isub).profhitg,sp100);
        stackdat.jack(isub).profcbgsub = prof15;
        stackdat.jack(isub).profcbg100 = prof100;
        
        [prof15,prof100] = profile_radial_binning...
            (stackdat.jack(isub).profpsg,stackdat.jack(isub).profhitg,sp100);
        stackdat.jack(isub).profpsgsub = prof15;
        stackdat.jack(isub).profpsg100 = prof100;        
    end
    %% cov from jackknife
    dat_profcbg = zeros([Njack,numel(stackdat.r_arr)]);
    dat_profpsg = zeros([Njack,numel(stackdat.r_arr)]);
    dat_profcbgsub = zeros([Njack,numel(stackdat.rsub_arr)]);
    dat_profpsgsub = zeros([Njack,numel(stackdat.rsub_arr)]);
    dat_profcbg100 = zeros([1,Njack]);
    dat_profpsg100 = zeros([1,Njack]);
    
    for isub=1:Njack
        dat_profcbg(isub,:) = stackdat.jack(isub).profcbg;
        dat_profpsg(isub,:) = stackdat.jack(isub).profpsg;
        dat_profcbgsub(isub,:) = stackdat.jack(isub).profcbgsub;
        dat_profpsgsub(isub,:) = stackdat.jack(isub).profpsgsub;
        dat_profcbg100(isub) = stackdat.jack(isub).profcbg100;
        dat_profpsg100(isub) = stackdat.jack(isub).profpsg100;
    end
    
    covcbg = zeros(numel(stackdat.r_arr));
    covpsg = zeros(numel(stackdat.r_arr));
    for isub=1:numel(stackdat.r_arr)
        for jsub=1:numel(stackdat.r_arr)
            dati = dat_profcbg(:,isub);
            datj = dat_profcbg(:,jsub);
            covcbg(isub,jsub) = mean(dati.*datj) - mean(dati)*mean(datj);
            dati = dat_profpsg(:,isub);
            datj = dat_profpsg(:,jsub);
            covpsg(isub,jsub) = mean(dati.*datj) - mean(dati)*mean(datj);
        end
    end
    
    covcbgsub = zeros(numel(stackdat.rsub_arr));
    covpsgsub = zeros(numel(stackdat.rsub_arr));    
    for isub=1:numel(stackdat.rsub_arr)
        for jsub=1:numel(stackdat.rsub_arr)
            dati = dat_profcbgsub(:,isub);
            datj = dat_profcbgsub(:,jsub);
            covcbgsub(isub,jsub) = mean(dati.*datj) - mean(dati)*mean(datj);
            dati = dat_profpsg(:,isub);
            datj = dat_profpsg(:,jsub);
            covpsgsub(isub,jsub) = mean(dati.*datj) - mean(dati)*mean(datj);
        end
    end

    errcbg100 = mean(dat_profcbg100.^2) - mean(dat_profcbg100)^2;
    errpsg100 = mean(dat_profpsg100.^2) - mean(dat_profpsg100)^2;
    
    stackdat.datcov.profcbg = covcbg.*(Njack-1);
    stackdat.datcov.profpsg = covpsg.*(Njack-1);
    stackdat.datcov.profcbg_rho = normalize_cov(stackdat.datcov.profcbg);
    stackdat.datcov.profpsg_rho = normalize_cov(stackdat.datcov.profpsg);
    stackdat.datcov.profcbgsub = covcbgsub.*(Njack-1);
    stackdat.datcov.profpsgsub = covpsgsub.*(Njack-1);
    stackdat.datcov.profcbgsub_rho = normalize_cov(stackdat.datcov.profcbgsub);
    stackdat.datcov.profpsgsub_rho = normalize_cov(stackdat.datcov.profpsgsub);
    stackdat.datcov.profcbg100 = errcbg100.*(Njack-1);
    stackdat.datcov.profpsg100 = errpsg100.*(Njack-1);
    
    %% background stack 
    r_arr = stackdat.r_arr;
    rsub_arr = stackdat.rsub_arr;
    profcbg_arr = zeros([Nbg,numel(r_arr)]);
    profpsg_arr = zeros([Nbg,numel(r_arr)]);
    profcbgsub_arr = zeros([Nbg,numel(rsub_arr)]);
    profpsgsub_arr = zeros([Nbg,numel(rsub_arr)]);
    profcbg100_arr = zeros([1,Nbg]);
    profpsg100_arr = zeros([1,Nbg]);
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
    
    %%% BG cov %%%%%
    covcbg = zeros(numel(stackdat.r_arr));
    covpsg = zeros(numel(stackdat.r_arr));
    for isub=1:numel(stackdat.r_arr)
        for jsub=1:numel(stackdat.r_arr)
            dati = profcbg_arr(:,isub);
            datj = profcbg_arr(:,jsub);
            covcbg(isub,jsub) = mean(dati.*datj) - mean(dati)*mean(datj);
            dati = profpsg_arr(:,isub);
            datj = profpsg_arr(:,jsub);
            covpsg(isub,jsub) = mean(dati.*datj) - mean(dati)*mean(datj);
        end
    end

    covcbgsub = zeros(numel(stackdat.rsub_arr));
    covpsgsub = zeros(numel(stackdat.rsub_arr));    
    for isub=1:numel(stackdat.rsub_arr)
        for jsub=1:numel(stackdat.rsub_arr)
            dati = profcbgsub_arr(:,isub);
            datj = profcbgsub_arr(:,jsub);
            covcbgsub(isub,jsub) = mean(dati.*datj) - mean(dati)*mean(datj);
            dati = profpsgsub_arr(:,isub);
            datj = profpsgsub_arr(:,jsub);
            covpsgsub(isub,jsub) = mean(dati.*datj) - mean(dati)*mean(datj);
        end
    end
    
    errcbg100 = mean(profcbg100_arr.^2) - mean(profcbg100_arr)^2;
    errpsg100 = mean(profpsg100_arr.^2) - mean(profpsg100_arr)^2;

    stackdat.bgcov.profcbg = covcbg;
    stackdat.bgcov.profpsg = covpsg;
    stackdat.bgcov.profcbg_rho = normalize_cov(covcbg);
    stackdat.bgcov.profpsg_rho = normalize_cov(covpsg);
    stackdat.bgcov.profcbgsub = covcbgsub;
    stackdat.bgcov.profpsgsub = covpsgsub;
    stackdat.bgcov.profcbgsub_rho = normalize_cov(covcbgsub);
    stackdat.bgcov.profpsgsub_rho = normalize_cov(covpsgsub);
    stackdat.bgcov.profcbg100 = errcbg100;
    stackdat.bgcov.profpsg100 = errpsg100;
    %% BG-sub profile
    stackdat.bgsub.profcbg = stackdat.all.profcbg - stackdat.bg.profcbg;
    stackdat.bgsub.profpsg = stackdat.all.profpsg - stackdat.bg.profpsg;
    stackdat.bgsub.profcbgsub = stackdat.all.profcbgsub - stackdat.bg.profcbgsub;
    stackdat.bgsub.profpsgsub = stackdat.all.profpsgsub - stackdat.bg.profpsgsub;
    stackdat.bgsub.profcbg100 = stackdat.all.profcbg100 - stackdat.bg.profcbg100;
    stackdat.bgsub.profpsg100 = stackdat.all.profpsg100 - stackdat.bg.profpsg100;
    
    psf = psfdatall.comb(im).all.profcb;
    psfps = psfdatall.comb(im).all.profps;

    stackdat.bgsub.profcbpsf = psf.*stackdat.bgsub.profcbg(1);
    stackdat.bgsub.profpspsf = psfps.*stackdat.bgsub.profpsg(1);
    
    [psf15,psf100] = profile_radial_binning(psf,stackdat.all.profhitg,sp100);
    [psfps15,psfps100] = profile_radial_binning(psfps,stackdat.all.profhitg,sp100);
    stackdat.bgsub.profcbpsfsub = psf15.*stackdat.bgsub.profcbg(1);
    stackdat.bgsub.profpspsfsub = psfps15.*stackdat.bgsub.profpsg(1);
    stackdat.bgsub.profcbpsf100 = psf100.*stackdat.bgsub.profcbg(1);
    stackdat.bgsub.profpspsf100 = psfps100.*stackdat.bgsub.profpsg(1);
    
    %% psf cov (var)
%     stackdat.psfcov.profcbpsf = diag(psf_var.*stackdat.bgsub.profcbg(1));
%     stackdat.psfcov.profpspsf = diag(psfps_var.*stackdat.bgsub.profpsg(1));
%     [psf_var15,psf_var100] = profile_radial_binning(psf_var,...
%         stackdat.all.profhitg,sp100);
%     [psfps_var15,psfps_var100] = profile_radial_binning(psfps_var,...
%         stackdat.all.profhitg,sp100);
%     stackdat.psfcov.profcbpsfsub = diag(psf_var15.*stackdat.bgsub.profcbg(1));
%     stackdat.psfcov.profpspsfsub = diag(psfps_var15.*stackdat.bgsub.profpsg(1));
%     stackdat.psfcov.profcbpsf100 = diag(psf_var100.*stackdat.bgsub.profcbg(1));
%     stackdat.psfcov.profpspsf100 = diag(psfps_var100.*stackdat.bgsub.profpsg(1));
    scalecb = stackdat.bgsub.profcbg(1);
    scaleps = stackdat.bgsub.profpsg(1);
    stackdat.psfcov.profcbpsf = psfdatall.comb(im).datcov.profcb.*scalecb.^2;
    stackdat.psfcov.profpspsf = psfdatall.comb(im).datcov.profps.*scaleps.^2;
    stackdat.psfcov.profcbpsfsub = psfdatall.comb(im).datcov.profcbsub.*scalecb.^2;
    stackdat.psfcov.profpspsfsub = psfdatall.comb(im).datcov.profpssub.*scaleps.^2;
    stackdat.psfcov.profcbpsf100 = psfdatall.comb(im).datcov.profcb100.*scalecb.^2;
    stackdat.psfcov.profpspsf100 = psfdatall.comb(im).datcov.profps100.*scaleps.^2;    
     %% tot cov
    stackdat.cov.cb = stackdat.datcov.profcbg + ...
        stackdat.bgcov.profcbg + stackdat.psfcov.profcbpsf;
    
    stackdat.cov.ps = stackdat.datcov.profpsg + ...
        stackdat.bgcov.profpsg + stackdat.psfcov.profpspsf;
    
    stackdat.cov.cbsub = stackdat.datcov.profcbgsub + ...
        stackdat.bgcov.profcbgsub + stackdat.psfcov.profcbpsfsub;
    
    stackdat.cov.pssub = stackdat.datcov.profpsgsub + ...
        stackdat.bgcov.profpsgsub + stackdat.psfcov.profpspsfsub;

    stackdat.cov.cb100 = stackdat.datcov.profcbg100 + ...
        stackdat.bgcov.profcbg100 + stackdat.psfcov.profcbpsf100;
    
    stackdat.cov.ps100 = stackdat.datcov.profpsg100 + ...
        stackdat.bgcov.profpsg100 + stackdat.psfcov.profpspsf100;
    
    stackdat.cov.cb_rho = normalize_cov(stackdat.cov.cb);
    stackdat.cov.ps_rho = normalize_cov(stackdat.cov.ps);    
    stackdat.cov.cbsub_rho = normalize_cov(stackdat.cov.cbsub);
    stackdat.cov.pssub_rho = normalize_cov(stackdat.cov.pssub);
    stackdat.cov.cbsub_inv = inv(stackdat.cov.cbsub); 
    stackdat.cov.pssub_inv = inv(stackdat.cov.pssub); 
    
    %% excess profile
    stackdat.excess.cb = stackdat.bgsub.profcbg - stackdat.bgsub.profcbpsf;
    stackdat.excess.ps = stackdat.bgsub.profpsg - stackdat.bgsub.profpspsf;
    stackdat.excess.cbsub = stackdat.bgsub.profcbgsub - stackdat.bgsub.profcbpsfsub;
    stackdat.excess.pssub = stackdat.bgsub.profpsgsub - stackdat.bgsub.profpspsfsub;
    stackdat.excess.cb100 = stackdat.bgsub.profcbg100 - stackdat.bgsub.profcbpsf100;
    stackdat.excess.ps100 = stackdat.bgsub.profpsg100 - stackdat.bgsub.profpspsf100;

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
