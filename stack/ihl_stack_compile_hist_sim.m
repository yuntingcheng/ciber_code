function ihl_stack_compile_hist_sim(flight,inst,field,varargin)
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('field',@ischar);
  p.addOptional('hsc_idx',1,@isnumeric);
  p.addOptional('f_ihl',0,@isnumeric);
  p.addOptional('rvir',1,@isnumeric);
  p.addOptional('psf_model_stars',false,@islogical);
  
  p.parse(flight,inst,field,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  field    = p.Results.field;
  hsc_idx  = p.Results.hsc_idx;
  f_ihl    = p.Results.f_ihl;
  rvir     = p.Results.rvir;
  psf_model_stars     = p.Results.psf_model_stars;
  
  clear p varargin;
%%
ifield = 8;
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if strcmp(field,'SIDES')
    if rvir==1
    load(sprintf('%s/histdatsim_%d',loaddir,f_ihl*100),'histdat');
    else
    load(sprintf('%s/histdatsim_%d_rv%d',loaddir,f_ihl*100,rvir),...
        'histdat');
    end
elseif strcmp(field,'HSC')
    name = HSC_fields_info(hsc_idx);
    if rvir==1
    load(sprintf('%s/histdathsc_%s%d',loaddir,name,f_ihl*100),'histdat');
    else
    load(sprintf('%s/histdathsc_%s%d_rv%d',loaddir,name,f_ihl*100,rvir),...
        'histdat');
    end
end

psfdir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(psfdir,'fitpsfdat'),'fitpsfdat');

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

bkdir = (strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

nsim = 50;
dx = 1200;
nbins = 25;

plotbins = [];%[2,8,24];
m_min_arr = stackmapdat(ifield).m_min_arr;
m_max_arr = stackmapdat(ifield).m_max_arr;

profile = radial_prof(ones(2*dx+1),ones(2*dx+1),dx+1,dx+1,1,nbins);
rbinedges = profile.binedges * 0.7;
rbins = binedges2bins(rbinedges);

psfstack0 = fitpsfdat(ifield).psfmodel.prof_best_arr;
rpsfstack0 = fitpsfdat(ifield).psfmodel.r_arr;
psfstack = interp1(rpsfstack0,psfstack0,rbins,'linear','extrap');
psfstack(psfstack<0)=0;
psfstack = psfstack./psfstack(1);
ihlprofdat.r_arr = rbins;
ihlprofdat.binedges = rbinedges;
ihlprofdat.psf_arr = psfstack;

%%% get the PSF params %%%
% bestparam = fitpsfdat(ifield).bestparam;
% beta = bestparam(1);
% rc = bestparam(2);
% radmap = make_radius_map(zeros(2*dx+1),dx,dx).*0.7;
% psfmap = (1 + (radmap/rc).^2).^(-3.*beta./2);
% profile = radial_prof(psfmap,ones(2*dx+1),dx+1,dx+1,1,nbins);
% r_arr=profile.r*0.7;
% profpsf_arr=(profile.prof)./profile.prof(1);
% ihlprofdat.r_arr = r_arr;
% ihlprofdat.binedges = profile.binedges;
% ihlprofdat.psf_arr = profpsf_arr;
% %ihlprofdat.A = get_profile_mkk(flight,inst,ifield,dx,ihlprofdat.binedges);

sp100 = find(rbins>100);
%%
for im=1:3
    m_min = m_min_arr(im+9);
    m_max = m_max_arr(im+9);
    counts = histdat(im).counts;
    countg = histdat(im).countg;

    if counts<10 | countg<10
        continue
    end
    
    ihlprofdat.data(im).m_min = m_min;
    ihlprofdat.data(im).m_max = m_max;
    ihlprofdat.data(im).counts = counts;
    ihlprofdat.data(im).countg = countg;
    clipmax_arr=zeros([4,nbins]);
    clipmin_arr=zeros([4,nbins]);

    %%% find the index of S & G in bk stack file %%%
    if strcmp(field,'SIDES')
        if rvir==1
        load(strcat(bkdir,'bk_ps_bu/','histbkdatsim',...
            num2str(f_ihl*100),'_1'),'histbkdat');    
        else
        load(strcat(bkdir,'bk_ps_bu/','histbkdatsim',...
            num2str(f_ihl*100),'_rv',num2str(rvir),'_1'),'histbkdat');
        end
    elseif strcmp(field,'HSC')
        if rvir==1
        load(strcat(bkdir,'bk_ps_bu/','histbkdathsc_',name,...
            num2str(f_ihl*100),'_1'),'histbkdat');    
        else
        load(strcat(bkdir,'bk_ps_bu/','histbkdathsc_',name,...
            num2str(f_ihl*100),'_rv',num2str(rvir),'_1'),'histbkdat');
        end
    end
        
    for i=1:size(histbkdat,2)
        if histbkdat(i).N == counts
            bkidx_s = i;
        end
        if histbkdat(i).N == countg
            bkidx_g = i;
        end  
    end
    %%% get the stacking profile %%%
    profgcb = zeros([1,nbins]);
    errgcb = zeros([1,nbins]);
    profscb = zeros([1,nbins]);
    errscb = zeros([1,nbins]);
    profgps = zeros([1,nbins]);
    errgps = zeros([1,nbins]);
    profsps = zeros([1,nbins]);
    errsps = zeros([1,nbins]);

    profgcbbk = zeros([1,nbins]);
    errgcbbk = zeros([1,nbins]);
    profscbbk = zeros([1,nbins]);
    errscbbk = zeros([1,nbins]);
    profgpsbk = zeros([1,nbins]);
    errgpsbk = zeros([1,nbins]);
    profspsbk = zeros([1,nbins]);
    errspsbk = zeros([1,nbins]);
    
    for ibin=1:nbins
        fprintf('%s,im=%d, ibin=%d\n',dt.name,im,ibin);
        if ismember(ibin,plotbins)
            figure
            setwinsize(gcf,800,500);
        end

        for itype=1:4
            
            %%% get data %%%
            if itype==1
                bins = binedges2bins(histdat(im).Ibinedges_cb);
                d_all = histdat(im).histcbs(:,ibin,:);
                d_gs = histdat(im).histcbs(:,ibin,:) ...
                    +histdat(im).histcbg(:,ibin,:);
                tname = 'CIBER stars';
            elseif itype==2
                bins = binedges2bins(histdat(im).Ibinedges_cb);
                d_all = histdat(im).histcbg(:,ibin,:);
                d_gs = histdat(im).histcbs(:,ibin,:) ...
                    +histdat(im).histcbg(:,ibin,:);                
                tname = 'CIBER gals';
            elseif itype==3
                bins = binedges2bins(histdat(im).Ibinedges_ps);
                d_all = histdat(im).histpss(:,ibin,:);
                d_gs = histdat(im).histpss(:,ibin,:) ...
                    +histdat(im).histpsg(:,ibin,:);
                tname = 'Sim stars';
            elseif itype==4
                bins = binedges2bins(histdat(im).Ibinedges_ps);
                d_all = histdat(im).histpsg(:,ibin,:);
                d_gs = histdat(im).histpss(:,ibin,:) ...
                    +histdat(im).histpsg(:,ibin,:);             
                tname = 'Sim gals';
            end
            d = sum(d_all);
            d = reshape(d, size(bins));
            d_gs = sum(d_gs);
            d_gs = reshape(d_gs,size(bins));
            
            %%% get stacking profile %%%
            sp0 = find(d~=0);
            mean0 = sum(d.*bins)./sum(d);
            if rbinedges(ibin+1)<(get_mask_radius(2,8,m_min)/7+1)
                Q1 = quartile_from_hist(bins,d_gs,0.25);
                Q3 = quartile_from_hist(bins,d_gs,0.75);
                IQR = Q3 - Q1;
                clipmax = Q3+3*IQR;
                clipmin = Q1-3*IQR;
            else
                clipmax = inf;
                clipmin = -inf;
            end
            clipmax_arr(itype,ibin) = clipmax;
            clipmin_arr(itype,ibin) = clipmin;
            
            sp1 = find( (d~=0) & (bins < clipmax) & (bins > clipmin));
            mean1 = sum(d(sp1).*bins(sp1))./sum(d(sp1));
            
            %%% get stacking err from sub stack %%%
            Nsub = size(d_all,1);
            meansub_arr = zeros([1,Nsub]);
            for isub=1:Nsub
                d = d_all(isub,:,:);
                d = reshape(d, size(bins));
                spi = find((d~=0) & (bins < clipmax) & (bins > clipmin));
                meani = sum(d(spi).*bins(spi))./sum(d(spi));
                meansub_arr(isub) = meani;
            end
            Nsub = sum(~isnan(meansub_arr));
            
            %%% plot hist of stacking prof %%%
            if ismember(ibin,plotbins)
                subplot(2,2,itype)
                semilogy(bins(sp0),d(sp0),'b.');hold on
                semilogy(bins(sp1),d(sp1),'r.');
                %vline(mean0,'b-');
                vline(mean1,'r-');
                vline(clipmax,'r--');
                vline(clipmin,'r--');
                %vline(mean1 + std(meansub_arr)./sqrt(Nsub)./2 ,'k--');
                %vline(mean1 - std(meansub_arr)./sqrt(Nsub)./2 ,'k--');
                xlabel('I [nW/m^2/sr]');
                ylabel('counts');
                title(tname);
                xLimits = get(gca,'XLim');
                yLimits = get(gca,'YLim');
                if mean0 ~= mean1
                text(xLimits(1)+0.7*(xLimits(2) - xLimits(1)),...
                    0.5*yLimits(2),...
                    sprintf('<I>=%.2e',mean0),'color','b')    
                end
                text(xLimits(1)+0.7*(xLimits(2) - xLimits(1)),...
                    0.3*yLimits(2),...
                    sprintf('<I>=%.2e',mean1),'color','r')    
            end
            
            %%% get background prof %%%
            meanbk_arr = zeros([1,nsim]);
            
            for i=1:nsim
                if strcmp(field,'SIDES')
                    if rvir==1
                    load(strcat(bkdir,'bk_ps_bu/','histbkdatsim',...
                        num2str(f_ihl*100),'_',num2str(i)),'histbkdat');    
                    else
                    load(strcat(bkdir,'bk_ps_bu/','histbkdatsim',...
                        num2str(f_ihl*100),'_rv',num2str(rvir),...
                        '_',num2str(i)),'histbkdat');
                    end
                elseif strcmp(field,'HSC')
                    if rvir==1
                    load(strcat(bkdir,'bk_ps_bu/','histbkdathsc_',name,...
                        num2str(f_ihl*100),'_',num2str(i)),'histbkdat');    
                    else
                    load(strcat(bkdir,'bk_ps_bu/','histbkdathsc_',name,...
                        num2str(f_ihl*100),'_rv',num2str(rvir),...
                        '_',num2str(i)),'histbkdat');
                    end
                end

                if itype==1
                    d= histbkdat(bkidx_s).histcb(ibin,:);
                elseif itype==2
                    d= histbkdat(bkidx_g).histcb(ibin,:);
                elseif itype==3
                    d= histbkdat(bkidx_s).histps(ibin,:);
                elseif itype==4
                    d= histbkdat(bkidx_g).histps(ibin,:);
                end
                d = reshape(d, size(bins));
                spi = find((d~=0) & (bins < clipmax) & (bins > clipmin));
                if numel(spi)~=0
                    meani = sum(d(spi).*bins(spi))./sum(d(spi));
                else
                    meani = 0;
                end
                meanbk_arr(i) = meani;
            end
            
            %%% write stack profile data %%%
            if itype==1
                profscb(ibin) = mean1;
                errscb(ibin) = nanstd(meansub_arr)./sqrt(Nsub);
                profscbbk(ibin) = nanmean(meanbk_arr);
                errscbbk(ibin) = std(meanbk_arr);
            elseif itype==2
                profgcb(ibin) = mean1;
                errgcb(ibin) = nanstd(meansub_arr)./sqrt(Nsub);
                profgcbbk(ibin) = nanmean(meanbk_arr);
                errgcbbk(ibin) = std(meanbk_arr);
            elseif itype==3
                profsps(ibin) = mean1;
                errsps(ibin) = nanstd(meansub_arr)./sqrt(Nsub);
                profspsbk(ibin) = nanmean(meanbk_arr);
                errspsbk(ibin) = std(meanbk_arr);
            elseif itype==4
                profgps(ibin) = mean1;
                errgps(ibin) = nanstd(meansub_arr)./sqrt(Nsub);
                profgpsbk(ibin) = nanmean(meanbk_arr);
                errgpsbk(ibin) = std(meanbk_arr);
            end
            
                    
        end
        
        if ismember(ibin,plotbins)
            suptitle(sprintf('%dth radial bin, <r> = %.1e', ibin, rbins(ibin)));
        end
        
    end

    %%%%%%% average of r>100 arcsec %%%%%%%%%%%

    for itype=1:4

        %%% get data %%%
        if itype==1
            bins = binedges2bins(histdat(im).Ibinedges_cb);
            d_all = histdat(im).histcbs(:,sp100,:);
        elseif itype==2
            bins = binedges2bins(histdat(im).Ibinedges_cb);
            d_all = histdat(im).histcbg(:,sp100,:);
        elseif itype==3
            bins = binedges2bins(histdat(im).Ibinedges_ps);
            d_all = histdat(im).histpss(:,sp100,:);
        elseif itype==4
            bins = binedges2bins(histdat(im).Ibinedges_ps);
            d_all = histdat(im).histpsg(:,sp100,:);
        end
        d = sum(sum(d_all,2));
        d = reshape(d, size(bins));

        %%% get stacking profile %%%
        clipmax = inf;
        clipmin = -inf;
        sp1 = find( (d~=0) & (bins < clipmax) & (bins > clipmin));
        mean1 = sum(d(sp1).*bins(sp1))./sum(d(sp1));

        %%% get stacking err from sub stack %%%
        Nsub = size(d_all,1);
        meansub_arr = zeros([1,Nsub]);
        for isub=1:Nsub
            d = sum(d_all(isub,:,:),2);
            d = reshape(d, size(bins));
            spi = find((d~=0) & (bins < clipmax) & (bins > clipmin));
            meani = sum(d(spi).*bins(spi))./sum(d(spi));
            meansub_arr(isub) = meani;
        end
        Nsub = sum(~isnan(meansub_arr));

        %%% get background prof %%%
        meanbk_arr = zeros([1,nsim]);

        for i=1:nsim
            if strcmp(field,'SIDES')
                if rvir==1
                load(strcat(bkdir,'bk_ps_bu/','histbkdatsim',...
                    num2str(f_ihl*100),'_',num2str(i)),'histbkdat');    
                else
                load(strcat(bkdir,'bk_ps_bu/','histbkdatsim',...
                    num2str(f_ihl*100),'_rv',num2str(rvir),...
                    '_',num2str(i)),'histbkdat');
                end
            elseif strcmp(field,'HSC')
                if rvir==1
                load(strcat(bkdir,'bk_ps_bu/','histbkdathsc_',name,...
                    num2str(f_ihl*100),'_',num2str(i)),'histbkdat');    
                else
                load(strcat(bkdir,'bk_ps_bu/','histbkdathsc_',name,...
                    num2str(f_ihl*100),'_rv',num2str(rvir),...
                    '_',num2str(i)),'histbkdat');
                end
            end

            if itype==1
                d= histbkdat(bkidx_s).histcb(sp100,:);
            elseif itype==2
                d= histbkdat(bkidx_g).histcb(sp100,:);
            elseif itype==3
                d= histbkdat(bkidx_s).histps(sp100,:);
            elseif itype==4
                d= histbkdat(bkidx_g).histps(sp100,:);
            end
            d = reshape(sum(d), size(bins));
            spi = find((d~=0) & (bins < clipmax) & (bins > clipmin));
            if numel(spi)~=0
                meani = sum(d(spi).*bins(spi))./sum(d(spi));
            else
                meani = 0;
            end
            meanbk_arr(i) = meani;
        end

        %%% write stack profile data %%%
        if itype==1
            profscb100 = mean1;
            errscb100 = nanstd(meansub_arr)./sqrt(Nsub);
            profscbbk100 = nanmean(meanbk_arr);
            errscbbk100 = std(meanbk_arr);
        elseif itype==2
            profgcb100 = mean1;
            errgcb100 = nanstd(meansub_arr)./sqrt(Nsub);
            profgcbbk100 = nanmean(meanbk_arr);
            errgcbbk100 = std(meanbk_arr);
        elseif itype==3
            profsps100 = mean1;
            errsps100 = nanstd(meansub_arr)./sqrt(Nsub);
            profspsbk100 = nanmean(meanbk_arr);
            errspsbk100 = std(meanbk_arr);
        elseif itype==4
            profgps100 = mean1;
            errgps100 = nanstd(meansub_arr)./sqrt(Nsub);
            profgpsbk100 = nanmean(meanbk_arr);
            errgpsbk100 = std(meanbk_arr);
        end
    end
    
    %%%%%%%%%%%%%%%%%
    
    ihlprofdat.data(im).profscb = profscb;
    ihlprofdat.data(im).profscb_err= errscb;
    ihlprofdat.data(im).profgcb = profgcb;
    ihlprofdat.data(im).profgcb_err= errgcb;
    ihlprofdat.data(im).profsps = profsps;
    ihlprofdat.data(im).profsps_err= errsps;
    ihlprofdat.data(im).profgps = profgps;
    ihlprofdat.data(im).profgps_err= errgps;
    ihlprofdat.data(im).clipmax_arr = clipmax_arr;
    ihlprofdat.data(im).clipmin_arr = clipmin_arr;
    ihlprofdat.data100(im).profscb = profscb100;
    ihlprofdat.data100(im).profscb_err= errscb100;
    ihlprofdat.data100(im).profgcb = profgcb100;
    ihlprofdat.data100(im).profgcb_err = errgcb100;
    ihlprofdat.data100(im).profsps = profsps100;
    ihlprofdat.data100(im).profsps_err = errsps100;
    ihlprofdat.data100(im).profgps = profgps100;
    ihlprofdat.data100(im).profgps_err = errgps100;

    ihlprofdat.bk(im).profscb = profscbbk;
    ihlprofdat.bk(im).profscb_err= errscbbk;
    ihlprofdat.bk(im).profgcb = profgcbbk;
    ihlprofdat.bk(im).profgcb_err= errgcbbk;
    ihlprofdat.bk(im).profsps = profspsbk;
    ihlprofdat.bk(im).profsps_err= errspsbk;
    ihlprofdat.bk(im).profgps = profgpsbk;
    ihlprofdat.bk(im).profgps_err= errgpsbk;

    ihlprofdat.bk100(im).profscb = profscbbk100;
    ihlprofdat.bk100(im).profscb_err= errscbbk100;
    ihlprofdat.bk100(im).profgcb = profgcbbk100;
    ihlprofdat.bk100(im).profgcb_err= errgcbbk100;
    ihlprofdat.bk100(im).profsps = profspsbk100;
    ihlprofdat.bk100(im).profsps_err= errspsbk100;
    ihlprofdat.bk100(im).profgps = profgpsbk100;
    ihlprofdat.bk100(im).profgps_err= errgpsbk100;
end         
%%  get normalized profile
sp = [1];
for mc=1:3
    counts = histdat(mc).counts;
    countg = histdat(mc).countg;
    if counts<10 | countg<10
        continue
    end
    prof = ihlprofdat.data(mc).profgcb - ihlprofdat.bk(mc).profgcb;
    prof_err = ihlprofdat.data(mc).profgcb_err;
    prof_err1 = sqrt(ihlprofdat.data(mc).profgcb_err.^2 +...
        ihlprofdat.bk(mc).profgcb_err.^2);
    ihlprofdat.norm(mc).profgcb = prof;
    ihlprofdat.norm(mc).profgcb_err = prof_err;
    ihlprofdat.norm(mc).profgcb_err1 = prof_err1;

    prof = ihlprofdat.data(mc).profscb - ihlprofdat.bk(mc).profscb;
    prof_err = ihlprofdat.data(mc).profscb_err;
    prof_err1 = sqrt(ihlprofdat.data(mc).profscb_err.^2 +...
        ihlprofdat.bk(mc).profscb_err.^2);
    norm = mean(ihlprofdat.norm(mc).profgcb(sp))./mean(prof(sp));
    prof = prof.*norm;prof_err = prof_err.*norm;prof_err1 = prof_err1.*norm;
    ihlprofdat.norm(mc).profscb = prof;
    ihlprofdat.norm(mc).profscb_err = prof_err;
    ihlprofdat.norm(mc).profscb_err1 = prof_err1;
    
    prof = ihlprofdat.data(mc).profgps - ihlprofdat.bk(mc).profgps;
    prof_err = ihlprofdat.data(mc).profgps_err;
    prof_err1 = sqrt(ihlprofdat.data(mc).profgps_err.^2 +...
        ihlprofdat.bk(mc).profgps_err.^2);
    ihlprofdat.norm(mc).profgps = prof;
    ihlprofdat.norm(mc).profgps_err = prof_err;
    ihlprofdat.norm(mc).profgps_err1 = prof_err1;

    prof = ihlprofdat.data(mc).profsps - ihlprofdat.bk(mc).profsps;
    prof_err = ihlprofdat.data(mc).profsps_err;
    prof_err1 = sqrt(ihlprofdat.data(mc).profsps_err.^2 +...
        ihlprofdat.bk(mc).profsps_err.^2);
    norm = mean(ihlprofdat.norm(mc).profgps(sp))./mean(prof(sp));
    prof = prof.*norm;prof_err = prof_err.*norm;prof_err1 = prof_err1.*norm;
    ihlprofdat.norm(mc).profsps = prof;
    ihlprofdat.norm(mc).profsps_err = prof_err;
    ihlprofdat.norm(mc).profsps_err1 = prof_err1;
    
    if psf_model_stars
        norm = mean(ihlprofdat.norm(mc).profgcb(sp))./mean(ihlprofdat.psf_arr(sp));
        ihlprofdat.norm(mc).profscb = ihlprofdat.psf_arr.*norm;
        ihlprofdat.norm(mc).profscb_err = 0;
        ihlprofdat.norm(mc).profscb_err1 = 0;
        norm = mean(ihlprofdat.norm(mc).profgps(sp))./mean(ihlprofdat.psf_arr(sp));
        ihlprofdat.norm(mc).profsps = ihlprofdat.psf_arr.*norm;
        ihlprofdat.norm(mc).profsps_err = 0;
        ihlprofdat.norm(mc).profsps_err1 = 0;
    end
end
%%  get excess > 100 arcsec
sp = [1];
for mc=1:3
    counts = histdat(mc).counts;
    countg = histdat(mc).countg;
    if counts<10 | countg<10
        continue
    end
    profg = ihlprofdat.data100(mc).profgcb - ihlprofdat.bk100(mc).profgcb;
    profg_err = sqrt(ihlprofdat.data100(mc).profgcb_err.^2 +...
        ihlprofdat.bk100(mc).profgcb_err.^2);
    profs = ihlprofdat.data100(mc).profscb - ihlprofdat.bk100(mc).profscb;
    profs_err = sqrt(ihlprofdat.data100(mc).profscb_err.^2 +...
        ihlprofdat.bk100(mc).profscb_err.^2);
    norm = mean(ihlprofdat.norm(mc).profgcb(sp)) ...
        ./mean(ihlprofdat.norm(mc).profscb(sp));
    profs = profs.*norm;profs_err = profs_err.*norm;
    Ecb = profg - profs;
    Ecb_err = sqrt(profg_err.^2+profs_err.^2);
    if psf_model_stars
        Ecb = profg - mean(ihlprofdat.norm(mc).profscb(sp100));
        Ecb_err = profg_err;
    end
    profg = ihlprofdat.data100(mc).profgps - ihlprofdat.bk100(mc).profgps;
    profg_err = sqrt(ihlprofdat.data100(mc).profgps_err.^2 +...
        ihlprofdat.bk100(mc).profgps_err.^2);
    profs = ihlprofdat.data100(mc).profsps - ihlprofdat.bk100(mc).profsps;
    profs_err = sqrt(ihlprofdat.data100(mc).profsps_err.^2 +...
        ihlprofdat.bk100(mc).profsps_err.^2);
    norm = mean(ihlprofdat.norm(mc).profgps(sp)) ...
        ./mean(ihlprofdat.norm(mc).profsps(sp));
    profs = profs.*norm;profs_err = profs_err.*norm;
    Eps = profg - profs;
    Eps_err = sqrt(profg_err.^2+profs_err.^2);
    if psf_model_stars
        Eps = profg - mean(ihlprofdat.norm(mc).profsps(sp100));
        Eps_err = profg_err;
    end
    
    E = Ecb - Eps;
    E_err = sqrt(Ecb_err.^2+Eps_err.^2);
    
    ihlprofdat.excess(mc).diff100 = E;
    ihlprofdat.excess(mc).diff_err100 = E_err;   
end
%% get excess profile
for mc=1:3
    counts = histdat(mc).counts;
    countg = histdat(mc).countg;
    if counts<10 | countg<10
        continue
    end

    diffcb = ihlprofdat.norm(mc).profgcb - ihlprofdat.norm(mc).profscb;
    diffcb_err = sqrt(ihlprofdat.norm(mc).profgcb_err.^2 +...
                      ihlprofdat.norm(mc).profscb_err.^2);
    diffcb_err1 = sqrt(ihlprofdat.norm(mc).profgcb_err1.^2 +...
                      ihlprofdat.norm(mc).profscb_err1.^2);
    diffps = ihlprofdat.norm(mc).profgps - ihlprofdat.norm(mc).profsps;
    diffps_err = sqrt(ihlprofdat.norm(mc).profgps_err.^2 +...
                      ihlprofdat.norm(mc).profsps_err.^2);
    diffps_err1 = sqrt(ihlprofdat.norm(mc).profgps_err1.^2 +...
                      ihlprofdat.norm(mc).profsps_err1.^2);                  
    diff = diffcb - diffps;
    % there are some numerical issue, so force it to 0
    diffcb(1) = 0;
    diffps(1) = 0;
    diff(1) = 0;
    diff_err = sqrt(diffcb_err.^2 + diffps_err.^2);
    diff_err1 = sqrt(diffcb_err1.^2 + diffps_err1.^2);
    ihlprofdat.excess(mc).diffcb = diffcb;
    ihlprofdat.excess(mc).diffcb_err = diffcb_err;
    ihlprofdat.excess(mc).diffcb_err1 = diffcb_err1;
    ihlprofdat.excess(mc).diffps = diffps;
    ihlprofdat.excess(mc).diffps_err = diffps_err;
    ihlprofdat.excess(mc).diffps_err1 = diffps_err1;
    ihlprofdat.excess(mc).diff = diff;
    ihlprofdat.excess(mc).diff_err = diff_err;
    ihlprofdat.excess(mc).diff_err1 = diff_err1;
end
%%
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if strcmp(field,'SIDES')
    if rvir==1
        save(sprintf('%s/ihlprofdatsim_hist%d',loaddir,f_ihl*100),'ihlprofdat');
    else
        save(sprintf('%s/ihlprofdatsim_hist%d_rv%d',...
            loaddir,f_ihl*100,rvir),'ihlprofdat');
    end
elseif strcmp(field,'HSC')
    if rvir==1
        save(sprintf('%s/ihlprofdathsc_%s_hist%d',loaddir,name,f_ihl*100),...
            'ihlprofdat');
    else
        save(sprintf('%s/ihlprofdathsc_%s_hist%d_rv%d',...
            name,loaddir,f_ihl*100,rvir),'ihlprofdat');
    end
end
%%
return