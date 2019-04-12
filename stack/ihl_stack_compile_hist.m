function ihl_stack_compile_hist(flight,inst,ifield)
%%
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%shistdat_%s',loaddir,dt.name),'histdat');

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

%%% get the PSF params %%%
bestparam = fitpsfdat(ifield).bestparam;
beta = bestparam(1);
rc = bestparam(2);
radmap = make_radius_map(zeros(2*dx+1),dx,dx).*0.7;
psfmap = (1 + (radmap/rc).^2).^(-3.*beta./2);
profile = radial_prof(psfmap,ones(2*dx+1),dx+1,dx+1,1,nbins);
r_arr=profile.r*0.7;
profpsf_arr=(profile.prof)./profile.prof(1);
ihlprofdat.r_arr = r_arr;
ihlprofdat.binedges = profile.binedges;
ihlprofdat.psf_arr = profpsf_arr;
%ihlprofdat.A = get_profile_mkk(flight,inst,ifield,dx,ihlprofdat.binedges);

rbinedges = profile.binedges * 0.7;
rbins = binedges2bins(rbinedges);
%%
mc = 0;
for im=10:12
 
    mc = mc + 1;
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    counts = histdat(im).counts;
    countg = histdat(im).countg;
    ihlprofdat.data(mc).m_min = m_min;
    ihlprofdat.data(mc).m_max = m_max;
    ihlprofdat.data(mc).counts = counts;
    ihlprofdat.data(mc).countg = countg;
    
    %%% find the index of S & G in bk stack file %%%
%     load(strcat(bkdir,'bk_ps/',dt.name,'_histbkdat','_1'),'histbkdat');
%     
%     for i=1:size(histbkdat,2)
%         if histbkdat(i).N == counts
%             bkidx_s = i;
%         end
%         if histbkdat(i).N == countg
%             bkidx_g = i;
%         end  
%     end
    
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
                name = 'CIBER stars';
            elseif itype==2
                bins = binedges2bins(histdat(im).Ibinedges_cb);
                d_all = histdat(im).histcbg(:,ibin,:);
                name = 'CIBER gals';
            elseif itype==3
                bins = binedges2bins(histdat(im).Ibinedges_ps);
                d_all = histdat(im).histpss(:,ibin,:);
                name = 'Sim stars';
            elseif itype==4
                bins = binedges2bins(histdat(im).Ibinedges_ps);
                d_all = histdat(im).histpsg(:,ibin,:);
                name = 'Sim gals';
            end
            d = sum(d_all);
            d = reshape(d, size(bins));
            
            %%% get stacking profile %%%
            sp0 = find(d~=0);
            mean0 = sum(d.*bins)./sum(d);
            Q1 = quartile_from_hist(bins,d,0.25);
            Q3 = quartile_from_hist(bins,d,0.75);
            IQR = Q3 - Q1;
            clipmax = Q3+3*IQR;
            clipmin = Q1-3*IQR;
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
                vline(mean0,'b-');
                vline(mean1,'r-');
                vline(clipmax,'r--');
                vline(clipmin,'r--');
                vline(mean1 + std(meansub_arr)./sqrt(Nsub)./2 ,'k--');
                vline(mean1 - std(meansub_arr)./sqrt(Nsub)./2 ,'k--');
                xlabel('I [nW/m^2/sr]');
                ylabel('counts');
                title(name);
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
            
%             for i=1:nsim
%                 load(strcat(bkdir,'bk_ps/',dt.name,'_histbkdat',...
%                     '_',num2str(i)),'histbkdat');
%                 if itype==1
%                     d= histbkdat(bkidx_s).histcb(ibin,:);
%                 elseif itype==2
%                     d= histbkdat(bkidx_g).histcb(ibin,:);
%                 elseif itype==3
%                     d= histbkdat(bkidx_s).histps(ibin,:);
%                 elseif itype==4
%                     d= histbkdat(bkidx_g).histps(ibin,:);
%                 end
%                 d = reshape(d, size(bins));
%                 spi = find((d~=0) & (bins < clipmax) & (bins > clipmin));
%                 if numel(spi)~=0
%                     meani = sum(d(spi).*bins(spi))./sum(d(spi));
%                 else
%                     meani = 0;
%                 end
%                 meanbk_arr(i) = meani;
%             end
            
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
    
    ihlprofdat.data(mc).profscb = profscb;
    ihlprofdat.data(mc).profscb_err= errscb;
    ihlprofdat.data(mc).profgcb = profgcb;
    ihlprofdat.data(mc).profgcb_err= errgcb;
    ihlprofdat.data(mc).profsps = profsps;
    ihlprofdat.data(mc).profsps_err= errsps;
    ihlprofdat.data(mc).profgps = profgps;
    ihlprofdat.data(mc).profgps_err= errgps;

    ihlprofdat.bk(mc).profscb = profscbbk;
    ihlprofdat.bk(mc).profscb_err= errscbbk;
    ihlprofdat.bk(mc).profgcb = profgcbbk;
    ihlprofdat.bk(mc).profgcb_err= errgcbbk;
    ihlprofdat.bk(mc).profsps = profspsbk;
    ihlprofdat.bk(mc).profsps_err= errspsbk;
    ihlprofdat.bk(mc).profgps = profgpsbk;
    ihlprofdat.bk(mc).profgps_err= errgpsbk;

end         
%%  get normalized profile
mc = 0;
sp = [1];
for im=10:12
    mc = mc + 1;

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
    prof = prof.*norm;prof_err = prof_err.*norm;
    ihlprofdat.norm(mc).profsps = prof;
    ihlprofdat.norm(mc).profsps_err = prof_err;
    ihlprofdat.norm(mc).profsps_err1 = prof_err1;
       
end
%% get excess profile
mc = 0;
for im=1:3
    mc = mc + 1;
    diffcb = ihlprofdat.norm(mc).profgcb - ihlprofdat.norm(mc).profscb;
    diffcb_err = sqrt(ihlprofdat.norm(mc).profgcb_err.^2 +...
                      ihlprofdat.norm(mc).profscb_err.^2);
    diffcb_err1 = sqrt(ihlprofdat.norm(mc).profgcb_err1.^2 +...
                      ihlprofdat.norm(mc).profscb_err1.^2);
    diffps = ihlprofdat.norm(im).profgps - ihlprofdat.norm(im).profsps;
    diffps_err = sqrt(ihlprofdat.norm(mc).profgps_err.^2 +...
                      ihlprofdat.norm(mc).profsps_err.^2);
    diffps_err1 = sqrt(ihlprofdat.norm(mc).profgps_err1.^2 +...
                      ihlprofdat.norm(mc).profsps_err1.^2);                  
    diff = diffcb - diffps;
    % there are some numer issue, so force it to 0
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

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
save(sprintf('%s/%s_ihlprofdat_hist',loaddir,dt.name),'ihlprofdat');

return