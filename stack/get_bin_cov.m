flight=40030;
inst=1;
ifield=8;
%%
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%shistdat_%s',loaddir,dt.name),'histdat');
load(sprintf('%s/%s_ihlprofdat_hist',loaddir,dt.name),'ihlprofdat');

nbins = 25;

figure
setwinsize(gcf,900,1200);

mc = 0;
for im=10:12
    mc = mc + 1;
    m_min = ihlprofdat.data(mc).m_min;
    m_max = ihlprofdat.data(mc).m_min;
    clipmax_arr = ihlprofdat.data(mc).clipmax_arr;
    clipmin_arr = ihlprofdat.data(mc).clipmin_arr;
    subprof_arr = zeros([4,10,nbins]);
    %%% get the stacking profile %%%
    for ibin=1:nbins
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
            
            clipmax = clipmax_arr(itype,ibin);
            clipmin = clipmin_arr(itype,ibin);
            
            %%% get stacking err from sub stack %%%
            Nsub = size(d_all,1);
            for isub=1:Nsub
                d = d_all(isub,:,:);
                d = reshape(d, size(bins));
                spi = find((d~=0) & (bins < clipmax) & (bins > clipmin));
                meani = sum(d(spi).*bins(spi))./sum(d(spi));
                subprof_arr(itype,isub,ibin) = meani;
            end
        end
    end
    
    cov_all = zeros([nbins,nbins,4]);
    corr_all = zeros([nbins,nbins,4]);
    for i=1:nbins
    for j=1:nbins
        for itype=1:4
            profi = subprof_arr(itype,:,i);
            profj = subprof_arr(itype,:,j);
            sp = find((profj==profj) & (profi==profi));
            profi = profi(sp);
            profj = profj(sp);
            cov = mean(profi.*profj)-mean(profi)*mean(profj);
            cov_all(i,j,itype)=cov;
            sigi = sqrt(mean(profi.^2)-mean(profi).^2);
            sigj = sqrt(mean(profj.^2)-mean(profj).^2);
            corr_all(i,j,itype)=  cov./sigi./sigj;
        end
    end
    end
    subplot(4,3,mc)
    imageclip(corr_all(:,:,1));
    title('CIBER stars');
    caxis([-1,1]);
    subplot(4,3,mc+3)
    imageclip(corr_all(:,:,2));
    title('CIBER gals');
    caxis([-1,1]);
    subplot(4,3,mc+6)
    imageclip(corr_all(:,:,3));
    title('PanSTARRS stars');
    caxis([-1,1]);
    subplot(4,3,mc+9)
    imageclip(corr_all(:,:,4));
    title('PanSTARRS gals');
    caxis([-1,1]);
end
%%
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);

f_ihl=0;
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%shistdatsim_%s%d',loaddir,dt.name,f_ihl*100),'histdat');
load(sprintf('%s%s_ihlprofdatsim_hist%d',loaddir,dt.name,f_ihl*100),'ihlprofdat');

nbins = 25;

figure
setwinsize(gcf,900,1200);

mc = 0;
for im=10:12
    mc = mc + 1;
    m_min = ihlprofdat.data(mc).m_min;
    m_max = ihlprofdat.data(mc).m_min;
    clipmax_arr = ihlprofdat.data(mc).clipmax_arr;
    clipmin_arr = ihlprofdat.data(mc).clipmin_arr;
    subprof_arr = zeros([4,10,nbins]);
    %%% get the stacking profile %%%
    for ibin=1:nbins
        for itype=1:4
            
            %%% get data %%%
            if itype==1
                bins = binedges2bins(histdat(mc).Ibinedges_cb);
                d_all = histdat(mc).histcbs(:,ibin,:);
                name = 'CIBER stars';
            elseif itype==2
                bins = binedges2bins(histdat(mc).Ibinedges_cb);
                d_all = histdat(mc).histcbg(:,ibin,:);
                name = 'CIBER gals';
            elseif itype==3
                bins = binedges2bins(histdat(mc).Ibinedges_ps);
                d_all = histdat(mc).histpss(:,ibin,:);
                name = 'Sim stars';
            elseif itype==4
                bins = binedges2bins(histdat(mc).Ibinedges_ps);
                d_all = histdat(mc).histpsg(:,ibin,:);
                name = 'Sim gals';
            end
            
            clipmax = clipmax_arr(itype,ibin);
            clipmin = clipmin_arr(itype,ibin);
            
            %%% get stacking err from sub stack %%%
            Nsub = size(d_all,1);
            for isub=1:Nsub
                d = d_all(isub,:,:);
                d = reshape(d, size(bins));
                spi = find((d~=0) & (bins < clipmax) & (bins > clipmin));
                meani = sum(d(spi).*bins(spi))./sum(d(spi));
                subprof_arr(itype,isub,ibin) = meani;
            end
        end
    end
    
    cov_all = zeros([nbins,nbins,4]);
    corr_all = zeros([nbins,nbins,4]);
    for i=1:nbins
    for j=1:nbins
        for itype=1:4
            profi = subprof_arr(itype,:,i);
            profj = subprof_arr(itype,:,j);
            sp = find((profj==profj) & (profi==profi));
            profi = profi(sp);
            profj = profj(sp);
            cov = mean(profi.*profj)-mean(profi)*mean(profj);
            cov_all(i,j,itype)=cov;
            sigi = sqrt(mean(profi.^2)-mean(profi).^2);
            sigj = sqrt(mean(profj.^2)-mean(profj).^2);
            corr_all(i,j,itype)=  cov./sigi./sigj;
        end
    end
    end
    subplot(4,3,mc)
    imageclip(corr_all(:,:,1));
    title('CIBER stars');
    caxis([-1,1]);
    subplot(4,3,mc+3)
    imageclip(corr_all(:,:,2));
    title('CIBER gals');
    caxis([-1,1]);
    subplot(4,3,mc+6)
    imageclip(corr_all(:,:,3));
    title('PanSTARRS stars');
    caxis([-1,1]);
    subplot(4,3,mc+9)
    imageclip(corr_all(:,:,4));
    title('PanSTARRS gals');
    caxis([-1,1]);
end