function ihl_stack_compile_old(flight,inst,ifield)

mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);

stackdir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

psfdir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(psfdir,'fitpsfdat'),'fitpsfdat');

pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

npix = 1200;
nsim = 100;

m_min_arr = stackmapdat(ifield).m_min_arr;
m_max_arr = stackmapdat(ifield).m_max_arr;
counts_arr = stackmapdat(ifield).count_stacks_arr;
countg_arr = stackmapdat(ifield).count_stackg_arr;

%%% get the PSF %%%
bestparam = fitpsfdat(ifield).bestparam;
% A=bestparam(1);
% B=bestparam(2);
% sig=bestparam(3);
% r0=bestparam(4);
% alpha=bestparam(5);
beta = bestparam(1);
rc = bestparam(2);

radmap = make_radius_map(zeros(2*npix+1),npix,npix).*0.7;
%psfmap = A*exp(-radmap.^2./2./sig^2)+B./(1+(radmap./r0).^alpha);
psfmap = (1 + (radmap/rc).^2).^(-3.*beta./2);

profile = radial_prof(psfmap,ones(2*npix+1),npix+1,npix+1,1,25);
r_arr=profile.r*0.7;
profpsf_arr=(profile.prof)./profile.prof(1);

%%% write the basic info %%%
ihlprofdat.r_arr = r_arr;
ihlprofdat.psf_arr = profpsf_arr;
ihlprofdat.m_min_arr = m_min_arr;
ihlprofdat.m_max_arr = m_max_arr;
ihlprofdat.counts_arr = counts_arr;
ihlprofdat.countg_arr = countg_arr;

for im=10:12%1:numel(m_min_arr)

    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    counts = counts_arr(im);
    countg = countg_arr(im);
    countd = countg - counts;

    %%% get the background profile %%%
    if counts_arr(im)~=0
        profcb = zeros(size(r_arr));
        profcb_err = zeros(size(r_arr));
        profps = zeros(size(r_arr));
        profps_err = zeros(size(r_arr));
        for iter=1:nsim
            M = csvread(strcat(stackdir,'bk_ps/',dt.name,'_randprof',...
                '_',num2str(counts_arr(im)),'_',num2str(iter),'.csv'));

            profcb = profcb + M(:,2)';
            profps = profps + M(:,4)';
            profcb_err = profcb_err + M(:,2)'.^2;
            profps_err = profps_err + M(:,4)'.^2;
        end
        profcb_err = profcb_err - profcb;
        profps_err = profps_err - profps;   
        profcb = profcb./nsim;
        profps = profps./nsim;
        profcb_err = profcb_err./nsim;
        profps_err = profps_err./nsim;
    else
        profcb = zeros(1,numel(r_arr));
        profps = zeros(1,numel(r_arr));
        profcb_err = zeros(1,numel(r_arr));
        profps_err = zeros(1,numel(r_arr));        
    end
    ihlprofdat.bk(im).profscb = profcb;
    ihlprofdat.bk(im).profscb_err = profcb_err;
    ihlprofdat.bk(im).profsps = profps;
    ihlprofdat.bk(im).profsps_err = profps_err;

    if countg_arr(im)~=0
        profcb = zeros(size(r_arr));
        profcb_err = zeros(size(r_arr));
        profps = zeros(size(r_arr));
        profps_err = zeros(size(r_arr));

        for iter=1:nsim
            M = csvread(strcat(stackdir,'bk_ps/',dt.name,'_randprof',...
                '_',num2str(countg_arr(im)),'_',num2str(iter),'.csv'));
            profcb = profcb + M(:,2)';
            profps = profps + M(:,4)';
            profcb_err = profcb_err + M(:,2)'.^2;
            profps_err = profps_err + M(:,4)'.^2;
        end
        
        profcb_err = profcb_err - profcb;
        profps_err = profps_err - profps;   
        profcb = profcb./nsim;
        profps = profps./nsim;
        profcb_err = profcb_err./nsim;
        profps_err = profps_err./nsim;
    else
        profcb = zeros(1,numel(r_arr));
        profps = zeros(1,numel(r_arr));
        profcb_err = zeros(1,numel(r_arr));
        profps_err = zeros(1,numel(r_arr));        
    end
    ihlprofdat.bk(im).profgcb = profcb;
    ihlprofdat.bk(im).profgcb_err = profcb_err;
    ihlprofdat.bk(im).profgps = profps;
    ihlprofdat.bk(im).profgps_err = profps_err;

    
    
    %%% get the stacking maps %%%
    datascb=fitsread(strcat(stackdir,'/ciber_ps/',dt.name,...
        '_stampers',num2str(m_min),'_',num2str(m_max),'.fits'));
    datagcb=fitsread(strcat(stackdir,'/ciber_ps/',dt.name,...
        '_stamperg',num2str(m_min),'_',num2str(m_max),'.fits'));                    
    datasps=fitsread(strcat(stackdir,'/panstarrs_ps/',dt.name,...
        '_stampers',num2str(m_min),'_',num2str(m_max),'.fits'));
    datagps=fitsread(strcat(stackdir,'/panstarrs_ps/',dt.name,...
        '_stamperg',num2str(m_min),'_',num2str(m_max),'.fits'));                    
    mdatas=load(strcat(stackdir,'/ciber_ps/',dt.name,...
        '_hitmaps',num2str(m_min),'_',num2str(m_max),'.mat'));
    mdatas = mdatas.hitmap;
    mdatag=load(strcat(stackdir,'/ciber_ps/',dt.name,...
        '_hitmapg',num2str(m_min),'_',num2str(m_max),'.mat'));
    mdatag = mdatag.hitmap;
    stackscb=datascb./mdatas;
    stackgcb=datagcb./mdatag;
    stacksps=datasps./mdatas;
    stackgps=datagps./mdatag;
    
    %{
    %%% get the mock stacking map %%%
    if 1==2% disable the mock%countd > 0
        mockscb=fitsread(strcat(stackdir,'/ciber_ps/',dt.name,...
            '_stampers',num2str(m_min),'_',num2str(m_max),'mock.fits'));
        mocksps=fitsread(strcat(stackdir,'/panstarrs_ps/',dt.name,...
            '_stampers',num2str(m_min),'_',num2str(m_max),'mock.fits'));
        mmocks=load(strcat(stackdir,'/ciber_ps/',dt.name,...
            '_hitmaps',num2str(m_min),'_',num2str(m_max),'mock.mat'));
        mmocks = mmocks.hitmap;
        stackmockscb=mockscb./mmocks;
        stackmocksps=mocksps./mmocks;
        
        %%% plot the stacking maps %%%
        figure
        setwinsize(gcf,800,600)
        subplot(2,2,1)
        imageclip(stackscb(npix-100:npix+100,npix-100:npix+100));
        title(strcat('dataCB',{' '}, 'stack',{' '}, num2str(counts),{' '},'stars'));
        subplot(2,2,2)
        imageclip(stackmockscb(npix-100:npix+100,npix-100:npix+100));
        title(strcat('mockCB',{' '}, 'stack',{' '}, num2str(countd),{' '},'stars'));
        subplot(2,2,3)
        imageclip(stacksps(npix-100:npix+100,npix-100:npix+100));
        title(strcat('dataPS',{' '}, 'stack',{' '}, num2str(counts),{' '},'stars'));
        subplot(2,2,4)
        imageclip(stackmocksps(npix-100:npix+100,npix-100:npix+100));
        title(strcat('mockPS',{' '}, 'stack',{' '}, num2str(countd),{' '},'stars'));

        suptitle(strcat('mock stacking PanSTARRS',{' '},num2str(m_min),...
            '<mAB(y band)<',num2str(m_max)));

%         savename=strcat(pltsavedir,dt.name,'_mockstackmaps',...
%             num2str(m_min),'_',num2str(m_max));
%         print(savename,'-dpng');
        drawnow 
        
        stackscb=(datascb + mockscb)./(mdatas + mmocks);
        stacksps=(datasps + mocksps)./(mdatas + mmocks);
    
    else
        mmocks = zeros(size(mdatas));
    end
    %}
    mmocks = zeros(size(mdatas));
    
    %%% plot the stacking maps %%%
    figure
    setwinsize(gcf,800,600)
    subplot(2,2,1)
    imageclip(stackscb(npix-100:npix+100,npix-100:npix+100));
    title(strcat('CIBER',{' '}, 'stack',{' '}, num2str(counts),{' '},'stars'));
    subplot(2,2,2)
    imageclip(stackgcb(npix-100:npix+100,npix-100:npix+100));
    title(strcat('CIBER',{' '}, 'stack',{' '}, num2str(countg),{' '},'gals'));
    subplot(2,2,3)
    imageclip(stacksps(npix-100:npix+100,npix-100:npix+100));
    title(strcat('PanSTARRS',{' '}, 'stack',{' '}, num2str(counts),{' '},'stars'));
    subplot(2,2,4)
    imageclip(stackgps(npix-100:npix+100,npix-100:npix+100));
    title(strcat('PanSTARRS',{' '}, 'stack',{' '}, num2str(countg),{' '},'gals'));

    suptitle(strcat('PanSTARRS',{' '},num2str(m_min),...
        '<mAB(y band)<',num2str(m_max)));
    
    if ismember(im,10:13)
        savename=strcat(pltsavedir,dt.name,'_stackmaps',...
            num2str(m_min),'_',num2str(m_max));
        print(savename,'-dpng');
    end
    close
    
    %%% get the stacking profile %%%
    prof = radial_prof(stackscb,ones(2*npix+1),npix+1,npix+1,1,25,...
        'sig',5,'iter_clip',0,'weight',mdatas + mmocks);
    bksub = prof.prof - ihlprofdat.bk(im).profscb;
    err = sqrt(prof.err.^2 + ihlprofdat.bk(im).profscb_err.^2);
%     if countd > 0
%         bksub = prof.prof - ihlprofdat.bk(im).profgcb;
%         err = sqrt(prof.err.^2 + ihlprofdat.bk(im).profgcb_err.^2);
%     end
    ihlprofdat.data(im).profscb = bksub;
    ihlprofdat.data(im).profscb_err= err;
    norm = stackscb(npix+1,npix+1);  
    ihlprofdat.norm(im).profscb = bksub./norm;
    ihlprofdat.norm(im).profscb_err= err./norm;

    prof = radial_prof(stacksps,ones(2*npix+1),npix+1,npix+1,1,25,...
        'sig',5,'iter_clip',0,'weight',mdatas + mmocks);
    bksub = prof.prof - ihlprofdat.bk(im).profsps;
    err = sqrt(prof.err.^2 + ihlprofdat.bk(im).profsps_err.^2);
%     if countd > 0
%         bksub = prof.prof - ihlprofdat.bk(im).profgps;
%         err = sqrt(prof.err.^2 + ihlprofdat.bk(im).profgps_err.^2);
%     end
    ihlprofdat.data(im).profsps = bksub;
    ihlprofdat.data(im).profsps_err= err;
    norm = stacksps(npix+1,npix+1);  
    ihlprofdat.norm(im).profsps = bksub./norm;
    ihlprofdat.norm(im).profsps_err= err./norm;
    
    prof = radial_prof(stackgcb,ones(2*npix+1),npix+1,npix+1,1,25,...
        'sig',5,'iter_clip',0,'weight',mdatag);
    bksub = prof.prof - ihlprofdat.bk(im).profgcb;
    err = sqrt(prof.err.^2 + ihlprofdat.bk(im).profgcb_err.^2);
    ihlprofdat.data(im).profgcb = bksub;
    ihlprofdat.data(im).profgcb_err= err;
    norm = stackgcb(npix+1,npix+1);  
    ihlprofdat.norm(im).profgcb = bksub./norm;
    ihlprofdat.norm(im).profgcb_err= err./norm;
    
    prof = radial_prof(stackgps,ones(2*npix+1),npix+1,npix+1,1,25,...
        'sig',5,'iter_clip',0,'weight',mdatag);
    bksub = prof.prof - ihlprofdat.bk(im).profgps;
    err = sqrt(prof.err.^2 + ihlprofdat.bk(im).profgps_err.^2);
    ihlprofdat.data(im).profgps = bksub;
    ihlprofdat.data(im).profgps_err= err;
    norm = stackgps(npix+1,npix+1);  
    ihlprofdat.norm(im).profgps = bksub./norm;
    ihlprofdat.norm(im).profgps_err= err./norm;
    
    %%% get the excess profile %%%
    diffcb = ihlprofdat.norm(im).profgcb - ihlprofdat.norm(im).profscb;
    diffcb_err = sqrt(ihlprofdat.norm(im).profgcb_err.^2 +...
                      ihlprofdat.norm(im).profscb_err.^2);
    diffps = ihlprofdat.norm(im).profgps - ihlprofdat.norm(im).profsps;
    diffps_err = sqrt(ihlprofdat.norm(im).profgps_err.^2 +...
                      ihlprofdat.norm(im).profsps_err.^2);
    diff = diffcb - diffps;
    diff_err = sqrt(diffcb_err.^2 + diffps_err.^2);
    ihlprofdat.excess(im).diffcb = diffcb;
    ihlprofdat.excess(im).diffcb_err = diffcb_err;
    ihlprofdat.excess(im).diffps = diffps;
    ihlprofdat.excess(im).diffsp_err = diffps_err;
    ihlprofdat.excess(im).diff = diff;
    ihlprofdat.excess(im).diff_err = diff_err;
    
end

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
save(sprintf('%s/%s_ihlprofdat',loaddir,dt.name),'ihlprofdat');

return