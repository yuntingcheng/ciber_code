function ihl_stack_compileFF(flight,inst,ifield)

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

m_min_arr = stackmapdat(ifield).m_min_arr;
m_max_arr = stackmapdat(ifield).m_max_arr;
counts_arr = stackmapdat(ifield).count_stacks_arr;
countg_arr = stackmapdat(ifield).count_stackg_arr;

%%% get the PSF %%%
bestparam = fitpsfdat(ifield).bestparam;
beta = bestparam(1);
rc = bestparam(2);

radmap = make_radius_map(zeros(2*npix+1),npix,npix).*0.7;
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

    %%% get the background profile %%%
%     profcb = zeros(size(r_arr));
%     profcb_err = zeros(size(r_arr));
%     profps = zeros(size(r_arr));
%     profps_err = zeros(size(r_arr));
%     for iter=1:nsim
%         M = csvread(strcat(stackdir,'bk_ps/',dt.name,'_randprof',...
%             '_',num2str(counts_arr(im)),'_',num2str(iter),'.csv'));
% 
%         profcb = profcb + M(:,2)';
%         profps = profps + M(:,4)';
%         profcb_err = profcb_err + M(:,2)'.^2;
%         profps_err = profps_err + M(:,4)'.^2;
%     end
%     profcb_err = profcb_err - profcb;
%     profps_err = profps_err - profps;   
%     profcb = profcb./nsim;
%     profps = profps./nsim;
%     profcb_err = profcb_err./nsim;
%     profps_err = profps_err./nsim;

%     ihlprofdat.bk(im).profscb = profcb;
%     ihlprofdat.bk(im).profscb_err = profcb_err;
%     ihlprofdat.bk(im).profsps = profps;
%     ihlprofdat.bk(im).profsps_err = profps_err;

%     profcb = zeros(size(r_arr));
%     profcb_err = zeros(size(r_arr));
%     profps = zeros(size(r_arr));
%     profps_err = zeros(size(r_arr));
% 
%     for iter=1:nsim
%         M = csvread(strcat(stackdir,'bk_ps/',dt.name,'_randprof',...
%             '_',num2str(countg_arr(im)),'_',num2str(iter),'.csv'));
%         profcb = profcb + M(:,2)';
%         profps = profps + M(:,4)';
%         profcb_err = profcb_err + M(:,2)'.^2;
%         profps_err = profps_err + M(:,4)'.^2;
%     end
% 
%     profcb_err = profcb_err - profcb;
%     profps_err = profps_err - profps;   
%     profcb = profcb./nsim;
%     profps = profps./nsim;
%     profcb_err = profcb_err./nsim;
%     profps_err = profps_err./nsim;
% 
%     ihlprofdat.bk(im).profgcb = profcb;
%     ihlprofdat.bk(im).profgcb_err = profcb_err;
%     ihlprofdat.bk(im).profgps = profps;
%     ihlprofdat.bk(im).profgps_err = profps_err;

    %%% get the stacking maps %%%
    datascb=fitsread(strcat(stackdir,'/FFerr_test/',dt.name,...
        '_stampers',num2str(m_min),'_',num2str(m_max),'FF.fits'));
    datagcb=fitsread(strcat(stackdir,'/FFerr_test/',dt.name,...
        '_stamperg',num2str(m_min),'_',num2str(m_max),'FF.fits'));                    
    datasps=fitsread(strcat(stackdir,'/FFerr_test/',dt.name,...
        '_stampers',num2str(m_min),'_',num2str(m_max),'.fits'));
    datagps=fitsread(strcat(stackdir,'/FFerr_test/',dt.name,...
        '_stamperg',num2str(m_min),'_',num2str(m_max),'.fits'));                    
    mdatas=load(strcat(stackdir,'/FFerr_test/',dt.name,...
        '_hitmaps',num2str(m_min),'_',num2str(m_max),'.mat'));
    mdatas = mdatas.hitmap;
    mdatag=load(strcat(stackdir,'/FFerr_test/',dt.name,...
        '_hitmapg',num2str(m_min),'_',num2str(m_max),'.mat'));
    mdatag = mdatag.hitmap;
    stackscb=datascb./mdatas;
    stackgcb=datagcb./mdatag;
    stacksps=datasps./mdatas;
    stackgps=datagps./mdatag;
    
    mmocks = zeros(size(mdatas));
    
    %%% plot the stacking maps %%%
    figure
    setwinsize(gcf,800,600)
    subplot(2,2,1)
    imageclip(stackscb);
    title(strcat('PanSTARRS + FF error',{' '},...
        'stack',{' '}, num2str(counts),{' '},'stars'));
    vs = caxis;
    subplot(2,2,2)
    imageclip(stackgcb);
    title(strcat('PanSTARRS + FF error',{' '},...
        'stack',{' '}, num2str(countg),{' '},'gals'));
    vg = caxis;  
    subplot(2,2,3)
    imageclip(stacksps);
    caxis(vs);
    title(strcat('PanSTARRS',{' '}, 'stack',{' '}, num2str(counts),{' '},'stars'));
    subplot(2,2,4)
    imageclip(stackgps);
    caxis(vg);
    title(strcat('PanSTARRS',{' '}, 'stack',{' '}, num2str(countg),{' '},'gals'));

    suptitle(strcat({' '},num2str(m_min),...
        '<mAB(y band)<',num2str(m_max)));
    if ismember(im,10:13)
        savename=strcat(pltsavedir,dt.name,'_stackmaps',...
            num2str(m_min),'_',num2str(m_max),'FF');
        print(savename,'-dpng');
    end

    %%% get the stacking profile %%%
    prof = radial_prof(stackscb,ones(2*npix+1),npix+1,npix+1,1,25,...
        'sig',5,'iter_clip',0,'weight',mdatas + mmocks);
%     bksub = prof.prof - ihlprofdat.bk(im).profscb;
%     err = sqrt(prof.err.^2 + ihlprofdat.bk(im).profscb_err.^2);
    bksub = prof.prof;
    err = prof.err;
    ihlprofdat.data(im).profscb = bksub;
    ihlprofdat.data(im).profscb_err= err;
    norm = stackscb(npix+1,npix+1);  
    ihlprofdat.norm(im).profscb = bksub./norm;
    ihlprofdat.norm(im).profscb_err= err./norm;

    prof = radial_prof(stacksps,ones(2*npix+1),npix+1,npix+1,1,25,...
        'sig',5,'iter_clip',0,'weight',mdatas + mmocks);
%     bksub = prof.prof - ihlprofdat.bk(im).profsps;
%     err = sqrt(prof.err.^2 + ihlprofdat.bk(im).profsps_err.^2);
    bksub = prof.prof;
    err = prof.err;
    ihlprofdat.data(im).profsps = bksub;
    ihlprofdat.data(im).profsps_err= err;
    norm = stacksps(npix+1,npix+1);  
    ihlprofdat.norm(im).profsps = bksub./norm;
    ihlprofdat.norm(im).profsps_err= err./norm;
    
    prof = radial_prof(stackgcb,ones(2*npix+1),npix+1,npix+1,1,25,...
        'sig',5,'iter_clip',0,'weight',mdatag);
%     bksub = prof.prof - ihlprofdat.bk(im).profgcb;
%     err = sqrt(prof.err.^2 + ihlprofdat.bk(im).profgcb_err.^2);
    bksub = prof.prof;
    err = prof.err;
    ihlprofdat.data(im).profgcb = bksub;
    ihlprofdat.data(im).profgcb_err= err;
    norm = stackgcb(npix+1,npix+1);  
    ihlprofdat.norm(im).profgcb = bksub./norm;
    ihlprofdat.norm(im).profgcb_err= err./norm;
    
    prof = radial_prof(stackgps,ones(2*npix+1),npix+1,npix+1,1,25,...
        'sig',5,'iter_clip',0,'weight',mdatag);
%     bksub = prof.prof - ihlprofdat.bk(im).profgps;
%     err = sqrt(prof.err.^2 + ihlprofdat.bk(im).profgps_err.^2);
    bksub = prof.prof;
    err = prof.err;
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
save(sprintf('%s/%s_ihlprofdatFF',loaddir,dt.name),'ihlprofdat');

return