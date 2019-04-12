function ihl_stack_compile(flight,inst,ifield)
%%
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);

stackdir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));
pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

psfdir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(psfdir,'fitpsfdat'),'fitpsfdat');

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/%s_mcerrdat',loaddir,dt.name),'mcerrdat');
load(sprintf('%s/stackcountdat_%s',loaddir,dt.name),'stackcountdat');

nsim = 50;
dx = 1200;
nbins = 25;

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

%%% get the stack counts %%%
counts_arr = zeros(1,3);
countg_arr = zeros(1,3);
for im = 10:12        
    counts_arr(im-9) = stackcountdat(im).counts;
    countg_arr(im-9) = stackcountdat(im).countg;
end

%%% get the background profile %%%
countbg = max([counts_arr, countg_arr])*3;
profcb = zeros(size(r_arr));
profcb_err = zeros(size(r_arr));
profps = zeros(size(r_arr));
profps_err = zeros(size(r_arr));
for iter=1:nsim
    M = csvread(strcat(stackdir,'bk_ps/',dt.name,'_randprof',...
        '_',num2str(countbg),'_',num2str(iter),'.csv'));
    profcb = profcb + M(:,2)';
    profps = profps + M(:,4)';
    profcb_err = profcb_err + M(:,2)'.^2;
    profps_err = profps_err + M(:,4)'.^2;
end
profcb = profcb./nsim;
profps = profps./nsim;
profcb_err = sqrt(profcb_err./nsim - profcb.^2);
profps_err = sqrt(profps_err./nsim - profps.^2);

ihlprofdat.bkprofcb = profcb;
ihlprofdat.bkprofcb_err = profcb_err;
ihlprofdat.bkprofps = profps;
ihlprofdat.bkprofps_err = profps_err;
%%
mc = 0;
for im=10:12
%%  
    mc = mc + 1;
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    dat = mcerrdat(im).dat;
    counts = dat.counts;
    countg = dat.countg;
    ihlprofdat.data(mc).m_min = m_min;
    ihlprofdat.data(mc).m_max = m_max;
    ihlprofdat.data(mc).counts = counts;
    ihlprofdat.data(mc).countg = countg;
    
    %%% get error bar from substack %%%
    profscb_err = dat.profscb_std .* sqrt(sum(~isnan(dat.profscb_mat))./counts);
    profgcb_err = dat.profgcb_std .* sqrt(sum(~isnan(dat.profgcb_mat))./countg);
    profsps_err = dat.profsps_std .* sqrt(sum(~isnan(dat.profsps_mat))./counts);
    profgps_err = dat.profgps_std .* sqrt(sum(~isnan(dat.profgps_mat))./countg);
    
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
   %%   
   %{
    %%% plot the stacking maps %%%
    figure
    setwinsize(gcf,800,600)
    subplot(2,2,1)
    imageclip(stackscb(dx-100:dx+100,dx-100:dx+100));
    title(strcat('CIBER',{' '}, 'stack',{' '}, num2str(counts),{' '},'stars'));
    subplot(2,2,2)
    imageclip(stackgcb(dx-100:dx+100,dx-100:dx+100));
    title(strcat('CIBER',{' '}, 'stack',{' '}, num2str(countg),{' '},'gals'));
    subplot(2,2,3)
    imageclip(stacksps(dx-100:dx+100,dx-100:dx+100));
    title(strcat('PanSTARRS',{' '}, 'stack',{' '}, num2str(counts),{' '},'stars'));
    subplot(2,2,4)
    imageclip(stackgps(dx-100:dx+100,dx-100:dx+100));
    title(strcat('PanSTARRS',{' '}, 'stack',{' '}, num2str(countg),{' '},'gals'));

    suptitle(strcat('PanSTARRS',{' '},num2str(m_min),...
        ' < mAB < ',num2str(m_max)));
    
    savename=strcat(pltsavedir,dt.name,'_stackmaps',...
       num2str(m_min),'_',num2str(m_max));
    print(savename,'-dpng');close
    
    figure
    setwinsize(gcf,800,600)
    subplot(2,2,1)
    imageclip(stackscb);
    v = caxis;
    title(strcat('CIBER',{' '}, 'stack',{' '}, num2str(counts),{' '},'stars'));
    subplot(2,2,2)
    imageclip(stackgcb);
    title(strcat('CIBER',{' '}, 'stack',{' '}, num2str(countg),{' '},'gals'));
    caxis(v);
    subplot(2,2,3)
    imageclip(stacksps);
    v = caxis;
    title(strcat('PanSTARRS',{' '}, 'stack',{' '}, num2str(counts),{' '},'stars'));
    subplot(2,2,4)
    imageclip(stackgps);
    caxis(v)
    title(strcat('PanSTARRS',{' '}, 'stack',{' '}, num2str(countg),{' '},'gals'));
    suptitle(strcat('PanSTARRS',{' '},num2str(m_min),...
        ' < mAB < ',num2str(m_max)));
    savename=strcat(pltsavedir,dt.name,'_stackmaps',...
       num2str(m_min),'_',num2str(m_max),'_full');
    print(savename,'-dpng');close
    %}
    %% 
    
    %%% get the stacking profile %%%
    prof = radial_prof(stackscb,ones(2*dx+1),dx+1,dx+1,1,nbins,...
        0,'weight',mdatas);
    bksub = prof.prof - ihlprofdat.bkprofcb;
    err = sqrt(profscb_err.^2 + ihlprofdat.bkprofcb_err.^2);
    ihlprofdat.data(mc).profscb = bksub;
    ihlprofdat.data(mc).profscb_err= err;
    
    prof = radial_prof(stacksps,ones(2*dx+1),dx+1,dx+1,1,nbins,...
        0,'weight',mdatas);
    bksub = prof.prof - ihlprofdat.bkprofps;
    err = sqrt(profsps_err.^2 + ihlprofdat.bkprofps_err.^2);
    ihlprofdat.data(mc).profsps = bksub;
    ihlprofdat.data(mc).profsps_err= err;
    
    prof = radial_prof(stackgcb,ones(2*dx+1),dx+1,dx+1,1,nbins,...
        0,'weight',mdatas);
    bksub = prof.prof - ihlprofdat.bkprofcb;
    err = sqrt(profgcb_err.^2 + ihlprofdat.bkprofcb_err.^2);
    ihlprofdat.data(mc).profgcb = bksub;
    ihlprofdat.data(mc).profgcb_err= err;

    prof = radial_prof(stackgps,ones(2*dx+1),dx+1,dx+1,1,nbins,...
        0,'weight',mdatas);
    bksub = prof.prof - ihlprofdat.bkprofps;
    err = sqrt(profgps_err.^2 + ihlprofdat.bkprofps_err.^2);
    ihlprofdat.data(mc).profgps = bksub;
    ihlprofdat.data(mc).profgps_err= err;
    
end
%%  get normalized profile
mc = 0;
sp = [1];
for im=10:12
    mc = mc + 1;

    prof = ihlprofdat.data(mc).profgcb;
    prof_err = sqrt(ihlprofdat.data(mc).profgcb_err.^2);
    ihlprofdat.norm(mc).profgcb = prof;
    ihlprofdat.norm(mc).profgcb_err = prof_err;

    prof = ihlprofdat.data(mc).profscb;
    prof_err = sqrt(ihlprofdat.data(mc).profscb_err.^2);
    norm = mean(ihlprofdat.norm(mc).profgcb(sp))./mean(prof(sp));
    prof = prof.*norm;prof_err = prof_err.*norm;
    ihlprofdat.norm(mc).profscb = prof;
    ihlprofdat.norm(mc).profscb_err = prof_err;
    
    prof = ihlprofdat.data(mc).profgps;
    prof_err = sqrt(ihlprofdat.data(mc).profgps_err.^2);
    ihlprofdat.norm(mc).profgps = prof;
    ihlprofdat.norm(mc).profgps_err = prof_err;

    prof = ihlprofdat.data(mc).profsps;
    prof_err = sqrt(ihlprofdat.data(mc).profsps_err.^2);
    norm = mean(ihlprofdat.norm(mc).profgps(sp))./mean(prof(sp));
    prof = prof.*norm;prof_err = prof_err.*norm;
    ihlprofdat.norm(mc).profsps = prof;
    ihlprofdat.norm(mc).profsps_err = prof_err;
       
end
%% get excess profile
mc = 0;
for im=1:3
    mc = mc + 1;
    diffcb = ihlprofdat.norm(mc).profgcb - ihlprofdat.norm(mc).profscb;
    diffcb_err = sqrt(ihlprofdat.norm(mc).profgcb_err.^2 +...
                      ihlprofdat.norm(mc).profscb_err.^2);
    diffps = ihlprofdat.norm(im).profgps - ihlprofdat.norm(im).profsps;
    diffps_err = sqrt(ihlprofdat.norm(mc).profgps_err.^2 +...
                      ihlprofdat.norm(mc).profsps_err.^2);
    diff = diffcb - diffps;
    % there are some numer issue, so force it to 0
    diff(1) = 0;
    diff_err = sqrt(diffcb_err.^2 + diffps_err.^2);
    ihlprofdat.excess(mc).diffcb = diffcb;
    ihlprofdat.excess(mc).diffcb_err = diffcb_err;
    ihlprofdat.excess(mc).diffps = diffps;
    ihlprofdat.excess(mc).diffps_err = diffps_err;
    ihlprofdat.excess(mc).diff = diff;
    ihlprofdat.excess(mc).diff_err = diff_err;
end

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
save(sprintf('%s/%s_ihlprofdat',loaddir,dt.name),'ihlprofdat');

return