function ihl_stack_compile_sim(flight,inst,ifield, set)
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
load(sprintf('%sstackmapdatsim',loaddir),'stackmapdatsim');
load(sprintf('%s/%s_mcerrdatsim%d',loaddir,dt.name,set),'mcerrdat');
load(sprintf('%s/stackcountdatsim_%s_%d',loaddir,dt.name,set),'stackcountdat');

dx = 1200;
nbins = 25;

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
for im = 1:3        
    counts_arr(im) = stackcountdat(im).counts;
    countg_arr(im) = stackcountdat(im).countg;
end

%%
mc = 0;
for im=1:3
%%  
    mc = mc + 1;
    m_min = im + 15;
    m_max = m_min + 1;
    counts = counts_arr(im);
    countg = countg_arr(im);
    ihlprofdat.data(mc).m_min = m_min;
    ihlprofdat.data(mc).m_max = m_max;
    ihlprofdat.data(mc).counts = counts;
    ihlprofdat.data(mc).countg = countg;
    
    %%% get error bar from substack %%%
    dat = mcerrdat(im).dat;
    profscb_err = dat.profscb_std .* sqrt(sum(~isnan(dat.profscb_mat))./counts);
    profgcb_err = dat.profgcb_std .* sqrt(sum(~isnan(dat.profgcb_mat))./countg);
    profsps_err = dat.profsps_std .* sqrt(sum(~isnan(dat.profsps_mat))./counts);
    profgps_err = dat.profgps_std .* sqrt(sum(~isnan(dat.profgps_mat))./countg);

    %%% get the stacking maps %%%
    datascb=fitsread(strcat(stackdir,'/ciber_sim/',dt.name,...
        '_stampers',num2str(m_min),'_',num2str(m_max),'_',num2str(set),'.fits'));
    datagcb=fitsread(strcat(stackdir,'/ciber_sim/',dt.name,...
        '_stamperg',num2str(m_min),'_',num2str(m_max),'_',num2str(set),'.fits'));                    
    datasps=fitsread(strcat(stackdir,'/panstarrs_sim/',dt.name,...
        '_stampers',num2str(m_min),'_',num2str(m_max),'_',num2str(set),'.fits'));
    datagps=fitsread(strcat(stackdir,'/panstarrs_sim/',dt.name,...
        '_stamperg',num2str(m_min),'_',num2str(m_max),'_',num2str(set),'.fits'));                    
    mdatas=load(strcat(stackdir,'/ciber_sim/',dt.name,...
        '_hitmaps',num2str(m_min),'_',num2str(m_max),'_',num2str(set),'.mat'));
    mdatas = mdatas.hitmap;
    mdatag=load(strcat(stackdir,'/ciber_sim/',dt.name,...
        '_hitmapg',num2str(m_min),'_',num2str(m_max),'_',num2str(set),'.mat'));
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
    %print(savename,'-dpng');close
    drawnow;close 
    
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
    %print(savename,'-dpng');close
    drawnow;close
    %}
    %%
    %%% get the stacking profile %%%
    prof = radial_prof(stackscb,ones(2*dx+1),dx+1,dx+1,1,nbins,...
        0,'weight',mdatas);
    ihlprofdat.data(mc).profscb = prof.prof;
    ihlprofdat.data(mc).profscb_err= profscb_err;
    
    prof = radial_prof(stacksps,ones(2*dx+1),dx+1,dx+1,1,nbins,...
        0,'weight',mdatas);
    ihlprofdat.data(mc).profsps = prof.prof;
    ihlprofdat.data(mc).profsps_err= profsps_err;
    
    prof = radial_prof(stackgcb,ones(2*dx+1),dx+1,dx+1,1,nbins,...
        0,'weight',mdatas);
    ihlprofdat.data(mc).profgcb = prof.prof;
    ihlprofdat.data(mc).profgcb_err= profgcb_err;

    prof = radial_prof(stackgps,ones(2*dx+1),dx+1,dx+1,1,nbins,...
        0,'weight',mdatas);
    ihlprofdat.data(mc).profgps = prof.prof;
    ihlprofdat.data(mc).profgps_err= profgps_err;
    
end
%%  get normalized profile
mc = 0;
sp = [1];
for im=1:3
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
    diff_err = sqrt(diffcb_err.^2 + diffps_err.^2);
    ihlprofdat.excess(mc).diffcb = diffcb;
    ihlprofdat.excess(mc).diffcb_err = diffcb_err;
    ihlprofdat.excess(mc).diffps = diffps;
    ihlprofdat.excess(mc).diffps_err = diffps_err;
    ihlprofdat.excess(mc).diff = diff;
    ihlprofdat.excess(mc).diff_err = diff_err;
end

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
save(sprintf('%s/%s_ihlprofdatsim%d',loaddir,dt.name,set),'ihlprofdat');

return