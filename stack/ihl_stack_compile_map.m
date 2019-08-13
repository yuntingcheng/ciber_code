function ihl_stack_compile_map(flight,inst,ifield,varargin)
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('ifield',@isnumeric);
  p.addOptional('masklim',false,@islogical);
  p.addOptional('rmin',nan,@isnumeric);
  p.addOptional('sample_type','jack_random',@ischar);
  p.addOptional('smooth',false,@islogical);
  p.addOptional('sm_scale',10,@isnumeric);%[pixel]
  
  p.parse(flight,inst,ifield,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  ifield   = p.Results.ifield;
  masklim = p.Results.masklim;
  rmin     = p.Results.rmin;
  sample_type=p.Results.sample_type;
  smooth   = p.Results.smooth;
  sm_scale = p.Results.sm_scale;
  
  clear p varargin;
%%
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');

psfdir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(psfdir,'fitpsfdat'),'fitpsfdat');

bkdir = (strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/bk_ps/'));

nsim = 50;
dx = 1200;
nbins = 25;

m_min_arr = 16:19;
m_max_arr = m_min_arr+1;

%%% get the PSF params %%%
bestparam = fitpsfdat(ifield).bestparam;
beta = bestparam(1);
rc = bestparam(2);
radmap = make_radius_map(zeros(2*dx+1),dx,dx).*0.7;
psfmap = (1 + (radmap/rc).^2).^(-3.*beta./2);
profile = radial_prof(psfmap,ones(2*dx+1),dx+1,dx+1,1,nbins);
profpsf_arr=(profile.prof)./profile.prof(1);
ihlprofdat.psf_arr = profpsf_arr;
r_arr=profile.r*0.7;
ihlprofdat.r_arr = r_arr;
ihlprofdat.rbinedges_arr = profile.binedges;

%%
for im=1:numel(m_min_arr)

    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    if rmin==2
        load(sprintf('%s/stackdat_%s_%d_%d_%s_rmin2',...
            loaddir,dt.name,m_min,m_max,sample_type),'stackdat');
    elseif smooth
        load(sprintf('%s/stackdat_%s_%d_%d_%s_sm%d',...
            loaddir,dt.name,m_min,m_max,sample_type,sm_scale),'stackdat');
    elseif masklim
        load(sprintf('%s/stackdat_%s_%d_%d_%s_masklim',...
            loaddir,dt.name,m_min,m_max,sample_type),'stackdat');        
    elseif isnan(rmin)
        load(sprintf('%s/stackdat_%s_%d_%d_%s',...
            loaddir,dt.name,m_min,m_max,sample_type),'stackdat');
    end
    sp100 = find(ihlprofdat.r_arr>100);
    
    ihlprofdat.data(im).m_min = m_min;
    ihlprofdat.data(im).m_max = m_max;
    ihlprofdat.data(im).counts = stackdat.all.counts;
    ihlprofdat.data(im).countg = stackdat.all.countg;
    
    ihlprofdat.data(im).profscb = stackdat.all.profcbs;
    ihlprofdat.data(im).profscb_err= stackdat.errjack.profcbs;
    ihlprofdat.data(im).profgcb = stackdat.all.profcbg;
    ihlprofdat.data(im).profgcb_err= stackdat.errjack.profcbg;
    ihlprofdat.data(im).profsps = stackdat.all.profpss;
    ihlprofdat.data(im).profsps_err= stackdat.errjack.profpss;
    ihlprofdat.data(im).profgps = stackdat.all.profpsg;
    ihlprofdat.data(im).profgps_err= stackdat.errjack.profpsg;
    ihlprofdat.data100(im).profscb = stackdat.all.profcbs100;
    ihlprofdat.data100(im).profscb_err= stackdat.errjack.profcbs100;
    ihlprofdat.data100(im).profgcb = stackdat.all.profcbg100;
    ihlprofdat.data100(im).profgcb_err = stackdat.errjack.profcbg100;
    ihlprofdat.data100(im).profsps = stackdat.all.profpss100;
    ihlprofdat.data100(im).profsps_err = stackdat.errjack.profpss100;
    ihlprofdat.data100(im).profgps = stackdat.all.profpsg100;
    ihlprofdat.data100(im).profgps_err = stackdat.errjack.profpsg100;


    %%% get the bk stacking profile %%%
    r_arr = ihlprofdat.r_arr;
    profcbs_arr = zeros([nsim,numel(r_arr)]);
    profcbg_arr = zeros([nsim,numel(r_arr)]);
    profpss_arr = zeros([nsim,numel(r_arr)]);
    profpsg_arr = zeros([nsim,numel(r_arr)]);
    profcbs100 = zeros([1,nsim]);
    profcbg100 = zeros([1,nsim]);
    profpss100 = zeros([1,nsim]);
    profpsg100 = zeros([1,nsim]);

    for iter=1:nsim
        if rmin==2
            load(sprintf('%s/stackdat_%s_rmin2_bk%d',...
                bkdir,dt.name,iter),'stackdat');
        elseif smooth
            load(sprintf('%s/stackdat_%s_sm%d_bk%d',...
                bkdir,dt.name,sm_scale,iter),'stackdat');
        elseif masklim
            load(sprintf('%s/stackdat_%s_masklim_bk%d',...
                bkdir,dt.name,iter),'stackdat');  
        elseif isnan(rmin)
            load(sprintf('%s/stackdat_%s_bk%d',...
                bkdir,dt.name,iter),'stackdat');
        end
        profcbs_arr(iter,:) = stackdat(im).profcbs_arr;
        profcbg_arr(iter,:) = stackdat(im).profcbg_arr;
        profpss_arr(iter,:) = stackdat(im).profpss_arr;
        profpsg_arr(iter,:) = stackdat(im).profpsg_arr;
        profcbs100(iter) = sum(stackdat(im).profcbs_arr(sp100).*...
            stackdat(im).hitmaps_arr(sp100))./sum(stackdat(im).hitmaps_arr(sp100));
        profcbg100(iter) = sum(stackdat(im).profcbg_arr(sp100).*...
            stackdat(im).hitmapg_arr(sp100))./sum(stackdat(im).hitmapg_arr(sp100));
        profpss100(iter) = sum(stackdat(im).profpss_arr(sp100).*...
            stackdat(im).hitmaps_arr(sp100))./sum(stackdat(im).hitmaps_arr(sp100));
        profpsg100(iter) = sum(stackdat(im).profpsg_arr(sp100).*...
            stackdat(im).hitmapg_arr(sp100))./sum(stackdat(im).hitmapg_arr(sp100));
    end
    
    profcbs_err = nanstd(profcbs_arr);
    profcbg_err = nanstd(profcbg_arr);
    profpss_err = nanstd(profpss_arr);
    profpsg_err = nanstd(profpsg_arr);
    profcbs_arr = nanmean(profcbs_arr);
    profcbg_arr = nanmean(profcbg_arr);
    profpss_arr = nanmean(profpss_arr);
    profpsg_arr = nanmean(profpsg_arr);
    profcbs_arr = spline(r_arr,profcbs_arr,r_arr);
    profcbg_arr = spline(r_arr,profcbg_arr,r_arr);
    profpss_arr = spline(r_arr,profpss_arr,r_arr);
    profpsg_arr = spline(r_arr,profpsg_arr,r_arr);
    profcbs_err = spline(r_arr,profcbs_err,r_arr);
    profcbg_err = spline(r_arr,profcbg_err,r_arr);
    profpss_err = spline(r_arr,profpss_err,r_arr);
    profpsg_err = spline(r_arr,profpsg_err,r_arr);
    
    profcbs_100err = nanstd(profcbs100);
    profcbg_100err = nanstd(profcbg100);
    profpss_100err = nanstd(profpss100);
    profpsg_100err = nanstd(profpsg100);
    profcbs_100 = nanmean(profcbs100);
    profcbg_100 = nanmean(profcbg100);
    profpss_100 = nanmean(profpss100);
    profpsg_100 = nanmean(profpsg100);

    ihlprofdat.bk(im).profscb = profcbs_arr;
    ihlprofdat.bk(im).profscb_err= profcbs_err;
    ihlprofdat.bk(im).profgcb = profcbg_arr;
    ihlprofdat.bk(im).profgcb_err= profcbg_err;
    ihlprofdat.bk(im).profsps = profpss_arr;
    ihlprofdat.bk(im).profsps_err= profpss_err;
    ihlprofdat.bk(im).profgps = profpsg_arr;
    ihlprofdat.bk(im).profgps_err= profpsg_err;
    ihlprofdat.bk100(im).profscb = profcbs_100;
    ihlprofdat.bk100(im).profscb_err= profcbs_100err;
    ihlprofdat.bk100(im).profgcb = profcbg_100;
    ihlprofdat.bk100(im).profgcb_err= profcbg_100err;
    ihlprofdat.bk100(im).profsps = profpss_100;
    ihlprofdat.bk100(im).profsps_err= profpss_100err;
    ihlprofdat.bk100(im).profgps = profpsg_100;
    ihlprofdat.bk100(im).profgps_err= profpsg_100err;
end         
%%  get normalized profile
sp = [1];
for im=1:4
    prof = ihlprofdat.data(im).profgcb - ihlprofdat.bk(im).profgcb;
    prof_err = sqrt(ihlprofdat.data(im).profgcb_err.^2 +...
        ihlprofdat.bk(im).profgcb_err.^2);
    ihlprofdat.norm(im).profgcb = prof;
    ihlprofdat.norm(im).profgcb_err = prof_err;

    prof = ihlprofdat.data(im).profscb - ihlprofdat.bk(im).profscb;
    prof_err = sqrt(ihlprofdat.data(im).profscb_err.^2 +...
        ihlprofdat.bk(im).profscb_err.^2);
    norm = mean(ihlprofdat.norm(im).profgcb(sp))./mean(prof(sp));
    prof = prof.*norm;prof_err = prof_err.*norm;
    ihlprofdat.norm(im).profscb = prof;
    ihlprofdat.norm(im).profscb_err = prof_err;
    
    prof = ihlprofdat.data(im).profgps - ihlprofdat.bk(im).profgps;
    prof_err = sqrt(ihlprofdat.data(im).profgps_err.^2 +...
        ihlprofdat.bk(im).profgps_err.^2);
    ihlprofdat.norm(im).profgps = prof;
    ihlprofdat.norm(im).profgps_err = prof_err;

    prof = ihlprofdat.data(im).profsps - ihlprofdat.bk(im).profsps;
    prof_err = sqrt(ihlprofdat.data(im).profsps_err.^2 +...
        ihlprofdat.bk(im).profsps_err.^2);
    norm = mean(ihlprofdat.norm(im).profgps(sp))./mean(prof(sp));
    prof = prof.*norm;prof_err = prof_err.*norm;
    ihlprofdat.norm(im).profsps = prof;
    ihlprofdat.norm(im).profsps_err = prof_err;
    
end
%%  get excess > 100 arcsec
sp = [1];
for im=1:4
    profg = ihlprofdat.data100(im).profgcb - ihlprofdat.bk100(im).profgcb;
    profg_err = sqrt(ihlprofdat.data100(im).profgcb_err.^2 +...
        ihlprofdat.bk100(im).profgcb_err.^2);
    profs = ihlprofdat.data100(im).profscb - ihlprofdat.bk100(im).profscb;
    profs_err = sqrt(ihlprofdat.data100(im).profscb_err.^2 +...
        ihlprofdat.bk100(im).profscb_err.^2);
    norm = mean(ihlprofdat.norm(im).profgcb(sp)) ...
        ./mean(ihlprofdat.norm(im).profscb(sp));
    profs = profs.*norm;profs_err = profs_err.*norm;
    Ecb = profg - profs;
    Ecb_err = sqrt(profg_err.^2+profs_err.^2);
    
    profg = ihlprofdat.data100(im).profgps - ihlprofdat.bk100(im).profgps;
    profg_err = sqrt(ihlprofdat.data100(im).profgps_err.^2 +...
        ihlprofdat.bk100(im).profgps_err.^2);
    profs = ihlprofdat.data100(im).profsps - ihlprofdat.bk100(im).profsps;
    profs_err = sqrt(ihlprofdat.data100(im).profsps_err.^2 +...
        ihlprofdat.bk100(im).profsps_err.^2);
    norm = mean(ihlprofdat.norm(im).profgps(sp)) ...
        ./mean(ihlprofdat.norm(im).profsps(sp));
    profs = profs.*norm;profs_err = profs_err.*norm;
    Eps = profg - profs;
    Eps_err = sqrt(profg_err.^2+profs_err.^2);
    
    E = Ecb - Eps;
    E_err = sqrt(Ecb_err.^2+Eps_err.^2);
    
    ihlprofdat.excess(im).diff100 = E;
    ihlprofdat.excess(im).diff_err100 = E_err;
    
end
%% get excess profile
for im=1:4
    diffcb = ihlprofdat.norm(im).profgcb - ihlprofdat.norm(im).profscb;
    diffcb_err = sqrt(ihlprofdat.norm(im).profgcb_err.^2 +...
                      ihlprofdat.norm(im).profscb_err.^2);
    diffps = ihlprofdat.norm(im).profgps - ihlprofdat.norm(im).profsps;
    diffps_err = sqrt(ihlprofdat.norm(im).profgps_err.^2 +...
                      ihlprofdat.norm(im).profsps_err.^2);
    diff = diffcb - diffps;
    % there are some numerical issue, so force it to 0
    diffcb(1) = 0;
    diffps(1) = 0;
    diff(1) = 0;
    diff_err = sqrt(diffcb_err.^2 + diffps_err.^2);
    ihlprofdat.excess(im).diffcb = diffcb;
    ihlprofdat.excess(im).diffcb_err = diffcb_err;
    ihlprofdat.excess(im).diffps = diffps;
    ihlprofdat.excess(im).diffps_err = diffps_err;
    ihlprofdat.excess(im).diff = diff;
    ihlprofdat.excess(im).diff_err = diff_err;
end
%%
savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
if rmin==2
    save(sprintf('%s/ihlprofdat_%s_%s_rmin2',...
        savedir,dt.name,sample_type),'ihlprofdat');
elseif smooth
    save(sprintf('%s/ihlprofdat_%s_%s_sm%d',...
        savedir,dt.name,sample_type,sm_scale),'ihlprofdat');
elseif masklim
    save(sprintf('%s/ihlprofdat_%s_%s_masklim',...
        savedir,dt.name,sample_type),'ihlprofdat');        
elseif isnan(rmin)
    save(sprintf('%s/ihlprofdat_%s_%s',...
        savedir,dt.name,sample_type),'ihlprofdat');
end

return