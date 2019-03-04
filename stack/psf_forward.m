flight=40030;
mypaths=get_paths(flight);
inst = 2;
Nsims=100;% # of random position in the sims
pixsize = 0.7;
dx=1200;

savedir=(strcat(mypaths.ciberdir,'doc/20171130_psfstack/forward_model/'));
loaddir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');
psfdir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf/TM',...
        num2str(inst),'/');

for ifield=8:-1:4
%%
dt = get_dark_times(flight,inst,ifield);

%%% get the CIBER stacked PSF map and profile %%%
stamper = fitsread(strcat(psfdir,dt.name,'_stamper.fits'));
hitmap = fitsread(strcat(psfdir,dt.name,'_hitmap.fits'));
mapcb = stamper./hitmap;
radmap = make_radius_map(mapcb,dx+1,dx+1) .* pixsize;

r_arr = fitpsfdat(ifield).r_arr;
pcb_arr = fitpsfdat(ifield).p_arr;
ecb_arr = fitpsfdat(ifield).e_arr;
mcb_arr = fitpsfdat(ifield).m_arr;
beta_cb = fitpsfdat(ifield).bestparam(1);
rc_cb = fitpsfdat(ifield).bestparam(2);
C = fitpsfdat(ifield).bestparam(3);


mapcb = mapcb - C;
%%
%%% iterate on beta, C %%%%
beta_arr = 0.5:0.5:3;
rc_arr = 1:0.5:6;

beta_best = 0.;
rc_best = 0.;
loss_best = 1e10;
stamp_best=zeros(2*dx+1,2*dx+1);
psim_best_arr = pcb_arr;
psfmap_best = mapcb;

% src position random in a pixel between [npix+0.5,npix+0.5+10]
src_coord=dx+0.5+10*rand(2,Nsims);
%%
for beta = beta_arr
    disp(sprintf('field%d, beta = %.1f',ifield,beta));
for rc = rc_arr
    %%% do the random center sim for Nsims times %%%
    stamp=zeros(2*dx+1,2*dx+1);
    for isim = 1:Nsims
        xsrc=src_coord(1,isim);
        ysrc=src_coord(2,isim);
        radmap = make_radius_map(zeros(2*dx+10),xsrc,ysrc).*0.7;
        psfmap = (1 + (radmap./rc).^2).^(-3.*beta./2);    
        psfmap_coarse=rebin_map_coarse(psfmap,10);
        psfmap_fine=imresize(psfmap_coarse,10,'method','nearest');
        stamp=stamp+psfmap_fine(round(xsrc)-dx:round(xsrc)+dx,...
                                round(ysrc)-dx:round(ysrc)+dx);
    end
    stamp=stamp./Nsims;
    
    %%% calculate the profile %%%
    profile = radial_prof(stamp,ones(2*dx+1),dx+1,dx+1,1,32);
    psim_arr=(profile.prof)./profile.prof(1);
    r_arr = profile.r.*pixsize;
    m_arr = (1 + (r_arr/rc).^2).^(-3.*beta./2);
    m_arr = m_arr./m_arr(1);
    %%% calculate the chi2 %%%%
    sp = find(r_arr<15);
    %loss = (log10(pcb_arr) - log10(psim_arr)).^2 ./ (ecb_arr.^2);
    loss = (pcb_arr - psim_arr).^2 ./ (ecb_arr.^2);
    %loss = (pcb_arr - psim_arr).^2;
    loss = sum(loss(sp));
    if loss < loss_best
        loss_best = loss;
        beta_best = beta;
        rc_best = rc;
        stamp_best = stamp;
        psim_best_arr = psim_arr;
        m_best_arr = m_arr;
    end
end
end

disp(sprintf('TM%d, field%d, beta = %.1f, rc = %.1f',...
    inst,ifield,beta_best,rc_best));

%%% save the best fit params
psfmodel.beta_best = beta_best;
psfmodel.rc_best = rc_best;
psfmodel.r_arr = r_arr;
psfmodel.prof_best_arr = psim_best_arr;

fitpsfdat(ifield).psfmodel = psfmodel;
save(strcat(loaddir,'fitpsfdat'),'fitpsfdat');

%%
%%% plot the results
fig=figure;
setwinsize(gcf,800,400)

subplot(1,2,1)
errorbar(r_arr,pcb_arr,ecb_arr,'ok','linewidth',1.5,...
    'Displayname','CIBER data stack');hold on
plot(r_arr,psim_best_arr + C,'b','linewidth',2,...
    'Displayname','best fit \beta-model stack');
plot(r_arr,m_best_arr + C,'r','linewidth',2,...
    'Displayname','best fit \beta-model PSF');

set(gca,'XScale','log','YScale','lin');
ylim([-0.1,1.1])
xlim([4e-1,7e2])
xlabel('arcsec', 'fontsize', 18)
ylabel('Normalized Profile', 'fontsize', 18)
h=legend('show','Location','northeast');
set(h,'fontsize',12)
legend boxoff

subplot(1,2,2)
errorbar(r_arr,pcb_arr,ecb_arr,'ok','linewidth',1.5,...
    'Displayname','CIBER data stack');hold on
plot(r_arr,psim_best_arr + C,'b','linewidth',2,...
    'Displayname','best fit model stack');
plot(r_arr,m_best_arr + C,'r','linewidth',2,...
    'Displayname','best fit model stack');
set(gca,'XScale','log','YScale','log');
ylim([1e-5,1.2])
xlim([0.5,1e3])
xlabel('r [arcsec]', 'fontsize', 18)
ylabel('Normalized Profile', 'fontsize', 18)

drawnow
savename = strcat(savedir,'TM',num2str(inst),'_',dt.name,'_psfmodel');
print(savename,'-dpng');

end

%% calculate normalization
flight=40030;
inst = 2;
mypaths=get_paths(flight);
pixsize = 0.7;
dx=1200;

loaddir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');

for ifield=4:8
    beta=fitpsfdat(ifield).psfmodel.beta_best;
    rc=fitpsfdat(ifield).psfmodel.rc_best;
    
    % normalize A,B
    radmap = make_radius_map(zeros(2*dx+1),dx+1,dx+1).*0.7;
    psfmap =(1 + (radmap/rc).^2).^(-3.*beta./2);
    norm=1/sum(psfmap(:));
    
    fitpsfdat(ifield).psfmodel.norm = norm;
end
save(strcat(loaddir,'fitpsfdat'),'fitpsfdat');
