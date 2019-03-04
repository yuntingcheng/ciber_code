%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit the spherical average PSF profile with an analytic function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
pixsize=0.7;
dx = 1200;
mypaths=get_paths(flight);

psfdir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf/TM',...
        num2str(inst),'/');

savedir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
%% plot the PSF map
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    name=dt.name;

    stamper = fitsread(strcat(psfdir,name,'_stamper.fits'));
    hitmap = fitsread(strcat(psfdir,name,'_hitmap.fits'));
    psfmap = stamper./hitmap;
    
    figure
    imageclip(psfmap./psfmap(dx+1,dx+1));
    
    title(dt.name,'fontsize',18)
    savename=strcat(savedir,'psfmap_i',num2str(ifield));
    print(savename,'-dpng');%close    

end
%% fit psf 1D profile
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    name=dt.name;

    stamper = fitsread(strcat(psfdir,name,'_stamper.fits'));
    hitmap = fitsread(strcat(psfdir,name,'_hitmap.fits'));
    psfmap = stamper./hitmap;
    radmap = make_radius_map(psfmap,dx+1,dx+1) .* pixsize;
    profile = radial_prof(psfmap,ones(2*dx+1),dx+1,dx+1,1,32);
    r_arr=profile.r*pixsize;
    p_arr=profile.prof./profile.prof(1);
    e_arr=profile.err./profile.prof(1);

    C = mean(psfmap(find(radmap > 200))) / profile.prof(1);
    beta_arr = 2:0.1:4;
    rc_arr = 5:0.2:12;
    err2best = 1e8;
    chi2best = 1e8;
    for beta = beta_arr
    for rc = rc_arr
        sp = r_arr<15;
        m_arr = (1 + (r_arr/rc).^2).^(-3.*beta./2);
        chi2=sum((m_arr(sp)-p_arr(sp)).^2./e_arr(sp).^2)./(numel(r_arr)-2);
        err2=sum((m_arr-p_arr).^2);

        if chi2<chi2best
            err2best = err2;
            chi2best=chi2;
            bestparam=[beta,rc,C,chi2];
            mb_arr=m_arr;
        end
    end
    end
    
    fitpsfdat(ifield).r_arr=r_arr;
    fitpsfdat(ifield).p_arr=p_arr;
    fitpsfdat(ifield).e_arr=e_arr;
    fitpsfdat(ifield).bestparam=bestparam;%[beta,rc,C,chi2best]
    fitpsfdat(ifield).m_arr=mb_arr;% best model
end
save(strcat(savedir,'fitpsfdat'),'fitpsfdat');
%% plot the fitting results
load(strcat(savedir,'fitpsfdat'),'fitpsfdat');
for ifield=4:8
    figure
    dt=get_dark_times(flight,inst,ifield);
    r_arr=fitpsfdat(ifield).r_arr;
    p_arr=fitpsfdat(ifield).p_arr;
    e_arr=fitpsfdat(ifield).e_arr;
    m_arr=fitpsfdat(ifield).m_arr;
    
    C = fitpsfdat(ifield).bestparam(3);
    chi2 = fitpsfdat(ifield).bestparam(4);
    loglog(r_arr,m_arr,'r--','linewidth',2);hold on
    loglog(r_arr,ones(size(r_arr))*C,'m--','linewidth',2);hold on
    loglog(r_arr,m_arr + C,'b','linewidth',2);hold on
    errorbar(r_arr,p_arr,e_arr,'ok','linewidth',1.5,'DisplayName','Data');
    
    title(dt.name,'fontsize',18)
    xlim([0.5,5e2])
    ylim([1e-5,1.2])
    xlabel('r [arcsec]','fontsize',18);
    ylabel('Normalized Stacking Profile','fontsize',18);
    h=legend({'\beta-model fit','constant background',...
        '\beta-model + background','Stacking data'},'Location','northeast');
    set(h,'fontsize',12)
    legend boxoff
    
    savename=strcat(savedir,'psf_fit_i',num2str(ifield));
    print(savename,'-dpng');%close    

end
%% plot all fields in one fig
load(strcat(savedir,'fitpsfdat'),'fitpsfdat');
figure
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    r_arr=fitpsfdat(ifield).r_arr;
    p_arr=fitpsfdat(ifield).p_arr;

    loglog(r_arr,p_arr,'color',get_color(ifield-3),...
        'linewidth',2,'DisplayName',dt.name);hold on
    xlim([0.5,5e2])
    ylim([1e-5,1.2])
    xlabel('r(arcsec)','fontsize',18);
    ylabel('PSF','fontsize',18);

end
h=legend('show','Location','northeast');
set(h,'fontsize',12)
legend boxoff

savename=strcat(savedir,'psf_fit');
print(savename,'-dpng');%close  

%% plot I(<r) 
figure
load(strcat(savedir,'fitpsfdat'),'fitpsfdat');
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);

    beta = fitpsfdat(ifield).bestparam(1);
    rc = fitpsfdat(ifield).bestparam(2);
    C = fitpsfdat(ifield).bestparam(3);
    npix=2000;
    radmap = make_radius_map(zeros(2*npix+1),npix+1,npix+1).*pixsize;
    psfmap = (1 + (radmap./rc).^2).^(-3.*beta./2);
    redge_arr=logspace(log10(0.7),log10(max(radmap(:))),30);
    r_arr=(redge_arr(2:end)+redge_arr(1:end-1))/2;
    cprof_arr=zeros(size(r_arr));
    for ir=1:numel(r_arr)
        sp=find(radmap<=redge_arr(ir+1));
        cprof_arr(ir)=sum(psfmap(sp));
    end
    %%% plot to see when it converge
    semilogx(r_arr,cprof_arr./cprof_arr(end),'-',...
        'linewidth',2,'DisplayName',dt.name);hold on
    drawnow

end

for ipix=1:2:5
plot([ipix/2*7 ipix/2*7],[0 2],'k','Displayname',sprintf('%dx%d',ipix,ipix))
end

xlim([5e-1,1e2])
ylim([0,1.1])
h=legend('show','Location','southeast');
set(h,'fontsize',10)
legend boxoff
xlabel('$r(arcsec)$','interpreter','latex','fontsize',18)
ylabel('$I(<r)$','interpreter','latex','fontsize',18)
savename=strcat(savedir,'psf_int');
print(savename,'-dpng');%close  

%% normalize
load(strcat(savedir,'fitpsfdat'),'fitpsfdat');
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    beta = fitpsfdat(ifield).bestparam(1);
    rc = fitpsfdat(ifield).bestparam(2);
    C = fitpsfdat(ifield).bestparam(3);
    npix=2000;
    radmap = make_radius_map(zeros(2*npix+1),npix+1,npix+1).*pixsize;
    psfmap = (1 + (radmap./rc).^2).^(-3.*beta./2);
    norm=sum(psfmap(:));    
    fitpsfdat(ifield).norm=norm;

end
save(strcat(savedir,'fitpsfdat'),'fitpsfdat');
%% fraction of PSF inside a pixel 
load(strcat(savedir,'fitpsfdat'),'fitpsfdat');
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    beta = fitpsfdat(ifield).bestparam(1);
    rc = fitpsfdat(ifield).bestparam(2);
    C = fitpsfdat(ifield).bestparam(3);
    norm = fitpsfdat(ifield).norm;
    npix=4.5;
    radmap = make_radius_map(zeros(2*npix+1),npix+1,npix+1).*pixsize;
    psfmap =  (1./norm).*(1 + (radmap./rc).^2).^(-3.*beta./2);
    disp(sprintf('field%d, psf_frac=%.3f',ifield,sum(psfmap(:))));
    fitpsfdat(ifield).pix_frac=sum(psfmap(:));
end
save(strcat(savedir,'fitpsfdat'),'fitpsfdat');