%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit the spherical average PSF profile with an analytic function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=2;
pixsize=0.7;
mypaths=get_paths(flight);

psfdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf/yt/inst',...
        num2str(inst),'/j0_14/');

savedir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
%% fit psf 1D profile
load(strcat(savedir,'fitpsfdat'),'fitpsfdat');
for ifield=4:8
dt=get_dark_times(flight,inst,ifield);
name=dt.name;

a=fitsread(strcat(psfdir,name,'_A.fits'));
b=fitsread(strcat(psfdir,name,'_B.fits'));
c=fitsread(strcat(psfdir,name,'_C.fits'));
d=fitsread(strcat(psfdir,name,'_D.fits'));
ah=fitsread(strcat(psfdir,name,'_Ahitmap.fits'));
bh=fitsread(strcat(psfdir,name,'_Bhitmap.fits'));
ch=fitsread(strcat(psfdir,name,'_Chitmap.fits'));
dh=fitsread(strcat(psfdir,name,'_Dhitmap.fits'));

psfmap=(a.*ah+b.*bh+c.*ch+d.*dh)./(ah+bh+ch+dh);

profile = radial_prof_pk(psfmap,ones(801),401,401,300);
r_arr=profile.r*pixsize;
p_arr=profile.prof./profile.prof(1);
e_arr=profile.err./profile.prof(1);   

A_arr=0.1:0.1:3;
B_arr=0.1:0.1:3;
sig_arr=3:0.1:6;
r0_arr=3:0.1:5;
alpha_arr=2.5:0.1:5;

chi2best=1e8;
bestparam=zeros(1,5);
for A=A_arr
for B=B_arr
for sig=sig_arr
for r0=r0_arr
for alpha=alpha_arr
    m_arr=A*exp(-r_arr.^2./2./sig^2)+B./(1+(r_arr./r0).^alpha);
    sp=find(r_arr<20);
    chi2=sum((m_arr(sp)-p_arr(sp)).^2./e_arr(sp).^2)./(numel(sp)-5);
    if chi2<chi2best
        chi2best=chi2;
        bestparam=[A,B,sig,r0,alpha,chi2best];
        mb_arr=m_arr;
    end
    
end
end
end
end
end
fitpsfdat(ifield).r_arr=r_arr;
fitpsfdat(ifield).p_arr=p_arr;
fitpsfdat(ifield).e_arr=e_arr;
fitpsfdat(ifield).bestparam=bestparam;%[A,B,sig,r0,alpha,chi2best]
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
    
    bestparam=fitpsfdat(ifield).bestparam;
    A=bestparam(1);
    B=bestparam(2);
    sig=bestparam(3);
    r0=bestparam(4);
    alpha=bestparam(5);
    chi2=bestparam(6);
    
    m1_arr=A*exp(-r_arr.^2./2./sig^2);
    m2_arr=B./(1+(r_arr./r0).^alpha);
    loglog(r_arr,m_arr,'r','DisplayName','Model Fit');hold on
    loglog(r_arr,m1_arr,'b--','DisplayName','exponential')
    loglog(r_arr,m2_arr,'color',[0,0.8,0.2],'linestyle','--',...
        'DisplayName','power law')
    errorbar(r_arr,p_arr,e_arr,'o-k','DisplayName','Data');
    
    title(dt.name,'fontsize',18)
    xlim([0.5,5e2])
    ylim([1e-5,1.2])
    xlabel('r(arcsec)','fontsize',18);
    ylabel('PSF','fontsize',18);
    h=legend('show','Location','northeast');
    set(h,'fontsize',12)
    legend boxoff

    strchi2 = ['Chi2/DoF (r<20)= ',num2str(chi2,'%.2f')];
    text(50,0.06,strchi2,'fontsize',12);
    
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
    e_arr=fitpsfdat(ifield).e_arr;
    m_arr=fitpsfdat(ifield).m_arr;
    
    bestparam=fitpsfdat(ifield).bestparam;
    A=bestparam(1);
    B=bestparam(2);
    sig=bestparam(3);
    r0=bestparam(4);
    alpha=bestparam(5);
    chi2=bestparam(6);
    
    m1_arr=A*exp(-r_arr.^2./2./sig^2);
    m2_arr=B./(1+(r_arr./r0).^alpha);
    
    loglog(r_arr,p_arr,'color',get_color(ifield-3),'DisplayName',dt.name);hold on
    
end
title(strcat('TM',num2str(inst)),'fontsize',18)
xlim([0.5,5e2])
ylim([1e-5,1.2])
xlabel('r(arcsec)','fontsize',18);
ylabel('PSF','fontsize',18);
h=legend('show','Location','northeast');
set(h,'fontsize',12)
legend boxoff

savename=strcat(savedir,'psf_fit');
print(savename,'-dpng');%close    


%% normalize
load(strcat(savedir,'fitpsfdat'),'fitpsfdat');
for ifield=4:8
    bestparam=fitpsfdat(ifield).bestparam;
    A=bestparam(1);
    B=bestparam(2);
    sig=bestparam(3);
    r0=bestparam(4);
    alpha=bestparam(5);
    chi2best=bestparam(6);
    npix=2000;
    radmap = make_radius_map(zeros(2*npix+1),npix,npix).*pixsize;
    psfmap = A*exp(-radmap.^2./2./sig^2)+B./(1+(radmap./r0).^alpha);
    redge_arr=logspace(log10(200),log10(max(radmap(:))),30);
    r_arr=(redge_arr(2:end)+redge_arr(1:end-1))/2;
    cprof_arr=zeros(size(r_arr));
    for ir=1:numel(r_arr)
        sp=find(radmap<=redge_arr(ir+1));
        cprof_arr(ir)=sum(psfmap(sp));
    end
    %%% plot to see when it converge
    figure
    semilogx(r_arr,cprof_arr,'-o');
    drawnow
    
    norm=cprof_arr(end);
    bestparam_norm=[A/norm,B/norm,sig,r0,alpha,chi2best];
    fitpsfdat(ifield).bestparam_norm=bestparam_norm;

end
save(strcat(savedir,'fitpsfdat'),'fitpsfdat');

%% plot I(<r) 
figure
load(strcat(savedir,'fitpsfdat'),'fitpsfdat');
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);

    bestparam=fitpsfdat(ifield).bestparam;
    A=bestparam(1);
    B=bestparam(2);
    sig=bestparam(3);
    r0=bestparam(4);
    alpha=bestparam(5);
    chi2best=bestparam(6);
    npix=2000;
    radmap = make_radius_map(zeros(2*npix+1),npix,npix).*pixsize;
    psfmap = A*exp(-radmap.^2./2./sig^2)+B./(1+(radmap./r0).^alpha);
    redge_arr=logspace(log10(0.7),log10(max(radmap(:))),30);
    r_arr=(redge_arr(2:end)+redge_arr(1:end-1))/2;
    cprof_arr=zeros(size(r_arr));
    for ir=1:numel(r_arr)
        sp=find(radmap<=redge_arr(ir+1));
        cprof_arr(ir)=sum(psfmap(sp));
    end
    %%% plot to see when it converge
    semilogx(r_arr,cprof_arr./cprof_arr(end),'-','DisplayName',dt.name);hold on
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
%% write the best fit param to txt file for IDL read
load(strcat(savedir,'fitpsfdat'),'fitpsfdat');
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    
    bestparam=fitpsfdat(ifield).bestparam_norm;    
    dlmwrite(strcat(savedir,dt.name,'_bestparam.txt'),bestparam);
end
%% fraction of PSF inside a pixel 
load(strcat(savedir,'fitpsfdat'),'fitpsfdat');
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);

    bestparam=fitpsfdat(ifield).bestparam_norm;
    A=bestparam(1);
    B=bestparam(2);
    sig=bestparam(3);
    r0=bestparam(4);
    alpha=bestparam(5);
    chi2best=bestparam(6);
    npix=4.5;
    radmap = make_radius_map(zeros(2*npix+1),npix+1,npix+1).*pixsize;
    psfmap = A*exp(-radmap.^2./2./sig^2)+B./(1+(radmap./r0).^alpha);
    disp(sprintf('field%d, psf_frac=%.3f',ifield,sum(psfmap(:))));
    fitpsfdat(ifield).pix_frac=sum(psfmap(:));
end
save(strcat(savedir,'fitpsfdat'),'fitpsfdat');
