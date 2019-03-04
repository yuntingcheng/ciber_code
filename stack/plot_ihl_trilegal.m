flight=40030;
inst=1;
ifield=8;
quad='A';
dt=get_dark_times(flight,inst,ifield);

npix=800;
%% get the PSF
loaddir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');

load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');

bestparam=fitpsfdat(ifield).bestparam_norm;
A=bestparam(1);
B=bestparam(2);
sig=bestparam(3);
r0=bestparam(4);
alpha=bestparam(5);

radmap = make_radius_map(zeros(2*npix+1),npix+1,npix+1).*0.7;
psfmap = A*exp(-radmap.^2./2./sig^2)+B./(1+(radmap./r0).^alpha);

profile = radial_prof(psfmap,ones(2*npix+1),npix+1,npix+1,1,25);
profpsf_arr=(profile.prof)./profile.prof(1);
%%
corr=0;

figure
setwinsize(gcf,1000,800)
for m=13:20
    loaddir='/Users/ytcheng/ciber/doc/20171018_stackihl/trilegal/';
    if corr==0
        stamper=fitsread(strcat(loaddir,dt.name,'_',quad,...
            '_stamper',num2str(m),'msub.fits'));
    else
        stamper=fitsread(strcat(loaddir,dt.name,'_',quad,...
            '_stamper',num2str(m),'msub_corr.fits'));        
    end
    
    norm = stamper(npix+1,npix+1);  
    profile = radial_prof(stamper,ones(2*npix+1),npix+1,npix+1,1,25);
    r_arr=profile.r*0.7;
    prof_arr=(profile.prof)./norm;
    err_arr=profile.err./norm;
    
    subplot(2,2,1)
    errorbar(r_arr,prof_arr,err_arr,'.-','color',get_color(m));hold on
    subplot(2,2,3)
    errorbar(r_arr,prof_arr,err_arr,'.-','color',get_color(m));hold on

    loaddir='/Users/ytcheng/ciber/doc/20171018_stackihl/trilegal/';
    if corr==0
    stamper=fitsread(strcat(loaddir,dt.name,'_',quad,...
        '_stamper',num2str(m),'.fits'));
    else
    stamper=fitsread(strcat(loaddir,dt.name,'_',quad,...
        '_stamper',num2str(m),'_corr.fits'));        
    end
    
    norm = stamper(npix+1,npix+1);  
    profile = radial_prof(stamper,ones(2*npix+1),npix+1,npix+1,1,25);
    r_arr=profile.r*0.7;
    prof_arr=(profile.prof)./norm;
    err_arr=profile.err./norm;
    subplot(2,2,2)
    errorbar(r_arr,prof_arr,err_arr,'.-','color',get_color(m));hold on
    subplot(2,2,4)
    errorbar(r_arr,prof_arr,err_arr,'.-','color',get_color(m),...
        'Displayname',strcat(num2str(m),'<m<',num2str(m+1)));hold on
    
end
subplot(2,2,1)
semilogx(r_arr,profpsf_arr,'k--','linewidth',2);
set(gca, 'XScale', 'log','YScale', 'lin');
ylim([-0.1,1.1])
xlim([4e-1,7e2])
xlabel('arcsec')
ylabel('<I_{stack}>')
title('mean substracted')

subplot(2,2,3)
semilogx(r_arr,profpsf_arr,'k--','linewidth',2);
set(gca, 'XScale', 'log','YScale', 'log');
ylim([1e-3,1.2])
xlim([4e-1,7e2])
xlabel('arcsec')
ylabel('<I_{stack}>')

subplot(2,2,2)
semilogx(r_arr,profpsf_arr,'k--','linewidth',2);
set(gca, 'XScale', 'log','YScale', 'lin');
ylim([-0.1,1.1])
xlim([4e-1,7e2])
xlabel('arcsec')
ylabel('<I_{stack}>')
title('No mean substracted')

subplot(2,2,4)
semilogx(r_arr,profpsf_arr,'k--','linewidth',2,'Displayname','PSF model');
set(gca, 'XScale', 'log','YScale', 'log');
ylim([1e-3,1.2])
xlim([4e-1,7e2])
xlabel('arcsec')
ylabel('<I_{stack}>')
h=legend('show','Location','southwest');
set(h,'fontsize',10)
legend boxoff

axes; 
set(gca,'Visible','off'); 
if corr==0
    h=title('Trilegal random sources','fontsize',20);

else
    h=title('Trilegal correlated sources','fontsize',20);
end

set(h,'Visible','on');

if corr==0
    savename=strcat(loaddir,dt.name,'_rprof_trilegal');
else
    savename=strcat(loaddir,dt.name,'_rprof_trilegal_corr');
end

print(savename,'-dpng');%close
