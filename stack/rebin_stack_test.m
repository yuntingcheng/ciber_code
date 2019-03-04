flight=40030;
inst=1;
ifield=8;
npix=800;
Nsamp=100;
savedir=(strcat('/Users/ytcheng/ciber/doc/20171130_psfstack/plots/'));

loaddir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');

load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');

bestparam=fitpsfdat(ifield).bestparam_norm;
A=bestparam(1);
B=bestparam(2);
sig=bestparam(3);
r0=bestparam(4);
alpha=bestparam(5);

scale=10;
sigs=sig.*scale;
r0s=r0.*scale;
radmap = make_radius_map(zeros(1601),801,801).*0.7;
psfmap = A*exp(-radmap.^2./2./sigs^2)+B./(1+(radmap./r0s).^alpha);
As=A./sum(psfmap(:));
Bs=B./sum(psfmap(:));
%% sim stamp
%%% stamp with rebin
stamp=zeros(2*npix+1,2*npix+1);
src_coord=npix+0.5+10*rand(2,Nsamp);

for i=1:Nsamp
    
    xsrc=src_coord(1,i);
    ysrc=src_coord(2,i);
    radmap = make_radius_map(zeros(2*npix+10),xsrc,ysrc).*0.7;
    psfmap = As*exp(-radmap.^2./2./sigs^2)+Bs./(1+(radmap./r0s).^alpha);

    psfmap_coarse=rebin_map_coarse(psfmap,10);
    psfmap_fine=imresize(psfmap_coarse,10,'method','nearest');

    stamp=stamp+psfmap_fine(round(xsrc)-npix:round(xsrc)+npix,...
                            round(ysrc)-npix:round(ysrc)+npix);
end
stamp=stamp./Nsamp;
%%
%%% stamp without rebin
stamp_exact=zeros(2*npix+1,2*npix+1);

for i=1:Nsamp
    
    xsrc=src_coord(1,i);
    ysrc=src_coord(2,i);
    radmap = make_radius_map(zeros(2*npix+10),xsrc,ysrc).*0.7;
    psfmap = As*exp(-radmap.^2./2./sigs^2)+Bs./(1+(radmap./r0s).^alpha);

    stamp_exact=stamp_exact+psfmap(round(xsrc)-npix:round(xsrc)+npix,...
                            round(ysrc)-npix:round(ysrc)+npix);
end
stamp_exact=stamp_exact./Nsamp;
%%%
radmap = make_radius_map(zeros(2*npix+1),npix+1,npix+1).*0.7;
psfmap = As*exp(-radmap.^2./2./sigs^2)+Bs./(1+(radmap./r0s).^alpha);
%%

v=[-1e-4/(scale^2),1e-4/(scale^2)];

figure
setwinsize(gcf,1000,400)
subplot(1,2,1)
imageclip(stamp(npix+1-50*scale:npix+1+50*scale,npix+1-50*scale:npix+1+50*scale)-...
    psfmap(npix+1-50*scale:npix+1+50*scale,npix+1-50*scale:npix+1+50*scale));
title('rebinned stacking-input PSF');
caxis(v);
v=caxis;
subplot(1,2,2)
imageclip(stamp_exact(npix+1-50*scale:npix+1+50*scale,npix+1-50*scale:npix+1+50*scale)-...
    psfmap(npix+1-50*scale:npix+1+50*scale,npix+1-50*scale:npix+1+50*scale));
caxis(v);
title('ideal stacking-input PSF');

savename=strcat(savedir,'stack_test2D_diff',num2str(scale));
print(savename,'-dpng');
%%
profile = radial_prof(psfmap,ones(2*npix+1),npix+1,npix+1,1,25);
profpsf_arr=(profile.prof)./profile.prof(1);

profile = radial_prof(stamp,ones(2*npix+1),npix+1,npix+1,1,25);
profstamp_arr=(profile.prof)./profile.prof(1);

profile = radial_prof(stamp_exact,ones(2*npix+1),npix+1,npix+1,1,25);
profstampex_arr=(profile.prof)./profile.prof(1);

r_arr=profile.r*0.7;

figure
setwinsize(gcf,1000,400)
subplot(1,2,1)
semilogx(r_arr,profpsf_arr,'k--','linewidth',2,'Displayname','input PSF');hold on
semilogx(r_arr,profstamp_arr,'r','Displayname','rebinned stacking');
semilogx(r_arr,profstampex_arr,'b','Displayname','ideal stacking');
ylim([-0.1,1.1])
xlim([4e-1,7e2])
xlabel('arcsec')
ylabel('<I_{stack}>')

subplot(1,2,2)
loglog(r_arr,profpsf_arr,'k--','linewidth',2,'Displayname','input PSF');hold on
loglog(r_arr,profstamp_arr,'r','Displayname','rebinned stacking');
loglog(r_arr,profstampex_arr,'b','Displayname','ideal stacking');
xlim([4e-1,7e2])
ylim([1e-4,2e0])
xlabel('arcsec')
ylabel('<I_{stack}>')
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff

%savename=strcat(savedir,'stack_test2D_scale',num2str(scale));
%print(savename,'-dpng');