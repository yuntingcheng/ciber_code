flight=40030;
inst=1;
ifield=8;
dt=get_dark_times(flight,inst,ifield);
npix=400;

savedir=(strcat('/Users/ytcheng/ciber/doc/20171130_psfstack/plots/'));

loaddir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
%%

load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');

bestparam=fitpsfdat(ifield).bestparam_norm;
A=bestparam(1);
B=bestparam(2);
sig=bestparam(3);
r0=bestparam(4);
alpha=bestparam(5);

radmap = make_radius_map(zeros(2*npix+1),npix,npix).*0.7;
psfmap = A*exp(-radmap.^2./2./sig^2)+B./(1+(radmap./r0).^alpha);
profile = radial_prof(psfmap,ones(2*npix+1),npix+1,npix+1,1,15);
psffit_arr=(profile.prof)./profile.prof(1);
%% get the CIBER stacked PSF map
psfdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf/yt/inst',...
        num2str(inst),'/j0_14/');

quad_arr=['A','B','C','D'];
radmap = make_radius_map(zeros(801),401,401).*0.7;
qtot=zeros(801);
for iquad=1:numel(quad_arr)
quad=quad_arr(iquad);
q=fitsread(strcat(psfdir,dt.name,'_',quad,'.fits'));
qh=fitsread(strcat(psfdir,dt.name,'_',quad,'hitmap.fits'));
qmap=q./qh;
bk=mean(qmap(find(radmap>100)));
qmap=qmap-bk;

qtot=qtot+qmap;
end
qtot=qtot./sum(qtot(:));

profile = radial_prof(qtot,ones(size(qmap)),401,401,1,15);
psf_arr=(profile.prof)./profile.prof(1);
epsf_arr=profile.err./profile.prof(1);
%% sim stamp Gaussian
%{
Nsamp=100;
stamp=zeros(2*npix+1,2*npix+1);
src_coord=npix+0.5+10*rand(2,Nsamp);

scale=0.7;
radmap = make_radius_map(zeros(801),401,401).*0.7;
psfmap = exp(-radmap.^2./2./(sig.*scale)^2);
Ag=1/sum(psfmap(:));

for i=1:Nsamp
    
    xsrc=src_coord(1,i);
    ysrc=src_coord(2,i);
    radmap = make_radius_map(zeros(2*npix+10),xsrc,ysrc).*0.7;
    psfmap = Ag*exp(-radmap.^2./2./(sig.*scale)^2);
    psfmap_coarse=rebin_map_coarse(psfmap,10);
    psfmap_fine=imresize(psfmap_coarse,10,'method','nearest');

    stamp=stamp+psfmap_fine(round(xsrc)-npix:round(xsrc)+npix,...
                            round(ysrc)-npix:round(ysrc)+npix);
end
stamp=stamp./Nsamp;

savename=strcat(savedir,'fit_gauss');
%}
%% sim stamp Gaussian+power law

Nsamp=100;
% uniform between [npix+0.5,npix+0.5+10]
src_coord=npix+0.5+10*rand(2,Nsamp);
stamp=zeros(2*npix+1,2*npix+1);

sA=0.55;
sB=0.7;
alphas=alpha*0.9;
radmap = make_radius_map(zeros(801),401,401).*0.7;
psfmap = A*exp(-radmap.^2./2./(sig*sA)^2)+...
    B./(1+(radmap./(r0*sB)).^alphas);
As=A/sum(psfmap(:));
Bs=B/sum(psfmap(:));

for i=1:Nsamp
    xsrc=src_coord(1,i);
    ysrc=src_coord(2,i);
    radmap = make_radius_map(zeros(2*npix+10),xsrc,ysrc).*0.7;
    psfmap = As*exp(-radmap.^2./2./(sig*sA)^2)+...
             Bs./(1+(radmap./(r0*sB)).^alphas);
    psfmap_coarse=rebin_map_coarse(psfmap,10);
    psfmap_fine=imresize(psfmap_coarse,10,'method','nearest');

    stamp=stamp+psfmap_fine(round(xsrc)-npix:round(xsrc)+npix,...
                            round(ysrc)-npix:round(ysrc)+npix);
end
savename=strcat(savedir,'fit_gausspw');
stamp=stamp./Nsamp;
%%
profile = radial_prof(stamp,ones(2*npix+1),npix+1,npix+1,1,15);
profstamp_arr=(profile.prof)./profile.prof(1);


r_arr=profile.r*0.7;

fig=figure;
setwinsize(gcf,1300,400)

v=[-1e-4,1e-4];
subplot(1,3,1)
x_arr=0.7*(-50:50);
imagesc(x_arr,x_arr,stamp(npix+1-50:npix+1+50,npix+1-50:npix+1+50)-...
    qtot(npix+1-50:npix+1+50,npix+1-50:npix+1+50));
title('rebinned stack-CIBER stack');
xlabel('arcsec')
ylabel('arcsec')
caxis(v);
v=caxis;

subplot(1,3,2)
errorbar(r_arr,psf_arr,epsf_arr,'ro-','Displayname','CIBER stack');hold on
semilogx(r_arr,profstamp_arr,'color',[0.00  0.50  0.00],'linewidth',2,...
    'Displayname','trial Gaussian model');
semilogx(r_arr,psffit_arr,'k','Displayname','analytic fit');

set(gca,'XScale','log','YScale','lin');
ylim([-0.1,1.1])
xlim([4e-1,7e2])
xlabel('arcsec')
ylabel('<I_{stack}>')
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff

subplot(1,3,3)
errorbar(r_arr,psf_arr,epsf_arr,'ro-');hold on
errorbar(r_arr,-psf_arr,epsf_arr,'bo-');hold on
plot(r_arr,profstamp_arr,'color',[0.00  0.50  0.00],'linewidth',2);
plot(r_arr,psffit_arr,'k');
set(gca,'XScale','log','YScale','log');
ylim([1e-5,1.1])
xlim([4e-1,7e2])
xlabel('arcsec')
ylabel('<I_{stack}>')

print(savename,'-dpng');