%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the fit_res results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=2;
npix=400;
pixsize=0.7;
quad_arr=['A','B','C','D'];
mypaths=get_paths(flight);

psfdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(psfdir,'fitpsfdat'),'fitpsfdat');
loaddir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/resmap/TM',...
    num2str(inst),'/');
savedir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/resmap/fitdat/TM',...
    num2str(inst),'/');
%%% 2M sim map with new offseted PSF
loaddir1=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/resmap/offset/TM',...
    num2str(inst),'/');
%% plot the data and best fit residual map
load(strcat(savedir,'bestpar'),'bestpar');

fitrad=15;
radmap = make_radius_map(zeros(2*npix+1),npix,npix).*pixsize;
sp=find(radmap<fitrad);
mask=zeros(size(radmap));
mask(sp)=1;

for ifield=4:8
figure
setwinsize(gcf,1200,800)

for iquad=1:4

dt=get_dark_times(flight,inst,ifield);
quad=quad_arr(iquad);
resdat0=fitsread(strcat(loaddir,dt.name,'_',quad,'_','resmap0.fits'));
resdat=fitsread(strcat(loaddir,dt.name,'_',quad,'_','resmap.fits'));
resdat1=fitsread(strcat(loaddir1,dt.name,'_',quad,'_','resmap.fits'));
bestparam_norm=fitpsfdat(ifield).bestparam_norm;
A=bestparam_norm(1);
B=bestparam_norm(2);
sig=bestparam_norm(3);
r0=bestparam_norm(4);
alpha=bestparam_norm(5);
chi2best=bestparam_norm(6);

%%% 2MASS sim map PSF
psfmap = A*exp(-radmap.^2./2./sig^2)+B./(1+(radmap./r0).^alpha);

dx=bestpar(ifield).quad(iquad).dx;
dy=bestpar(ifield).quad(iquad).dy;
rs=bestpar(ifield).quad(iquad).rs;

radmap1 = make_radius_map(zeros(2*npix+1),npix-dy,npix-dx).*pixsize;
psfmap1 = A*exp(-radmap1.^2./2./(sig*rs)^2)+...
    B./(1+(radmap1./(r0*rs)).^alpha);
psfmapn = A*exp(-radmap.^2./2./(sig*rs)^2)+B./(1+(radmap./(r0*rs)).^alpha);
norm1=sum(psfmapn(:));
psfmap1=psfmap1./norm1;

resmap=psfmap1-psfmap;
resmod=resmap./psfmap(401,401);

v=[-0.4 0.4];
pltmin=npix+1-25;pltmax=npix+1+25;
subplot(4,4,iquad)
imageclip(mask(pltmin:pltmax,pltmin:pltmax).*resdat0(pltmin:pltmax,pltmin:pltmax));
caxis(v)
title(strcat(quad,'\_residual'))
subplot(4,4,iquad+4)
imageclip(mask(pltmin:pltmax,pltmin:pltmax).*resdat(pltmin:pltmax,pltmin:pltmax));
title(strcat(quad,'\_residual(scaled)'))
caxis(v)
subplot(4,4,iquad+8)
imageclip(mask(pltmin:pltmax,pltmin:pltmax).*resmod(pltmin:pltmax,pltmin:pltmax));
caxis(v)
title(strcat(quad,'\_residual model'))
subplot(4,4,iquad+12)
imageclip(mask(pltmin:pltmax,pltmin:pltmax).*resdat1(pltmin:pltmax,pltmin:pltmax));
caxis(v)
title(strcat(quad,'\_residual new'))

end
figtitle(strcat('TM',num2str(inst),'\_',dt.name),'fontweight','bold');

imname=strcat(savedir,'resfit_TM',num2str(inst),'_',dt.name);
print(imname,'-dpng');close
end