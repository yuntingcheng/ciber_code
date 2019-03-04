%%% get CIBER pix size psf
flight=40030;
inst=1;
pixscale=7;

for ifield=8%1:8
dt=get_dark_times(flight,inst,ifield);
load(strcat('/Volumes/HD1TB/CIBER/data/',num2str(flight),...
    '/dr/TM',num2str(inst),'_',dt.name,'_dr150206.mat'));
%% average over phi to get the 1D PSF
psf = data.psf.psf;
prof = radial_prof_pk(psf,psf*0+1,51,51,100);
%% find the best fit param for fitting function
x = prof.r.*data.psf.pixscale;

sig_arr=3:0.1:6;
amp_arr=5:1:80;
chi2_arr=zeros(numel(sig_arr),numel(amp_arr));
bestpar=[sig_arr(1),amp_arr(1)];
bestchi2=1e8;

for isig=1:numel(sig_arr)
    sig=sig_arr(isig);
    for iamp=1:numel(amp_arr)
        amp=amp_arr(iamp);
        modl1=exp(-(x.^2)./(2*sig^2));
        modl2 = amp*x.^(-3);
        low = find(x <= 8);
        modl2(low) =0;
        modl=modl1+modl2;
        chi2=nansum((log(prof.prof)-log(modl)).^2./(log(prof.err).^2));
        %chi2=nansum((log(prof.prof)-log(modl)).^2);
        chi2_arr(isig,iamp)=chi2;
        
        if chi2<bestchi2
            bestchi2=chi2;
            bestpar=[sig,amp];
        end
        
    end
end
%figure
%imageclip(chi2_arr);
%% 
npix=floor(data.psf.pixscale*size(data.psf.psf,1)/pixscale);
if mod(npix,2)==0;npix=npix-1;end
psfmap=zeros(npix,npix);
radius_map=make_radius_map(psfmap,(npix+1)/2,(npix+1)/2);
for i=1:npix
    for j=1:npix
        radius=radius_map(i,j)*pixscale;
        m1=exp(-(radius.^2)./(2*bestpar(1)^2));
        m2=bestpar(2)*radius.^(-3);
        if radius>x(8)
            psfmap(i,j)=m1+m2;
        else
            psfmap(i,j)=m1;
        end
    end
end

figure
setwinsize(gcf,1000,500)
subplot(1,2,1)

semilogx(x,prof.prof);hold on

modl1=exp(-(x.^2)./(2*bestpar(1)^2));
modl2 = bestpar(2)*x.^(-3);
low = find(x <= 8);
modl2(low) =0;
modl=modl1+modl2;
semilogx(x,modl,'-.');

cent=(npix+1)./2;
semilogx((1:cent-1).*pixscale,psfmap(cent,cent+1:npix),'o');

axis([1,1000,1e-5,2])
legend({'data','fitting func','CIBER repix'})

savedir='/Users/ytcheng/ciber/doc/20170320_FFsim/psf/';
save(strcat(savedir,'psf_',dt.name),'psfmap');
%%%%%%%%%%% test convolution %%%%%%%%%%%%%%%
map=randn(1024);
C=conv2(map,psfmap./sum(psfmap(:)));
Nout=(size(C,1)-size(map,1))/2;
convmap=C(Nout+1:Nout+1024,Nout+1:Nout+1024);
% get bl from MZ 
[Cl,l]=get_angular_spec(map,map,pixscale);
[Clc]=get_angular_spec(convmap,convmap,pixscale);
logbl=interp1(log10(data.psf.l),log10(data.psf.bl),log10(l),'spline');
blmz=10.^logbl;

subplot(1,2,2)
loglog(l,l.*(l+1).*Cl);hold on
loglog(l,l.*(l+1).*Clc);
loglog(l,l.*(l+1).*Clc./blmz);
legend({'input','input conv','bl correction'},'location','southeast')
end