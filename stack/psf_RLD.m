scale = 10; % subpixel is 'scale' times smaller than origial pixel
Npix = 30; % the stamp extend Npix pixels in each direction
Nsub = Npix*scale; % stacking stamp: (Nsub*2+1)**2 subpixels
sig = 0.8; % w.r.t. large pixel
sig_noise = 1e-1;
Nstack = 1000; % number of random source stacking

rr = make_radius_map(zeros(2*Nsub+1,2*Nsub+1),Nsub+1,Nsub+1);
psf_beam = exp(-rr.^2./2./(sig*scale).^2)./(2*pi*(sig*scale)^2);
[xx,yy] = meshgrid(-Nsub:Nsub,-Nsub:Nsub);
psf_pix = (10 - abs(xx)).*(10 - abs(yy));
psf_pix(abs(xx)>=10 | abs(yy) >= 10) = 0;
psf_pix = psf_pix ./ sum(psf_pix(:));
psf_tot = conv2(psf_beam, psf_pix,'same');

stack = zeros(2*Nsub+1);
src_coord=Nsub+0.5+10*rand(2,Nstack);
for isim = 1:Nstack
    xsrc=src_coord(1,isim);
    ysrc=src_coord(2,isim);
    radmap = make_radius_map(zeros(2*Nsub+10),xsrc,ysrc);
    psfmap = exp(-radmap.^2./2./(sig*scale).^2)./(2*pi*(sig*scale)^2); 
    psfmap_coarse=rebin_map_coarse(psfmap,10).*scale.^2;
    psfmap_coarse = psfmap_coarse + normrnd(0,sig_noise, size(psfmap_coarse));
    psfmap_fine=imresize(psfmap_coarse./scale.^2,10,'method','nearest');
    stack=stack+psfmap_fine(round(xsrc)-Nsub:round(xsrc)+Nsub,...
                            round(ysrc)-Nsub:round(ysrc)+Nsub);
end
stack=stack./Nstack;

figure
setwinsize(gcf,1200,300)
subplot(1,3,1)
imageclip(psf_beam(Nsub+1-30: Nsub+1+30,Nsub+1-30: Nsub+1+30));
caxis([0,0.005])
subplot(1,3,2)
imageclip(psf_pix(Nsub+1-30: Nsub+1+30,Nsub+1-30: Nsub+1+30));
caxis([0,0.005])
subplot(1,3,3)
imageclip(psf_tot(Nsub+1-30: Nsub+1+30,Nsub+1-30: Nsub+1+30));
caxis([0,0.005])

psf_rl = deconvlucy(stack,psf_pix);

figure
setwinsize(gcf,1200,300)
subplot(1,4,1)
imageclip(psf_beam(Nsub+1-30: Nsub+1+30,Nsub+1-30: Nsub+1+30));
caxis([0,0.005])
subplot(1,4,2)
imageclip(stack(Nsub+1-30: Nsub+1+30,Nsub+1-30: Nsub+1+30));
caxis([0,0.005])
subplot(1,4,3)
imageclip(stackcb(Nsub+1-30: Nsub+1+30,Nsub+1-30: Nsub+1+30));
caxis([0,0.005])
title('CIBER stack')
subplot(1,4,4)
imageclip(psf_rl(Nsub+1-30: Nsub+1+30,Nsub+1-30: Nsub+1+30));
caxis([0,0.005])
%%
flight=40030;
mypaths=get_paths(flight);
inst = 2;
ifield = 8;
dt = get_dark_times(flight,inst,ifield);
psfdir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf/TM',...
        num2str(inst),'/');
stamper = fitsread(strcat(psfdir,dt.name,'_stamper.fits'));
hitmap = fitsread(strcat(psfdir,dt.name,'_hitmap.fits'));
mapcb = stamper./hitmap;
mapcb = mapcb(1201-300:1201+300,1201-300:1201+300);
mapcb = mapcb./sum(mapcb(:));