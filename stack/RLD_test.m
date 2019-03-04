savedir=(strcat('/Users/ytcheng/ciber/doc/20171130_psfstack/plots/RLD/'));
npix=100;
Nsamp=2000;
w = zeros(2*npix+1,2*npix+1);
w(npix+1-5:npix+1+5,npix+1-5:npix+1+5) = 1;
%% Gaussian PSF

sig_arr = [0.1,0.3,0.5,1,2];
for isig=1:numel(sig_arr)
sig = sig_arr(isig);
radmap = make_radius_map(zeros(2*npix+1),npix+1,npix+1).*0.1;
psf = exp(-radmap.^2./2./sig.^2)./(2*pi*sig^2).*(0.1^2);

data=zeros(2*npix+1,2*npix+1);
src_coord=npix+0.5+10*rand(2,Nsamp);
for i=1:Nsamp
    xsrc=src_coord(1,i);
    ysrc=src_coord(2,i);
    radmap = make_radius_map(zeros(2*npix+10),xsrc,ysrc).*0.1;
    psfmap = exp(-radmap.^2./2./sig.^2)./(2*pi*sig^2).*(0.1^2);
    
    psfmap_coarse=rebin_map_coarse(psfmap,10);
    psfmap_fine=imresize(psfmap_coarse,10,'method','nearest');

    data=data+psfmap_fine(round(xsrc)-npix:round(xsrc)+npix,...
                            round(ysrc)-npix:round(ysrc)+npix);
end
data=data./Nsamp;

figure
profile = radial_prof(psf,ones(2*npix+1),npix+1,npix+1,1,20);
semilogx(profile.r.*0.1,profile.prof./profile.prof(1),'k--',...
    'linewidth',3,'Displayname','input PSF');hold on
profile = radial_prof(data,ones(2*npix+1),npix+1,npix+1,1,20);
semilogx(profile.r.*0.1,profile.prof./profile.prof(1),'k',...
    'linewidth',3,'Displayname','stacking');

u_t = data;
for niter=[1,10,100,500,1000]
    for i=1:niter
        c = conv2(u_t,w,'same');
        ratio = data./c;
        ratio(find(c==0))=0;
        u_t1 = u_t.*conv2(ratio,w,'same');
        u_t = u_t1;
    end
    profile = radial_prof(u_t,ones(2*npix+1),npix+1,npix+1,1,20);
    semilogx(profile.r.*0.1,profile.prof./profile.prof(1),...
        'linewidth',1.5,'Displayname',sprintf('R-L %d iters',niter));
    drawnow;
end
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff
title(strcat('\sigma = ',num2str(sig),' pixel size'),'fontsize',15);
ylim([-0.1,1.1])
xlim([7e-2,2e1])

xlabel('pixel size','fontsize',15);
ylabel('<I>','fontsize',15)

savename=strcat(savedir,'rld_sig',num2str(isig));
print(savename,'-dpng');
end