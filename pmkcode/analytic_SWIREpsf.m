pixscale = 7;
Npix=100;
[xx,yy] = meshgrid((-Npix/2:1:Npix/2).*pixscale);
rr = sqrt(xx.^2+yy.^2);

sig=10.0/2.35;
gau = exp(-(rr.^2)./(2*sig^2));
m2 = 350./rr.^(3.7);
cut=9.0;
sp = find(rr < cut);
m2(sp)=0;
sp = find(rr > cut);
gau(sp)=0;

psf = gau+m2;
norm=2.7199; % get from very large Npix
psf=psf./norm;
%%
map=randn(1024);
C=conv2(map,psf);
Nout=(size(C,1)-size(map,1))/2;
convmap=C(Nout+1:Nout+1024,Nout+1:Nout+1024);
[Cl,l]=get_angular_spec(map,map,pixscale);
[Clc]=get_angular_spec(convmap,convmap,pixscale);

figure
setwinsize(gcf,1000,500)
subplot(1,2,1)
inmap =psf-mean(psf(:));
[Clb,l] =get_angular_spec(inmap,inmap,pixscale);
loglog(l,Clb./max(Clb),'o-');hold on
loglog(l,Clc./Cl);hold on
legend({'PSF FFT','white ratio'},'location','southwest')

subplot(1,2,2)
loglog(l,l.*(l+1).*Cl);hold on
loglog(l,l.*(l+1).*Clc);
legend({'input','input conv'},'location','southeast')
%%
bl_arr=zeros(100,29);
for i=1:100
    i
map=randn(1024);
C=conv2(map,psf);
Nout=(size(C,1)-size(map,1))/2;
convmap=C(Nout+1:Nout+1024,Nout+1:Nout+1024);
[Cl,l]=get_angular_spec(map,map,pixscale);
[Clc]=get_angular_spec(convmap,convmap,pixscale);
bl_arr(i,:)=Clc./Cl;
end
%%
blSWIRE=squeeze(mean(bl_arr));
psfmap=psf;
loglog(l,blSWIRE);hold on
loglog(l,Clb./max(Clb),'o-');hold on
save('/Users/ytcheng/ciber/code/sim/blSWIRE','blSWIRE');
save('/Users/ytcheng/ciber/code/sim/psfSWIRE','psfmap');