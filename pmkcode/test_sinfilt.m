fitparams = get_default_fitparams;
fitparams.flight = 40030;
fr1 = get_data_frames(1,'raildark',fitparams);

[filtmap,dumbmap,sigmask,filtfrq,ts] = imager_filtts(fr1,'sigma',2);
%%
s=size(fr1);
n=s(1)/2

raw1 = fr1(1:n,:,:);
raw2 = fr1(n+1:2*n,:,:);
filt1 = filtfrq(1:n,:,:);
filt2 = filtfrq(n+1:2*n,:,:);
%%

[raw1m] = fastlinefit_frin(raw1,0,0,1000);
[raw2m] = fastlinefit_frin(raw2,0,0,1000);
[filt1m] = fastlinefit_frin(filt1,0,0,1000);
[filt2m] = fastlinefit_frin(filt2,0,0,1000);
%%

load('/home/pkorngut/projects/projects2015/Jul/Jul30_2015_imager_drcombine/fullpipe_Jul30_dr.mat');

f=8;
bigmask = alldat(f).bigmask;
bigmask(10:70,380:520)=0;
%%
rawdiff = raw1m-raw2m;
filtdiff = filt1m - filt2m;
nbins=100;

[Clraw,l,lth,dCl,binl,dl,Clraw2d] = ...
    get_angular_spec(rawdiff.*sigmask.*bigmask,rawdiff.*sigmask.*bigmask,7.0,nbins,1,ones(1024),'verbose',0,'superbin',0);

yraw = (l.*(l+1).*Clraw).^(1/2);

[Clfilt,l,lth,dCl,binl,dl,Clfilt2d] = ...
    get_angular_spec(filtdiff.*sigmask.*bigmask,filtdiff.*sigmask.*bigmask,7.0,nbins,1,ones(1024),'verbose',0,'superbin',0);
yfilt = (l.*(l+1).*Clfilt).^(1/2);

%%

subplot(1,2,1)
loglog(l,yraw,'Color','red','linewidth',2)
hold on
loglog(l,yfilt,'Color','green','linewidth',2)
hold off
set(gca,'FontSize',22,'Ytick',10.^(-2:0))
axis([1e2,1e5,1e-3,2])
xlabel('l')
ylabel('(l(l+1)Cl)^{1/2}')

dx=50;
cx=512;
ax =[cx-dx,cx+dx,cx-dx,cx+dx];
cax=[0,2e-9];

subplot(2,2,2)
imageclip(Clraw2d);
colormap jet
axis(ax)
caxis(cax)
title('unfiltered')

subplot(2,2,4)
imageclip(Clfilt2d);
title('filtered')
axis(ax)
caxis(cax)
%%
cax = ([-.6,.6]);

subplot(1,2,1)
imageclip(rawdiff.*sigmask);
title('Raw Difference')
caxis(cax)

subplot(1,2,2)
imageclip(filtdiff.*sigmask);
caxis(cax)
title('Filtered Difference')

%%
sig = 3;
smraw = fillpadsmooth(rawdiff,sigmask,sig);
smfilt = fillpadsmooth(filtdiff,sigmask,sig);
%%
cax = ([-.5,.5]);

subplot(1,2,1)
imageclip(smraw.*sigmask);
title('Raw Difference')
caxis(cax)

subplot(1,2,2)
imageclip(smfilt.*sigmask);
caxis(cax)
title('Filtered Difference')

