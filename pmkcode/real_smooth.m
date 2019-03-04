%%
im=fitsread_orient('/home/pkorngut/projects/projects2015/Jul/Jul22_2015_starmasks/instmasks/elat10instmask.fits');

load('/home/pkorngut/projects/projects2017/Sep/Sep29_2017_yunting_smooth/cbmap.mat');
meas =10.^a;
load('/home/pkorngut/projects/projects2017/Sep/Sep29_2017_yunting_smooth/ukmap.mat');
uk=10.^b;
mask = ones(1024);
death = find(a == 0);
mask(death)=0;
%%


%%
mask(:,110:140)=0;
sm = sigmaclipmask(meas,1.5,5);
mask=mask.*sm.*im;
death = find(abs(meas) > 200);
mask(death) = 0;

use = find(mask);
m = mean(meas(use));
meas=meas-m;
%%
ftmap = fft2(meas.*mask);
filt = ones(1024);
filt(1:50,:)=0;
filt(1024-50:1024,:)=0;
filtmodel = ftmap.*filt;
filtmap = imrotate(real(fft2(filtmodel)),180);



%%
sig = 15;
%smmeas = fillpadsmooth(meas,mask,sig);
smmeas = fillpadsmooth(filtmap,mask,sig);

smuk = fillpadsmooth(uk,mask,sig);
cax = [-60,60];

subplot(1,2,1)
imageclip(smmeas);
%caxis(cax)
subplot(1,2,2)
imageclip(smuk);
%caxis(cax)

%%
cax = [-60,60];

conts = 1:5:50;
subplot(1,1,1)
imagesc(smmeas);
caxis(cax*1e5)
axis square
colormap (hot)
hold on
contour(smuk,conts,'color','green','linewidth',2)