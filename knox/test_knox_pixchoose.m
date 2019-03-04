%%%%%%%%%%%%%%%%%%%%%%%%
%select the Cl2d pix to look at their histogram
%%%%%%%%%%%%%%%%%%%%%%%%
loaddir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/test_knox/rand_level/';
pixscale=7;
[Cl,~,~,l,~,dClknox,binl,dl] = Cl_from_Cl2d_clip(ones(1024),pixscale);

sig_arr=[0,0.25,0.5,0.75,1];
load(strcat(loaddir,'hist/dataw'),'data');
for isig=5%1:numel(sig_arr)
wCl2din=data(1).inavg_arr;
wsig=data(isig).sig;
winavg_arr=data(isig).inavg_arr;
winavg2_arr=data(isig).inavg2_arr;
woutavg_arr=data(isig).outavg_arr;
woutavg2_arr=data(isig).outavg2_arr;

winstd=sqrt(winavg2_arr-winavg_arr.^2);
woutstd=sqrt(woutavg2_arr-woutavg_arr.^2); 
end
%%
load(strcat(loaddir,'hist/data'),'data');
sig_arr=[0,0.25,0.5,0.75,1];
for isig=5%1:numel(sig_arr)
Cl2din=data(1).inavg_arr;
sig=data(isig).sig;
inavg_arr=data(isig).inavg_arr;
inavg2_arr=data(isig).inavg2_arr;
outavg_arr=data(isig).outavg_arr;
outavg2_arr=data(isig).outavg2_arr;
in_arr=data(isig).in_arr;
out_arr=data(isig).out_arr;

instd=sqrt(inavg2_arr-inavg_arr.^2);
outstd=sqrt(outavg2_arr-outavg_arr.^2); 
end
%% plot the mean hist
figure
setwinsize(gcf,1200,600)
subplot(2,3,1)
histogram(winavg_arr./wCl2din);
title('<Cl2Drnd>_{sim}/Cl2Din --white');
subplot(2,3,2)
histogram(woutavg_arr./winavg_arr);
title('<Cl2Dmap>_{sim}/<Cl2Drnd>_{sim} --white');
subplot(2,3,3)
histogram(woutavg_arr./wCl2din);
title('<Cl2Dmap>_{sim}/Cl2Din --white');
subplot(2,3,4)
histogram(inavg_arr./Cl2din);
title('<Cl2Drnd>_{sim}/Cl2Din --RN');
subplot(2,3,5)
histogram(outavg_arr./inavg_arr);
title('<Cl2Dmap>_{sim}/<Cl2Drnd>_{sim} --RN');
subplot(2,3,6)
histogram(outavg_arr./Cl2din);
title('<Cl2Dmap>_{sim}/Cl2Din --RN');
%savename=sprintf('%shist/mean_hist',savedir);
%print(savename,'-dpng');close
%% plot mean map
figure
setwinsize(gcf,1200,500)
subplot(1,2,1)
imageclip(inavg_arr./Cl2din);
title('<Cl2drnd>_{sim}/Cl2din --RN');
subplot(1,2,2)
imageclip(outavg_arr./Cl2din);
title('<Cl2dmap>_{sim}/Cl2din --RN');
%savename=sprintf('%shist/mean_map',savedir);
%print(savename,'-dpng');%close
%% plot std
figure
setwinsize(gcf,1000,1000)
subplot(2,2,1)
imageclip(instd./Cl2din);
title('sqrt(<Cl2drnd^2>_{sim}-<Cl2drnd>^2_{sim})/Cl2din --RN');
subplot(2,2,2)
histogram(instd./Cl2din);

subplot(2,2,3)
imageclip(outstd./Cl2din);
title('sqrt(<Cl2dmap^2>_{sim}-<Cl2dmap>^2_{sim})/Cl2din --RN');
subplot(2,2,4)
histogram(outstd./Cl2din);
%savename=sprintf('%shist/std_map',savedir);
%print(savename,'-dpng');%close
%% choose the pix
ell=get_l(1024,1024,pixscale,1);
jointmask=ones(1024);
xpixuse_arr=[];ypixuse_arr=[];vpixuse_arr=[];

for i=[20,26,29]
mask=zeros(1024);mask((ell >= binl(i)) & (ell <= binl(i+1)))=1;
[xplot,yplot]=find(mask);

jointmask=jointmask.*~mask;

masksort=mask;
masksort(513,:)=0;masksort(:,513)=0;
masksort(1,:)=0;masksort(:,1)=0;
masksort(1024,:)=0;masksort(:,1024)=0;

dat=Cl2din(find(masksort));
% #1: max(Cl2din)
[xpix,ypix]=find(Cl2din==max(dat));ind=find(Cl2din==max(dat));
xpix_arr(1)=xpix(1);ypix_arr(1)=ypix(1);value_arr(1)=max(dat);
masksort(ind)=0;

% #3: min(Cl2din)
[xpix,ypix]=find(Cl2din==min(dat));ind=find(Cl2din==min(dat));
xpix_arr(3)=xpix(1);ypix_arr(3)=ypix(1);value_arr(3)=min(dat);
masksort(ind)=0;

% numel(dat) must odd numbers so the median won't be avg of two values
if ~rem(numel(dat),2);dat=dat(1:end-1);end 
% #2: median(Cl2din)
[xpix,ypix]=find(Cl2din==median(dat));ind=find(Cl2din==median(dat));
xpix_arr(2)=xpix(1);ypix_arr(2)=ypix(1);value_arr(2)=median(dat);
masksort(ind)=0;

dat=outavg_arr./Cl2din;dat=dat(find(masksort));
% #4: max(outavg_arr./Cl2din)
[xpix,ypix]=find(outavg_arr./Cl2din==max(dat));
ind=find(outavg_arr./Cl2din==max(dat));
xpix_arr(4)=xpix(1);ypix_arr(4)=ypix(1);value_arr(4)=max(dat);
masksort(ind)=0;

if ~rem(numel(dat),2);dat=dat(1:end-1);end 
% #5: median(outavg_arr./Cl2din)
[xpix,ypix]=find(outavg_arr./Cl2din==median(dat));
ind=find(outavg_arr./Cl2din==median(dat));
xpix_arr(5)=xpix(1);ypix_arr(5)=ypix(1);value_arr(5)=median(dat);
masksort(ind)=0;

%6: max outstd./Cl2din
dat=outstd./Cl2din;dat=dat(find(masksort));
[xpix,ypix]=find(outstd./Cl2din==max(dat));
xpix_arr(6)=xpix(1);ypix_arr(6)=ypix(1);value_arr(6)=max(dat);

xpixuse_arr=[xpixuse_arr,xpix_arr];
ypixuse_arr=[ypixuse_arr,ypix_arr];
vpixuse_arr=[vpixuse_arr,value_arr];
end
%save(strcat(savedir,'hist/xpixuse_arr'),'xpixuse_arr');
%save(strcat(savedir,'hist/ypixuse_arr'),'ypixuse_arr');
%save(strcat(savedir,'hist/vpixuse_arr'),'vpixuse_arr');
%% show the 3 ell modes
figure
setwinsize(gcf,1200,500)
subplot(1,2,1)
imageclip(Cl2din);
subplot(1,2,2)
imageclip(Cl2din.*jointmask);