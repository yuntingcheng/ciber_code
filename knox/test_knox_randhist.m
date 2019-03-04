%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%See histogram before and after real space FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
pixscale=7;

darkstatdir=strcat('/Users/ytcheng/ciber/doc/',...
         '20160906_NoiseRealization/darkstat/',num2str(flight),'/');
load(strcat(darkstatdir,'TM',num2str(inst),'_mdarkps'),'mdarkps');
load(strcat(darkstatdir,'TM',num2str(inst),'_darkps'),'darkps');

savedir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/test_knox/rand_level/';

alldatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(alldatdir,'band',num2str(inst),'_alldat'),'alldat');

bigmask=alldat(8).bigmask;

%%% get knox error
[Cl,~,~,l,~,dClknox,binl,dl] = Cl_from_Cl2d_clip(ones(1024),pixscale);
dClknox=dClknox./Cl;dClknox(find(dClknox~=dClknox))=0;
%% 
load(strcat(savedir,'hist/xpixuse_arr'),'xpixuse_arr');
load(strcat(savedir,'hist/ypixuse_arr'),'ypixuse_arr');

sig_arr=[0,0.25,0.5,0.75,1];
nsim=10000;
xpix_arr=xpixuse_arr;
ypix_arr=ypixuse_arr;
npix=numel(xpix_arr);
%%
for isig=1:numel(sig_arr)
sig=sig_arr(isig);

in_arr=zeros(nsim,npix);
out_arr=zeros(nsim,npix);
inavg_arr=zeros(1024);
inavg2_arr=zeros(1024);
outavg_arr=zeros(1024);
outavg2_arr=zeros(1024);
for i=1:nsim
    if ~rem(i,500)
        disp(sprintf('isig=%d,i=%d',isig,i));
    end
    
    randmap=(normrnd(0,sig,1024).^2+normrnd(0,sig,1024).^2)./2-sig^2+1;
    Cl2din=randmap.*darkps(8).Clf2d_ave;
    %%%%%%%%%%%
    %Cl2din=randmap;% for white noise case
    %%%%%%%%%%%
    rnmap=ifft2(fftshift(sqrt(Cl2din)).*fft2(normrnd(0,1,1024)));
    rnmap=rnmap'./(pixscale/3600.0*pi/180.0);
    rnmap=real(rnmap)./abs(real(rnmap)).*abs(rnmap);
    [Cl,l,~,~,~,~,rCl2d]=get_angular_spec(rnmap,rnmap,pixscale);
    inavg_arr=inavg_arr+Cl2din./nsim;
    inavg2_arr=inavg2_arr+Cl2din.^2./nsim;
    outavg_arr=outavg_arr+rCl2d./nsim;
    outavg2_arr=outavg2_arr+rCl2d.^2./nsim;
    
    for j=1:npix
        in_arr(i,j)=Cl2din(xpix_arr(j),ypix_arr(j));
        out_arr(i,j)=rCl2d(xpix_arr(j),ypix_arr(j));
    end
end
data(isig).sig=sig;
data(isig).inavg_arr=inavg_arr;
data(isig).inavg2_arr=inavg2_arr;
data(isig).outavg_arr=outavg_arr;
data(isig).outavg2_arr=outavg2_arr;
data(isig).xpix_arr=xpix_arr;
data(isig).ypix_arr=ypix_arr;
data(isig).in_arr=in_arr;
data(isig).out_arr=out_arr;
end

%save(strcat(savedir,'hist/data'),'data');
%save(strcat(savedir,'hist/dataw'),'data');
save(strcat(savedir,'hist/data1'),'data');
%%
