%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test different level of random amplitude
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
sig_arr=[0,0.25,0.5,0.75,1];
%%
for isig=5%1:numel(sig_arr)
sig=sig_arr(isig)
%%% get the weight
rCl2d_arr=zeros(100,1024,1024);
for i=1:100
    i
    randmap=(normrnd(0,sig,1024,1024).^2+normrnd(0,sig,1024,1024).^2)./2-sig^2+1;
    Cl2din=randmap.*darkps(8).Clf2d_ave;
    rnmap=ifft2(fftshift(sqrt(Cl2din)).*fft2(normrnd(0,1,1024)));
    rnmap=rnmap'./(pixscale/3600.0*pi/180.0);
    rnmap=real(rnmap)./abs(real(rnmap)).*abs(rnmap);
    [Cl,l,~,~,~,~,rCl2d]=get_angular_spec(rnmap,rnmap,pixscale);
    rCl2d_arr(i,:,:)=rCl2d;
end
rCl2dave=squeeze(mean(rCl2d_arr));rCl2dstd=squeeze(std(rCl2d_arr));
clear rCl2d_arr;

%% 100 sims
rCl_arr=zeros(100,29);
rClw1_arr=zeros(100,29);
rClw2_arr=zeros(100,29);
rClw3_arr=zeros(100,29);

sCl_arr=zeros(100,29);
sClw1_arr=zeros(100,29);
sClw2_arr=zeros(100,29);
sClw3_arr=zeros(100,29);
for i=1:100
    i
    randmap=(normrnd(0,sig,1024,1024).^2+normrnd(0,sig,1024,1024).^2)./2-sig^2+1;
    Cl2din=randmap.*darkps(8).Clf2d_ave;
    rnmap=ifft2(fftshift(sqrt(Cl2din)).*fft2(normrnd(0,1,1024)));
    rnmap=rnmap'./(pixscale/3600.0*pi/180.0);
    rnmap=real(rnmap)./abs(real(rnmap)).*abs(rnmap);
    [rCl,l,~,~,~,~,rCl2d]=get_angular_spec(rnmap,rnmap,pixscale);
    [rClw1] = Cl_from_Cl2d_clip(rCl2d,pixscale,'w',(1./(rCl2dstd))');
    [rClw2] = Cl_from_Cl2d_clip(rCl2d,pixscale,'w',(1./(rCl2dave))');
    [rClw3] = Cl_from_Cl2d_clip(rCl2d,pixscale,'w',(1./(rCl2dstd+rCl2dave))');

    rCl_arr(i,:)=rCl;
    rClw1_arr(i,:)=rClw1;
    rClw2_arr(i,:)=rClw2;    
    rClw3_arr(i,:)=rClw3; 
    rCl2d_arr(i,:,:)=rCl2d;
    
    smap=ifft2(fftshift(sqrt(randmap)).*fft2(normrnd(0,1,1024)));
    smap=smap'./(pixscale/3600.0*pi/180.0);
    smap=real(smap)./abs(real(smap)).*abs(smap);
    
    [sCl,l,~,~,~,~,sCl2d]=get_angular_spec(smap,smap,pixscale);
    [sClw1] = Cl_from_Cl2d_clip(sCl2d,pixscale,'w',(1./(rCl2dstd))');
    [sClw2] = Cl_from_Cl2d_clip(sCl2d,pixscale,'w',(1./(rCl2dave))');
    [sClw3] = Cl_from_Cl2d_clip(sCl2d,pixscale,'w',(1./(rCl2dstd+rCl2dave))');
    sCl_arr(i,:)=sCl;
    sClw1_arr(i,:)=sClw1;
    sClw2_arr(i,:)=sClw2;    
    sClw3_arr(i,:)=sClw3; 
    
end

figure
loglog(l,dClknox,'k');hold on
loglog(l,std(rCl_arr)./mean(rCl_arr),'r.','markersize',10);
loglog(l,std(rClw1_arr)./mean(rClw1_arr),'ro','markersize',5);
loglog(l,std(rClw2_arr)./mean(rClw2_arr),'r+','markersize',5);
loglog(l,std(rClw3_arr)./mean(rClw3_arr),'rs','markersize',5);
loglog(l,std(sCl_arr)./mean(sCl_arr),'b.','markersize',10);
loglog(l,std(sClw1_arr)./mean(sClw1_arr),'bo','markersize',5);
loglog(l,std(sClw2_arr)./mean(sClw2_arr),'b+','markersize',5);
loglog(l,std(sClw3_arr)./mean(sClw3_arr),'bs','markersize',5);

xlabel('$\ell$','interpreter','latex','fontsize',18);
ylabel('$\delta C_\ell/\bar{C_\ell}$',...
        'interpreter','latex','fontsize',18);
title(sprintf('StD=%.2f',sig));
xlim([1e2,2e5]);ylim([1e-3,2])
legend({'Knox','Cl2d->rand->map->FFT','Cl2d->rand->map->FFT->weight(std)',...
'Cl2d->rand->map->FFT->weight(avg)','Cl2d->rand->map->FFT->weight(avg+std)'},...
'location','southwest');
savename=strcat(sprintf('%ssig%d',savedir,isig));
print(savename,'-dpng');%close
end

