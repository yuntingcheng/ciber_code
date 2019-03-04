%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use the weight from mean or std of 100 sim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
pixscale=7;

darkstatdir=strcat('/Users/ytcheng/ciber/doc/',...
         '20160906_NoiseRealization/darkstat/',num2str(flight),'/');
load(strcat(darkstatdir,'TM',num2str(inst),'_mdarkps'),'mdarkps');
load(strcat(darkstatdir,'TM',num2str(inst),'_darkps'),'darkps');

savedir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/test_knox/weight_from_sim/';

alldatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(alldatdir,'band',num2str(inst),'_alldat'),'alldat');

bigmask=alldat(8).bigmask;

%%% get knox error
[Cl,~,~,l,~,dClknox,binl,dl] = Cl_from_Cl2d_clip(ones(1024),pixscale);
dClknox=dClknox./Cl;dClknox(find(dClknox~=dClknox))=0;
%% generate weight
mCl2d_arr=zeros(100,1024,1024);
rCl2d_arr=zeros(100,1024,1024);
for i=1:100
    i
    rnmap = readnoise_realization(darkps(8).Clf2d_ave,pixscale,'norand',1);
    mrnmap=rnmap.*bigmask;mrnmap=mrnmap-mean(mrnmap(find(mrnmap)));
    mrnmap=mrnmap.*bigmask;
    [Cl,l,~,~,~,~,mCl2d]=get_angular_spec(mrnmap,mrnmap,pixscale);
    mCl2d_arr(i,:,:)=mCl2d;
    
    rnmap = readnoise_realization(darkps(8).Clf2d_ave,pixscale);
    [Cl,l,~,~,~,~,rCl2d]=get_angular_spec(rnmap,rnmap,pixscale);
    rCl2d_arr(i,:,:)=rCl2d;
end

mCl2dave=squeeze(mean(mCl2d_arr));mCl2dstd=squeeze(std(mCl2d_arr));
rCl2dave=squeeze(mean(rCl2d_arr));rCl2dstd=squeeze(std(rCl2d_arr));

clear mCl2d_arr;clear rCl2d_arr;
%% 
mCl_arr=zeros(100,29);
mClw1_arr=zeros(100,29);
mClw2_arr=zeros(100,29);
mClw3_arr=zeros(100,29);

rCl_arr=zeros(100,29);
rClw1_arr=zeros(100,29);
rClw2_arr=zeros(100,29);
rClw3_arr=zeros(100,29);

for i=1:100
    i
    rnmap = readnoise_realization(darkps(8).Clf2d_ave,pixscale,'norand',1);
    mrnmap=rnmap.*bigmask;mrnmap=mrnmap-mean(mrnmap(find(mrnmap)));
    mrnmap=mrnmap.*bigmask;
    [mCl,l,~,~,~,~,mCl2d]=get_angular_spec(mrnmap,mrnmap,pixscale);
    [mClw1] = Cl_from_Cl2d_clip(mCl2d,pixscale,'w',(1./(mCl2dstd))');
    [mClw2] = Cl_from_Cl2d_clip(mCl2d,pixscale,'w',(1./(mCl2dave))');
    [mClw3] = Cl_from_Cl2d_clip(mCl2d,pixscale,'w',(1./(mCl2dstd+mCl2dave))');

    mCl_arr(i,:)=mCl;
    mClw1_arr(i,:)=mClw1;
    mClw2_arr(i,:)=mClw2;
    mClw3_arr(i,:)=mClw3;
    
    rnmap = readnoise_realization(darkps(8).Clf2d_ave,pixscale);
    [rCl,l,~,~,~,~,rCl2d]=get_angular_spec(rnmap,rnmap,pixscale);
    [rClw1] = Cl_from_Cl2d_clip(rCl2d,pixscale,'w',(1./(rCl2dstd))');
    [rClw2] = Cl_from_Cl2d_clip(rCl2d,pixscale,'w',(1./(rCl2dave))');
    [rClw3] = Cl_from_Cl2d_clip(rCl2d,pixscale,'w',(1./(rCl2dstd+rCl2dave))');
    rCl_arr(i,:)=rCl;
    rClw1_arr(i,:)=rClw1;
    rClw2_arr(i,:)=rClw2;    
    rClw3_arr(i,:)=rClw3;    
end
%%
figure
loglog(l,dClknox,'k');hold on
loglog(l,std(mCl_arr)./mean(mCl_arr),'r.','markersize',10);
loglog(l,std(mClw1_arr)./mean(mClw1_arr),'ro','markersize',5);
loglog(l,std(mClw2_arr)./mean(mClw2_arr),'r+','markersize',5);
loglog(l,std(mClw3_arr)./mean(mClw3_arr),'rs','markersize',5);
loglog(l,std(rCl_arr)./mean(rCl_arr),'b.','markersize',10);
loglog(l,std(rClw1_arr)./mean(rClw1_arr),'bo','markersize',5);
loglog(l,std(rClw2_arr)./mean(rClw2_arr),'b+','markersize',5);
loglog(l,std(rClw3_arr)./mean(rClw3_arr),'bs','markersize',5);
xlabel('$\ell$','interpreter','latex','fontsize',18);
ylabel('$\delta C_\ell/\bar{C_\ell}$',...
        'interpreter','latex','fontsize',18)
xlim([1e2,2e5]);ylim([1e-3,2])
legend({'Knox','Cl2d->map->mask->FFT','Cl2d->map->mask->FFT->weight(std)',...
'Cl2d->map->mask->FFT->weight(avg)','Cl2d->map->mask->FFT->weight(avg+std)',...
'Cl2d->rand->map->FFT','Cl2d->rand->map->FFT->weight(std)',...
'Cl2d->rand->map->FFT->weight(avg)','Cl2d->rand->map->FFT->weight(avg+std)'},...
'location','southwest');
savename=strcat(sprintf('%sknox_wfromsim',savedir));
print(savename,'-dpng');%close