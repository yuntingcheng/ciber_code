%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See if the mask increase the 2DCl variance,
% so weighting can't work well (like random amp case).
% But the variance doesn't seem to increse...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
pixscale=7;

darkstatdir=strcat('/Users/ytcheng/ciber/doc/',...
         '20160906_NoiseRealization/darkstat/',num2str(flight),'/');
load(strcat(darkstatdir,'TM',num2str(inst),'_mdarkps'),'mdarkps');
load(strcat(darkstatdir,'TM',num2str(inst),'_darkps'),'darkps');

savedir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/test_knox/maskvar/';

alldatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(alldatdir,'band',num2str(inst),'_alldat'),'alldat');
%%
bigmask=alldat(8).bigmask;
weight=(1./darkps(8).Clf2d_std)';
weightm=(1./mdarkps(8).mClf2d_std)';
Cl_arr=zeros(100,29);
mCl_arr=zeros(100,29);
rCl_arr=zeros(100,29);
m1Cl_arr=zeros(100,29);
Clw_arr=zeros(100,29);
mClw_arr=zeros(100,29);
rClw_arr=zeros(100,29);
m1Clw_arr=zeros(100,29);

Cl2d_arr=zeros(100,1024,1024);
mCl2d_arr=zeros(100,1024,1024);
rCl2d_arr=zeros(100,1024,1024);
m1Cl2d_arr=zeros(100,1024,1024);
for i=1:100
    i
    rnmap = readnoise_realization(darkps(8).Clf2d_ave,pixscale,'norand',1);
    mrnmap=rnmap.*bigmask;mrnmap=mrnmap-mean(mrnmap(find(mrnmap)));
    mrnmap=mrnmap.*bigmask;
    rrnmap= readnoise_realization(darkps(8).Clf2d_ave,pixscale,'norand',0);
    m1rnmap=readnoise_realization(mdarkps(8).mClf2d_ave,pixscale,'norand',1);
    [Cl,l,~,~,~,~,Cl2d]=get_angular_spec(rnmap,rnmap,pixscale);
    [mCl,l,~,~,~,~,mCl2d]=get_angular_spec(mrnmap,mrnmap,pixscale);
    [rCl,l,~,~,~,~,rCl2d]=get_angular_spec(rrnmap,rrnmap,pixscale);
    [m1Cl,l,~,~,~,~,m1Cl2d]=get_angular_spec(m1rnmap,m1rnmap,pixscale);
    [Clw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight);
    [mClw] = Cl_from_Cl2d_clip(mCl2d,pixscale,'w',weightm);
    [rClw] = Cl_from_Cl2d_clip(rCl2d,pixscale,'w',weight);
    [m1Clw] = Cl_from_Cl2d_clip(m1Cl2d,pixscale,'w',weightm);
    Cl_arr(i,:)=Cl;
    mCl_arr(i,:)=mCl;
    rCl_arr(i,:)=rCl;
    m1Cl_arr(i,:)=m1Cl;
    
    Clw_arr(i,:)=Clw;
    mClw_arr(i,:)=mClw;
    rClw_arr(i,:)=rClw;
    m1Clw_arr(i,:)=m1Clw;
    
    Cl2d_arr(i,:,:)=Cl2d;
    mCl2d_arr(i,:,:)=mCl2d;
    rCl2d_arr(i,:,:)=rCl2d;
    m1Cl2d_arr(i,:,:)=m1Cl2d;
end

%% get knox error
[Cl,~,~,l,~,dCl,binl,dl] = Cl_from_Cl2d_clip(ones(1024),pixscale);;
dCl=dCl./Cl;dCl(find(dCl~=dCl))=0;
%%
p0=loglog(l,dCl,'k');hold on
p1=loglog(l,std(Cl_arr)./mean(Cl_arr),'k.','markersize',10);
loglog(l,std(Clw_arr)./mean(Clw_arr),'ko','markersize',5);
p2=loglog(l,std(rCl_arr)./mean(rCl_arr),'b.','markersize',10);
loglog(l,std(rClw_arr)./mean(rClw_arr),'bo','markersize',5);
p3=loglog(l,std(mCl_arr)./mean(mCl_arr),'g.','markersize',10);
loglog(l,std(mClw_arr)./mean(mClw_arr),'go','markersize',5);
p4=loglog(l,std(m1Cl_arr)./mean(m1Cl_arr),'r.','markersize',10);
loglog(l,std(m1Clw_arr)./mean(m1Clw_arr),'ro','markersize',5);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\delta C_\ell/\bar{C_\ell}$',...
        'interpreter','latex','fontsize',18)
xlim([1e2,2e5]);ylim([1e-3,2])
legend([p0,p1,p2,p3,p4],{'Knox','Cl2d->map->FFT',...
    'Cl2d->rand->map->FFT','Cl2d->map->mask->FFT',...
    'maskCl2d->map->FFT'},'location','southwest');
savename=strcat(sprintf('%sknox',savedir));
print(savename,'-dpng');%close

%%
d0=squeeze(std(Cl2d_arr)./mean(Cl2d_arr));clear Cl2d_arr
dm=squeeze(std(mCl2d_arr)./mean(mCl2d_arr));clear mCl2d_arr
dr=squeeze(std(rCl2d_arr)./mean(rCl2d_arr));clear rCl2d_arr
dm1=squeeze(std(m1Cl2d_arr)./mean(m1Cl2d_arr));clear m1Cl2d_arr
%%
ell=get_l(1024,1024,pixscale,1);

for i=9:29
figure    
mask=zeros(1024);mask((ell >= binl(i)) & (ell <= binl(i+1)))=1;
invmask=ones(1024);invmask(find(mask))=0;
[x,y]=find(mask);

setwinsize(gcf,1200,1200)

subplot(2,2,1)
a=dm1.*mask;a=a(find(a));
imageclip(dm1.*mask);
v=caxis;
imageclip(dm1.*mask+invmask);
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
caxis(v)
title({sprintf('%.2e<l<%.2e',binl(i),binl(i+1)),'mask2DCl'});

subplot(2,2,2)
b=dm.*mask;b=b(find(b));
imageclip(dm.*mask+invmask);
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
caxis(v)
title({sprintf('%.2e<l<%.2e',binl(i),binl(i+1)),'unmask2DCl->mask'});

subplot(2,2,3)
c=dr.*mask;c=c(find(c));
imageclip(dr.*mask+invmask.*mean(c));
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
caxis(v)
title({sprintf('%.2e<l<%.2e',binl(i),binl(i+1)),'unmask2DCl->rand->mask'});

subplot(2,2,4)
h1=histogram(a);hold on
h2=histogram(b,'BinEdges',h1.BinEdges);
h3=histogram(c,'BinEdges',h1.BinEdges);
legend([h1,h2,h3],...
{sprintf('mask2DCl,avg=%.2f,med=%.2f',mean(a),median(a)),...
sprintf('unmask2DCl,avg=%.2f,med=%.2f',mean(b),median(b)),...
sprintf('unmask2DClrand,avg=%.2f,med=%.2f',mean(c),median(c))},...
   'location','northeast')
savename=strcat(sprintf('%smodes/ell%d',savedir,i));
print(savename,'-dpng');close

end
%%%%%%%%%%%%%%%%%
