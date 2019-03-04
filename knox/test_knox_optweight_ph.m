%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Similar to test_knoxm2.m, but now adding ph noise 
% in total noise model, and run the simulated signal
% to see the var with different noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
pixscale=7;

darkstatdir=strcat('/Users/ytcheng/ciber/doc/',...
         '20160906_NoiseRealization/darkstat/',num2str(flight),'/');
load(strcat(darkstatdir,'TM',num2str(inst),'_mdarkps'),'mdarkps');
load(strcat(darkstatdir,'TM',num2str(inst),'_darkps'),'darkps');

savedir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/test_knox/optweight1/';

alldatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(alldatdir,'band',num2str(inst),'_alldat'),'alldat');
bigmask=alldat(8).bigmask;
%% get knox error
[Cl,~,~,l,~,dCl,binl,dl] = Cl_from_Cl2d_clip(ones(1024),pixscale);
dCl=dCl./Cl;dCl(find(dCl~=dCl))=0;
%% 
Cl2dave=mdarkps(8).mClf2d_ave;
Cl2dstd1=mdarkps(8).mClf2d_std;
Cl2dstd2=mdarkps(8).mClf2d_std.^2;
Cl2dstd4=mdarkps(8).mClf2d_std.^4;
Cl2dstd6=mdarkps(8).mClf2d_std.^6;
%% normalize Cl2dave & Cl2dvar
ell=get_l(1024,1024,pixscale,1);

for i=9:29
mask=zeros(1024);mask((ell >= binl(i)) & (ell <= binl(i+1)))=1;
[x,y]=find(mask);

a=Cl2dave.*mask;a=a(find(a));
Cl2dave(find(mask))=Cl2dave(find(mask))./median(a);
a=Cl2dstd1.*mask;a=a(find(a));
Cl2dstd1(find(mask))=Cl2dstd1(find(mask))./median(a);
a=Cl2dstd2.*mask;a=a(find(a));
Cl2dstd2(find(mask))=Cl2dstd2(find(mask))./median(a);
a=Cl2dstd4.*mask;a=a(find(a));
Cl2dstd4(find(mask))=Cl2dstd4(find(mask))./median(a);
a=Cl2dstd6.*mask;a=a(find(a));
Cl2dstd6(find(mask))=Cl2dstd6(find(mask))./median(a);
end
%% get the weights
weight0=1./Cl2dave';
weight2=1./Cl2dstd2';
weight4=1./Cl2dstd4';
weight6=1./Cl2dstd6';

weight22=1./(Cl2dave.*0.5+Cl2dstd2.*0.5)';
weight44=1./(Cl2dave.*0.5+Cl2dstd4.*0.5)';
weight66=1./(Cl2dave.*0.5+Cl2dstd6.*0.5)';
%% 100 noise realization PS w/ rand
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
frate=cp(inst).framerate;
G1=-2.68;G2=cal./G1;
[~,l,~,~,lbin]=get_angular_spec(randn(1024),randn(1024),pixscale);
logx=[log10(2e3) log10(1e2)];x=10.^logx;
logy=[log10(2) log10(12)];y=10.^logy;
logCl=interp1(logx,logy,log10(l),'linear','extrap');
Clshot=ones(size(logCl)).*1e3.*2.*pi./1e5./(1e5+1);
Clin=(((10.^logCl)).*2.*pi./l./(l+1))+Clshot;
zldir = '/Users/ytcheng/ciber/data/TM1_zl_pmk/';
zlname = strcat(zldir,'TM1_40030_zlgrad',alldat(8).name,'.fits');
zlgrad = fliplr(fitsread_orient(zlname));
drdir=strcat('/Users/ytcheng/ciber/data/',num2str(flight),'/dr/');
dataname = strcat(drdir,'TM',num2str(inst),'_',alldat(1).name,...
    '_dr150206.mat');
load(dataname);labflat=data.flat.labflat;

Cln_arr=zeros(100,29);
Clnw0_arr=zeros(100,29);
Clnw2_arr=zeros(100,29);
Clnw2m_arr=zeros(100,29);
Clnw4_arr=zeros(100,29);
Clnw6_arr=zeros(100,29);

Clnw22_arr=zeros(100,29);
Clnw44_arr=zeros(100,29);
Clnw66_arr=zeros(100,29);

Cls_arr=zeros(100,29);
Clsw0_arr=zeros(100,29);
Clsw2_arr=zeros(100,29);
Clsw4_arr=zeros(100,29);
Clsw6_arr=zeros(100,29);

Clsw22_arr=zeros(100,29);
Clsw44_arr=zeros(100,29);
Clsw66_arr=zeros(100,29);

for i=1:100
    i
    rnmap = readnoise_realization(darkps(8).Clf2d_ave,pixscale);
    sigmap = map_from_power_spec(lbin,Clin,1024,1024,pixscale,1);
    slopemap=(sigmap+zlgrad+rnmap).*labflat./cal;
    [phmap]=photonnoise_realization(slopemap,G1,alldat(8).nfr,frate);
    nmap=rnmap+phmap;
    nmap=nmap.*bigmask;nmap=nmap-mean(nmap(find(nmap)));nmap=nmap.*bigmask;

    [Cl,l,~,~,~,~,Cl2d]=get_angular_spec(nmap,nmap,pixscale);
    Cln_arr(i,:)=Cl;
    
    Clnw0_arr(i,:) = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight0);
    Clnw2_arr(i,:) = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight2);
    Clnw2m_arr(i,:) = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight2.*fmask');
    Clnw4_arr(i,:) = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight4);
    Clnw6_arr(i,:) = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight6);
    
    Clnw22_arr(i,:) = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight22);
    Clnw44_arr(i,:) = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight44);
    Clnw66_arr(i,:) = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight66);
    
    
    smap=sigmap.*bigmask;smap=smap-mean(smap(find(smap)));smap=smap.*bigmask;

    [Cl,l,~,~,~,~,Cl2d]=get_angular_spec(smap,smap,pixscale);
    Cls_arr(i,:)=Cl;
    
    Clsw0_arr(i,:) = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight0);
    Clsw2_arr(i,:) = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight2);
    Clsw4_arr(i,:) = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight4);
    Clsw6_arr(i,:) = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight6);
    
    Clsw22_arr(i,:) = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight22);
    Clsw44_arr(i,:) = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight44);
    Clsw66_arr(i,:) = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight66);

end
%%
figure
loglog(l,dCl,'k');hold on

pltn=loglog(l,std(Cln_arr)./mean(Cln_arr),'k.','markersize',10);
pltm=loglog(l,std(Clnw0_arr)./mean(Clnw0_arr),'ro','markersize',5);
pltv=loglog(l,std(Clnw2_arr)./mean(Clnw2_arr),'bo','markersize',5);
loglog(l,std(Clnw22_arr)./mean(Clnw22_arr),'b+','markersize',5);
loglog(l,std(Clnw2m_arr)./mean(Clnw2m_arr),'bx','markersize',5);
pltv2=loglog(l,std(Clnw4_arr)./mean(Clnw4_arr),'mo','markersize',5);
loglog(l,std(Clnw44_arr)./mean(Clnw44_arr),'m+','markersize',5);
pltv3=loglog(l,std(Clnw6_arr)./mean(Clnw66_arr),'go','markersize',5);
loglog(l,std(Clnw66_arr)./mean(Clnw66_arr),'g+','markersize',5);

legend([pltn,pltm,pltv,pltv2,pltv3],...
    {'no weight','mean','var','var^2','var^3'},...
    'Location','southwest','FontSize',15)
legend boxoff
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\delta C_\ell/\bar{C_\ell}$',...
        'interpreter','latex','fontsize',18)
xlim([1e2,2e5]);%ylim([3e-2,2e0]);
savename=strcat(savedir,'knox_sim_rand_CldCl');
print(savename,'-dpng');%close
%%
figure
pltn=loglog(l,mean(Cln_arr),'ko-');hold on
pltm=loglog(l,mean(Clnw0_arr),'ro-');
pltv=loglog(l,mean(Clnw2_arr),'bo-');
loglog(l,mean(Clnw22_arr),'b+-');
loglog(l,mean(Clnw2m_arr),'bx--');
pltv2=loglog(l,mean(Clnw4_arr),'mo-');
loglog(l,mean(Clnw44_arr),'m+-');
pltv3=loglog(l,mean(Clnw6_arr),'go-');
loglog(l,mean(Clnw66_arr),'g+-');

legend([pltn,pltm,pltv,pltv2,pltv3],...
    {'no weight','mean','var','var^2','var^3'},...
    'Location','southwest','FontSize',15)
legend boxoff

xlim([1e2,2e5]);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$C_\ell(nW m^{-2} sr^{-1})$',...
        'interpreter','latex','fontsize',18)
savename=strcat(savedir,'knox_sim_rand_Cl');
print(savename,'-dpng');%close

figure
pltn=loglog(l,std(Cln_arr),'ko-');hold on
pltm=loglog(l,std(Clnw0_arr),'ro-');
pltv=loglog(l,std(Clnw2_arr),'bo-');
loglog(l,std(Clnw22_arr),'b+-');
loglog(l,std(Clnw2m_arr),'bx--');
pltv2=loglog(l,std(Clnw4_arr),'mo-');
loglog(l,std(Clnw44_arr),'m+-');
pltv4=loglog(l,std(Clnw6_arr),'go-');
loglog(l,std(Clnw66_arr),'g+-');
xlim([1e2,2e5]);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\delta C_\ell(nW m^{-2} sr^{-1})$',...
        'interpreter','latex','fontsize',18)
legend([pltn,pltm,pltv,pltv2,pltv3],...
    {'no weight','mean','var','var^2','var^3'},...
    'Location','southwest','FontSize',15)
legend boxoff

savename=strcat(savedir,'knox_sim_rand_dCl');
print(savename,'-dpng');%close
%%

figure
pltn=loglog(l,mean(Cls_arr),'ko-');hold on
pltm=loglog(l,mean(Clsw0_arr),'ro-');
pltv=loglog(l,mean(Clsw2_arr),'bo-');
loglog(l,mean(Clsw22_arr),'b+-');
pltv2=loglog(l,mean(Clsw4_arr),'mo-');
loglog(l,mean(Clsw44_arr),'m+-');
pltv3=loglog(l,mean(Clsw6_arr),'go-');
loglog(l,mean(Clsw66_arr),'g+-');

legend([pltn,pltm,pltv,pltv2,pltv3],...
    {'no weight','mean','var','var^2','var^3'},...
    'Location','southwest','FontSize',15)
legend boxoff

xlim([1e2,2e5]);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$C_\ell(nW m^{-2} sr^{-1})$',...
        'interpreter','latex','fontsize',18)
savename=strcat(savedir,'knox_sim_rand_Cl_sig');
print(savename,'-dpng');%close

figure
pltn=loglog(l,std(Cls_arr),'ko-');hold on
pltm=loglog(l,std(Clsw0_arr),'ro-');
pltv=loglog(l,std(Clsw2_arr),'bo-');
loglog(l,std(Clsw22_arr),'b+-');
pltv2=loglog(l,std(Clsw4_arr),'mo-');
loglog(l,std(Clsw44_arr),'m+-');
pltv4=loglog(l,std(Clsw6_arr),'go-');
loglog(l,std(Clsw66_arr),'g+-');
xlim([1e2,2e5]);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\delta C_\ell(nW m^{-2} sr^{-1})$',...
        'interpreter','latex','fontsize',18)
legend([pltn,pltm,pltv,pltv2,pltv3],...
    {'no weight','mean','var','var^2','var^3'},...
    'Location','southwest','FontSize',15)
legend boxoff

savename=strcat(savedir,'knox_sim_rand_dCl_sig');
print(savename,'-dpng');%close


