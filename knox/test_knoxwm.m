flight=40030;
inst=1;
pixscale=7;

darkstatdir=strcat('/Users/ytcheng/ciber/doc/',...
         '20160906_NoiseRealization/darkstat/',num2str(flight),'/');
load(strcat(darkstatdir,'TM',num2str(inst),'_mdarkps'),'mdarkps');
load(strcat(darkstatdir,'TM',num2str(inst),'_darkps'),'darkps');

savedir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/test_knox/wmask/';

alldatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(alldatdir,'band',num2str(inst),'_alldat'),'alldat');
%% get knox error
[Cl,~,~,l,~,dCl,binl,dl] = Cl_from_Cl2d_clip(ones(1024),pixscale);
dCl=dCl./Cl;dCl(find(dCl~=dCl))=0;
%% Fourier mask

ell=get_l(1024,1024,pixscale,1);
logstdsm=fillpadsmooth(log10(mdarkps(8).mClf2d_std),ones(1024),5,'hsize_r',1);
fmask=zeros(1024);

figure
for i=9:29
mask=zeros(1024);mask((ell >= binl(i)) & (ell <= binl(i+1)))=1;
[x,y]=find(mask);

setwinsize(gcf,1500,300)

subplot(1,3,1)
a=logstdsm.*mask;a=a(find(a));m=median(a);
imageclip((logstdsm-m).*mask);
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
v=caxis;
title(sprintf('%.2e<l<%.2e',binl(i),binl(i+1)));

subplot(1,3,2)
cmask=ones(1024);cmask(logstdsm>median(a(:))+2*std(a(:)))=0;
cmask=cmask.*mask;
imageclip((logstdsm-m).*cmask);
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
caxis(v)

subplot(1,3,3)
imageclip(cmask);
xlim([min(x),max(x)]); ylim([min(y),max(y)]);

savename=strcat(sprintf('%smfmask/ell%d',savedir,i));
print(savename,'-dpng');%close

fmask(find(cmask))=1;

end
%%
setwinsize(gcf,1500,300)

a=log10(mdarkps(8).mClf2d_std).*fmask;m=median(a(find(a)));
subplot(1,3,1)
imageclip(log10(mdarkps(8).mClf2d_std)-m);
v=caxis;

subplot(1,3,2)
imageclip((log10(mdarkps(8).mClf2d_std)-m).*fmask);
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
caxis(v)

subplot(1,3,3)
imageclip(fmask);


savename=strcat(sprintf('%smfmask/elltot',savedir));
print(savename,'-dpng');%close
%%
weight=(1./mdarkps(8).mClf2d_std)';
mweight=(1./mdarkps(8).mClf2d_std)'.*fmask';
bigmask=alldat(8).bigmask;
%% 100 white noise sim

Cls_arr=zeros(100,29);
Clsw_arr=zeros(100,29);
Clsmw_arr=zeros(100,29);
Clsmwc5_arr=zeros(100,29);

for i=1:100
    i
    smap=randn(1024);    
    smap=smap.*bigmask;smap=smap-mean(smap(find(smap)));smap=smap.*bigmask;
    [Cl,l,~,~,~,~,Cl2d]=get_angular_spec(smap,smap,pixscale);
    Cls_arr(i,:)=Cl;

    [Clw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight);
    [Clmw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight);
    [Clmwc5] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight,'sigclip',5*ones(1,29));
    Clsw_arr(i,:)=Clw;
    Clsmw_arr(i,:)=Clmw;
    Clsmwc5_arr(i,:)=Clmwc5;
end

%% 100 noise realization PS w/o rand

Cln_arr=zeros(100,29);
Clnw_arr=zeros(100,29);
Clnmw_arr=zeros(100,29);
Clnmwc5_arr=zeros(100,29);

for i=1:100
    i
    rnmap = readnoise_realization(darkps(8).Clf2d_ave,pixscale,'norand',1);
    rnmap=rnmap.*bigmask;rnmap=rnmap-mean(rnmap(find(rnmap)));rnmap=rnmap.*bigmask;
    [Cl,l,~,~,~,~,Cl2d]=get_angular_spec(rnmap,rnmap,pixscale);
    [Clw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight);
    [Clmw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight);
    [Clnmwc5] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight,'sigclip',5*ones(1,29));
    Cln_arr(i,:)=Cl;
    Clnw_arr(i,:)=Clw;
    Clnmw_arr(i,:)=Clmw;
    Clnmwc5_arr(i,:)=Clnmwc5;
end
%%  plot dCl/Cl for 100 noise sim w/o rand
figure
loglog(l,dCl,'k');hold on

loglog(l,std(Cln_arr)./mean(Cln_arr),'r.','markersize',10);
loglog(l,std(Clnw_arr)./mean(Clnw_arr),'ro','markersize',5);
loglog(l,std(Clnmw_arr)./mean(Clnmw_arr),'r+','markersize',10);
%loglog(l,std(Clnmwc5_arr)./mean(Clnmwc5_arr),'r*');

loglog(l,std(Cls_arr)./mean(Cls_arr),'b.','markersize',10);
loglog(l,std(Clsw_arr)./mean(Clsw_arr),'bo','markersize',5);
loglog(l,std(Clsmw_arr)./mean(Clsmw_arr),'b+','markersize',10);
%loglog(l,std(Clsmwc5_arr)./mean(Clsmwc5_arr),'b*');

legend({'Knox','noise sim','noise sim weighted',...
    'noise sim weighted+Fourier mask'},...
    'Location','southwest','FontSize',15)
legend boxoff
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\delta C_\ell/\bar{C_\ell}$',...
        'interpreter','latex','fontsize',18)
xlim([1e2,2e5]);ylim([1e-3,2e0]);
savename=strcat(savedir,'knox_sim_weight');
print(savename,'-dpng');%close

%% see how weighting works mode by mode (w/o rand)
%figure
for i=9:29
figure    
mask=zeros(1024);mask((ell >= binl(i)) & (ell <= binl(i+1)))=1;
[x,y]=find(mask);

setwinsize(gcf,1500,300)

subplot(1,3,1)
a=log10(Cl2d).*mask;a=a(find(a));m1=median(a);
imageclip((log10(Cl2d)-m1).*mask);
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
v=caxis;
title(sprintf('%.2e<l<%.2e',binl(i),binl(i+1)));

subplot(1,3,2)
b=log10(Cl2d.*weight').*mask;b=b(find(b));m2=median(b);
imageclip((log10(Cl2d.*weight')-m2).*mask);
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
caxis(v)

subplot(1,3,3)
nwhist=histogram(a-m1);
hold on
whist=histogram(b-m2,'BinEdges',nwhist.BinEdges);
legend([nwhist,whist],...
{sprintf('no weight,StD=%.2f',std(a)),sprintf('weight,StD=%.2f',std(b))},...
    'location','northwest')
savename=strcat(sprintf('%sweightmode/ell%d',savedir,i));
print(savename,'-dpng');%close

end
%%%%%%%%%%%%%%%%%
%%
weight=(1./(mdarkps(8).mClf2d_std))';
mweight=(1./(mdarkps(8).mClf2d_std))'.*fmask';

weight1=(1./(mdarkps(8).mClf2d_ave))';
mweight1=(1./(mdarkps(8).mClf2d_ave))'.*fmask';

weight2=(1./(mdarkps(8).mClf2d_std.^2))';
mweight2=(1./(mdarkps(8).mClf2d_std.^2))'.*fmask';

weight3=(1./(mdarkps(8).mClf2d_ave+mdarkps(8).mClf2d_std))';
mweight3=(1./(mdarkps(8).mClf2d_ave+mdarkps(8).mClf2d_std))'.*fmask';
%% 100 white noise sim w/ chi2 rand

Cls_arr=zeros(100,29);
Clsw_arr=zeros(100,29);
Clsmw_arr=zeros(100,29);
Clsw1_arr=zeros(100,29);
Clsmw1_arr=zeros(100,29);
Clsw2_arr=zeros(100,29);
Clsmw2_arr=zeros(100,29);
Clsw3_arr=zeros(100,29);
Clsmw3_arr=zeros(100,29);

for i=1:100
    i  
    
    Cl2din=chi2rnd(ones(1024, 1024).*2)./2;
    smap=ifft2(fftshift(sqrt(Cl2din)).*fft2(normrnd(0,1,1024)));
    smap=smap'./(pixscale/3600.0*pi/180.0);
    smap=real(smap)./abs(real(smap)).*abs(smap);
    smap=smap.*bigmask;smap=smap-mean(smap(find(smap)));smap=smap.*bigmask;
    
    [Cl,l,~,~,~,~,Cl2d]=get_angular_spec(smap,smap,pixscale);
    Cls_arr(i,:)=Cl;

    [Clw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight);
    [Clmw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight);
    Clsw_arr(i,:)=Clw;
    Clsmw_arr(i,:)=Clmw;
    
    [Clw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight1);
    [Clmw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight1);
    Clsw1_arr(i,:)=Clw;
    Clsmw1_arr(i,:)=Clmw;
    
    [Clw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight2);
    [Clmw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight2);
    Clsw2_arr(i,:)=Clw;
    Clsmw2_arr(i,:)=Clmw;
    
    [Clw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight3);
    [Clmw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight3);
    Clsw3_arr(i,:)=Clw;
    Clsmw3_arr(i,:)=Clmw;
end
%% 100 noise realization PS w/ rand

Cln_arr=zeros(100,29);
Clnw_arr=zeros(100,29);
Clnmw_arr=zeros(100,29);
Clnw1_arr=zeros(100,29);
Clnmw1_arr=zeros(100,29);
Clnw2_arr=zeros(100,29);
Clnmw2_arr=zeros(100,29);
Clnw3_arr=zeros(100,29);
Clnmw3_arr=zeros(100,29);

for i=1:100
    i
    rnmap = readnoise_realization(darkps(8).Clf2d_ave,pixscale);
    rnmap=rnmap.*bigmask;rnmap=rnmap-mean(rnmap(find(rnmap)));rnmap=rnmap.*bigmask;

    [Cl,l,~,~,~,~,Cl2d]=get_angular_spec(rnmap,rnmap,pixscale);
    Cln_arr(i,:)=Cl;
    
    [Clw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight);
    [Clmw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight);
    Clnw_arr(i,:)=Clw;
    Clnmw_arr(i,:)=Clmw;
    
    [Clw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight1);
    [Clmw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight1);
    Clnw1_arr(i,:)=Clw;
    Clnmw1_arr(i,:)=Clmw;

    [Clw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight2);
    [Clmw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight2);
    Clnw2_arr(i,:)=Clw;
    Clnmw2_arr(i,:)=Clmw;

    [Clw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight3);
    [Clmw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight3);
    Clnw3_arr(i,:)=Clw;
    Clnmw3_arr(i,:)=Clmw;
    
end

%%  plot dCl/Cl for 100 noise sim w/ rand
figure
loglog(l,dCl,'k');hold on

loglog(l,std(Cln_arr)./mean(Cln_arr),'r.','markersize',10);
loglog(l,std(Clnw_arr)./mean(Clnw_arr),'ro','markersize',5);
loglog(l,std(Clnmw_arr)./mean(Clnmw_arr),'r+','markersize',10);

loglog(l,std(Cls_arr)./mean(Cls_arr),'b.','markersize',10);
loglog(l,std(Clsw_arr)./mean(Clsw_arr),'bo','markersize',5);
loglog(l,std(Clsmw_arr)./mean(Clsmw_arr),'b+','markersize',10);

legend({'Knox','noise sim','noise sim weighted',...
    'noise sim weighted+Fourier mask'},...
    'Location','southwest','FontSize',15)
legend boxoff
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\delta C_\ell/\bar{C_\ell}$',...
        'interpreter','latex','fontsize',18)
xlim([1e2,2e5]);ylim([1e-3,2e0]);
savename=strcat(savedir,'knox_sim_weight0_rand');
print(savename,'-dpng');%close
%%%%%%%
figure
loglog(l,dCl,'k');hold on

loglog(l,std(Cln_arr)./mean(Cln_arr),'r.','markersize',10);
loglog(l,std(Clnw1_arr)./mean(Clnw1_arr),'ro','markersize',5);
loglog(l,std(Clnmw1_arr)./mean(Clnmw1_arr),'r+','markersize',10);

loglog(l,std(Cls_arr)./mean(Cls_arr),'b.','markersize',10);
loglog(l,std(Clsw1_arr)./mean(Clsw1_arr),'bo','markersize',5);
loglog(l,std(Clsmw1_arr)./mean(Clsmw1_arr),'b+','markersize',10);

legend({'Knox','noise sim','noise sim weighted',...
    'noise sim weighted+Fourier mask'},...
    'Location','southwest','FontSize',15)
legend boxoff
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\delta C_\ell/\bar{C_\ell}$',...
        'interpreter','latex','fontsize',18)
xlim([1e2,2e5]);ylim([1e-3,2e0]);
savename=strcat(savedir,'knox_sim_weight1_rand');
print(savename,'-dpng');%close
%%%%%%%
figure
loglog(l,dCl,'k');hold on

loglog(l,std(Cln_arr)./mean(Cln_arr),'r.','markersize',10);
loglog(l,std(Clnw2_arr)./mean(Clnw2_arr),'ro','markersize',5);
loglog(l,std(Clnmw2_arr)./mean(Clnmw2_arr),'r+','markersize',10);
%loglog(l,std(Clnmwc51_arr)./mean(Clnmwc51_arr),'r*');

loglog(l,std(Cls_arr)./mean(Cls_arr),'b.','markersize',10);
loglog(l,std(Clsw2_arr)./mean(Clsw2_arr),'bo','markersize',5);
loglog(l,std(Clsmw2_arr)./mean(Clsmw2_arr),'b+','markersize',10);
%loglog(l,std(Clsmwc51_arr)./mean(Clsmwc51_arr),'b*');

legend({'Knox','noise sim','noise sim weighted',...
    'noise sim weighted+Fourier mask'},...
    'Location','southwest','FontSize',15)
legend boxoff
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\delta C_\ell/\bar{C_\ell}$',...
        'interpreter','latex','fontsize',18)
xlim([1e2,2e5]);ylim([1e-3,2e0]);
savename=strcat(savedir,'knox_sim_weight2_rand');
print(savename,'-dpng');%close

%%%%%%%
figure
loglog(l,dCl,'k');hold on

loglog(l,std(Cln_arr)./mean(Cln_arr),'r.','markersize',10);
loglog(l,std(Clnw3_arr)./mean(Clnw3_arr),'ro','markersize',5);
loglog(l,std(Clnmw3_arr)./mean(Clnmw3_arr),'r+','markersize',10);

loglog(l,std(Cls_arr)./mean(Cls_arr),'b.','markersize',10);
loglog(l,std(Clsw3_arr)./mean(Clsw3_arr),'bo','markersize',5);
loglog(l,std(Clsmw3_arr)./mean(Clsmw3_arr),'b+','markersize',10);

legend({'Knox','noise sim','noise sim weighted',...
    'noise sim weighted+Fourier mask'},...
    'Location','southwest','FontSize',15)
legend boxoff
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\delta C_\ell/\bar{C_\ell}$',...
        'interpreter','latex','fontsize',18)
xlim([1e2,2e5]);ylim([1e-3,2e0]);
savename=strcat(savedir,'knox_sim_weight3_rand');
print(savename,'-dpng');%close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% real RN data 
Cl2d_arr=mdarkps(8).mCld2d_arr;
Cl_arr=mdarkps(8).mCld_arr;

Clw_arr=zeros(size(Cl2d_arr,1),29);
Clmw_arr=zeros(size(Cl2d_arr,1),29);
%Clmwc5_arr=zeros(size(Cl2d_arr,1),29);
Clw1_arr=zeros(size(Cl2d_arr,1),29);
Clmw1_arr=zeros(size(Cl2d_arr,1),29);
%Clmwc51_arr=zeros(size(Cl2d_arr,1),29);
Clw2_arr=zeros(size(Cl2d_arr,1),29);
Clmw2_arr=zeros(size(Cl2d_arr,1),29);
Clw3_arr=zeros(size(Cl2d_arr,1),29);
Clmw3_arr=zeros(size(Cl2d_arr,1),29);

for i=1:size(Cl2d_arr,1)
  i  
Cl2d=squeeze(Cl2d_arr(i,:,:));
[Clw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight);
[Clmw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight);
%[Clmwc5] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight,'sigclip',5*ones(1,29));
Clw_arr(i,:)=Clw;
Clmw_arr(i,:)=Clmw;
%Clmwc5_arr(i,:)=Clmwc5;

[Clw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight1);
[Clmw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight1);
%[Clmwc5] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight1,'sigclip',5*ones(1,29));
Clw1_arr(i,:)=Clw;
Clmw1_arr(i,:)=Clmw;
%Clmwc51_arr(i,:)=Clmwc5;

[Clw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight2);
[Clmw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight2);
Clw2_arr(i,:)=Clw;
Clmw2_arr(i,:)=Clmw;

[Clw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',weight3);
[Clmw] = Cl_from_Cl2d_clip(Cl2d,pixscale,'w',mweight3);
Clw3_arr(i,:)=Clw;
Clmw3_arr(i,:)=Clmw;
end
%%  plot dCl/Cl real RN data

figure
loglog(l,dCl,'k');hold on

loglog(l,std(Cl_arr)./mean(Cl_arr),'r.','markersize',10);
loglog(l,std(Clw_arr)./mean(Clw_arr),'ro','markersize',5);
loglog(l,std(Clmw_arr)./mean(Clmw_arr),'r+','markersize',10);
%loglog(l,std(Clnmwc5_arr)./mean(Clnmwc5_arr),'r*');

loglog(l,std(Cls_arr)./mean(Cls_arr),'b.','markersize',10);
loglog(l,std(Clsw_arr)./mean(Clsw_arr),'bo','markersize',5);
loglog(l,std(Clsmw_arr)./mean(Clsmw_arr),'b+','markersize',10);
%loglog(l,std(Clsmwc5_arr)./mean(Clsmwc5_arr),'b*');

legend({'Knox','noise sim','noise sim weighted',...
    'noise sim weighted+Fourier mask'},...
    'Location','southwest','FontSize',15)
legend boxoff
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\delta C_\ell/\bar{C_\ell}$',...
        'interpreter','latex','fontsize',18)
xlim([1e2,2e5]);ylim([1e-3,2e0]);
savename=strcat(savedir,'knox_dat_weight0');
print(savename,'-dpng');%close
%%%%%%%
figure
loglog(l,dCl,'k');hold on

loglog(l,std(Cl_arr)./mean(Cl_arr),'r.','markersize',10);
loglog(l,std(Clw1_arr)./mean(Clw1_arr),'ro','markersize',5);
loglog(l,std(Clmw1_arr)./mean(Clmw1_arr),'r+','markersize',10);
%loglog(l,std(Clnmwc51_arr)./mean(Clnmwc51_arr),'r*');

loglog(l,std(Cls_arr)./mean(Cls_arr),'b.','markersize',10);
loglog(l,std(Clsw1_arr)./mean(Clsw1_arr),'bo','markersize',5);
loglog(l,std(Clsmw1_arr)./mean(Clsmw1_arr),'b+','markersize',10);
%loglog(l,std(Clsmwc51_arr)./mean(Clsmwc51_arr),'b*');

legend({'Knox','noise sim','noise sim weighted',...
    'noise sim weighted+Fourier mask'},...
    'Location','southwest','FontSize',15)
legend boxoff
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\delta C_\ell/\bar{C_\ell}$',...
        'interpreter','latex','fontsize',18)
xlim([1e2,2e5]);ylim([1e-3,2e0]);
savename=strcat(savedir,'knox_dat_weight1');
print(savename,'-dpng');%close
%%%%%%%
figure
loglog(l,dCl,'k');hold on

loglog(l,std(Cl_arr)./mean(Cl_arr),'r.','markersize',10);
loglog(l,std(Clw2_arr)./mean(Clw2_arr),'ro','markersize',5);
loglog(l,std(Clmw2_arr)./mean(Clmw2_arr),'r+','markersize',10);

loglog(l,std(Cls_arr)./mean(Cls_arr),'b.','markersize',10);
loglog(l,std(Clsw2_arr)./mean(Clsw2_arr),'bo','markersize',5);
loglog(l,std(Clsmw2_arr)./mean(Clsmw2_arr),'b+','markersize',10);

legend({'Knox','noise sim','noise sim weighted',...
    'noise sim weighted+Fourier mask'},...
    'Location','southwest','FontSize',15)
legend boxoff
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\delta C_\ell/\bar{C_\ell}$',...
        'interpreter','latex','fontsize',18)
xlim([1e2,2e5]);ylim([1e-3,2e0]);
savename=strcat(savedir,'knox_dat_weight2');
print(savename,'-dpng');%close
%%%%%%%
figure
loglog(l,dCl,'k');hold on

loglog(l,std(Cl_arr)./mean(Cl_arr),'r.','markersize',10);
loglog(l,std(Clw3_arr)./mean(Clw3_arr),'ro','markersize',5);
loglog(l,std(Clmw3_arr)./mean(Clmw3_arr),'r+','markersize',10);

loglog(l,std(Cls_arr)./mean(Cls_arr),'b.','markersize',10);
loglog(l,std(Clsw3_arr)./mean(Clsw3_arr),'bo','markersize',5);
loglog(l,std(Clsmw3_arr)./mean(Clsmw3_arr),'b+','markersize',10);

legend({'Knox','noise sim','noise sim weighted',...
    'noise sim weighted+Fourier mask'},...
    'Location','southwest','FontSize',15)
legend boxoff
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\delta C_\ell/\bar{C_\ell}$',...
        'interpreter','latex','fontsize',18)
xlim([1e2,2e5]);ylim([1e-3,2e0]);
savename=strcat(savedir,'knox_dat_weight3');
print(savename,'-dpng');%close
