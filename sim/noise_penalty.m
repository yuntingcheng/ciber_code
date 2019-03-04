flight=40030;
inst=1;
pixscale=7;

cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
frate=cp(inst).framerate;
G1=-2.616;G2=cal./G1;

DCdir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/40030/';
load(strcat(DCdir,'band',num2str(inst),'_DCtemplate'),'DCtemplate');

darkdir='/Users/ytcheng/ciber/doc/20150810_DarkCurrent/';
fields=get_fields(flight,inst);

alldatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(alldatdir,'band',num2str(inst),'_alldat'),'alldat');

savedir='/Users/ytcheng/ciber/doc/20170130_HalfNoise/';
ell = get_l(1024,1024,pixscale,1);
[~,~,~,l,~,dClknox,binl,dl,~,~,nell] = Cl_from_Cl2d_clip(ones(1024),pixscale);

%% get RN+ph Cl2d full

for ifield=1:numel(fields)
fname=fields(ifield).name;
fulldir=sprintf('%s/%d/TM%d/%s/full/data/',darkdir,flight,inst,fname);
scanfull=dir(strcat(fulldir,'dark*'));
ndarks=numel(scanfull);

nfr=fields(ifield).nfr;
bigmask=alldat(ifield).bigmask;
rawmap=alldat(ifield).rawmap;

Cl2dfull_arr=zeros(ndarks,1024,1024);
sig_arr=zeros(1,ndarks);
for i=1:ndarks
    load(strcat(fulldir,scanfull(i).name));
    rnmap=darkmap-DCtemplate;
    [~,mask]=get_skymap(rnmap,bigmask,3);
    rnmap=dc_offset_remove(rnmap,mask);
    rnmap=rnmap.*mask;
    rnmap=rnmap-mean(rnmap(find(rnmap)));rnmap=rnmap.*mask;
    
    phmap=photonnoise_realization(rawmap,G1,nfr,frate);
    phmap=phmap.*mask;
    phmap=phmap-mean(phmap(find(phmap)));phmap=phmap.*mask;
    
    nmap=rnmap+phmap;
    sig_arr(i)=std(nmap(find(nmap)));
    [~,~,~,~,~,~,Cl2d]=get_angular_spec(nmap,nmap,pixscale);
    Cl2dfull_arr(i,:,:)=Cl2d;
end
Cldata(ifield).Cl2dfull_arr=Cl2dfull_arr;
Cldata(ifield).sigfull_arr=sig_arr;
end
%% get RN+ph Cl2d half

for ifield=1:numel(fields)
fname=fields(ifield).name;
dir1=sprintf('%s/%d/TM%d/%s/first/data/',darkdir,flight,inst,fname);
dir2=sprintf('%s/%d/TM%d/%s/second/data/',darkdir,flight,inst,fname);
scan1=dir(strcat(dir1,'dark*'));
scan2=dir(strcat(dir2,'dark*'));
ndarks=numel(scan1);

nfr=fields(ifield).nfrhalf;
bigmask=alldat(ifield).bigmask;
rawmap=alldat(ifield).rawmap;

Cl2ddiff_arr=zeros(ndarks,1024,1024);
sig_arr=zeros(1,ndarks);
for i=1:ndarks
    load(strcat(dir1,scan1(i).name));darkmap1=darkmap;
    load(strcat(dir2,scan2(i).name));darkmap2=darkmap;
    rnmap=(darkmap1-darkmap2)./2;
    [~,mask]=get_skymap(rnmap,bigmask,3);
    rnmap=dc_offset_remove(rnmap,mask);
    rnmap=rnmap.*mask;
    rnmap=rnmap-mean(rnmap(find(rnmap)));rnmap=rnmap.*mask;
    
    phmap1=photonnoise_realization(rawmap,G1,nfr,frate);
    phmap2=photonnoise_realization(rawmap,G1,nfr,frate);
    phmap=(phmap1-phmap2)./2;
    phmap=phmap.*mask;
    phmap=phmap-mean(phmap(find(phmap)));phmap=phmap.*mask;
    
    nmap=rnmap+phmap;
    sig_arr(i)=std(nmap(find(nmap)));
    [~,~,~,~,~,~,Cl2d]=get_angular_spec(nmap,nmap,pixscale);
    Cl2ddiff_arr(i,:,:)=Cl2d;
end
Cldata(ifield).Cl2ddiff_arr=Cl2ddiff_arr;
Cldata(ifield).sigdiff_arr=sig_arr;
end
%% 1D Cl mean and var
for ifield=1:numel(fields)
Cl2dfull_arr=Cldata(ifield).Cl2dfull_arr;
Cl2ddiff_arr=Cldata(ifield).Cl2ddiff_arr;
weightf=1./(squeeze(std(Cl2dfull_arr)));
weightd=1./(squeeze(std(Cl2ddiff_arr)));

Clf_arr=zeros(size(Cl2dfull_arr,1),29);
Cld_arr=zeros(size(Cl2ddiff_arr,1),29);
for i=1:size(Cl2dfull_arr,1)
    Cl2df=squeeze(Cl2dfull_arr(i,:,:));
    Cl2dd=squeeze(Cl2ddiff_arr(i,:,:));
    [Clf,~,~,l] = Cl_from_Cl2d_clip(Cl2df,pixscale,'w',weightf);
    [Cld] = Cl_from_Cl2d_clip(Cl2dd,pixscale,'w',weightd);
    
    Clf_arr(i,:)=Clf;
    Cld_arr(i,:)=Cld;
end
Cldata(ifield).weightf=weightf;
Cldata(ifield).weightd=weightd;
Cldata(ifield).Clf_arr=Clf_arr;
Cldata(ifield).Cld_arr=Cld_arr;
Cldata(ifield).Clfavg=mean(Clf_arr);
Cldata(ifield).Clfstd=std(Clf_arr);
Cldata(ifield).Cldavg=mean(Cld_arr);
Cldata(ifield).Cldstd=std(Cld_arr);
end
%% plot compare PS and map var of diff and full
for ifield=1:numel(fields)
    Clfavg=Cldata(ifield).Clfavg;
    Cldavg=Cldata(ifield).Cldavg;
    ravg=Cldavg./Clfavg;
    
    Clfstd=Cldata(ifield).Clfstd;
    Cldstd=Cldata(ifield).Cldstd;
    rstd=Cldstd./Clfstd;
    
    sigfull_arr=Cldata(ifield).sigfull_arr;
    sigdiff_arr=Cldata(ifield).sigdiff_arr;
    rsig=mean(sigdiff_arr)./mean(sigfull_arr);
    
    figure
    pltCl=semilogx(l,ravg,'o-');hold on
    pltdCl=semilogx(l,rstd,'o-');
    pltsig=refline(0,rsig^2);pltsig.Color='k';
    plt1=refline(0,1);plt1.Color='k';set(plt1,'LineStyle','--');
    legend([pltCl,pltdCl,pltsig],...
    {'<avg(Cl)\_diff>_{lab}/<avg(Cl)\_full>_{lab}',...
     '<std(Cl)\_diff>_{lab}/<std(Cl)\_full>_{lab}',...
     '<var(map)\_diff>_{lab}/<var(map)\_full>_{lab}'},...
     'location','northwest');
    legend boxoff
    xlim([1e2,2e5]);
    xlabel('$\ell$','interpreter','latex','fontsize',18);
    ylabel('diff/full','fontsize',18);
    title(sprintf('%s,N=%d',fields(ifield).name,fields(ifield).nfr));
    savename=strcat(sprintf('%sNoisePenalty/inst%di%d',savedir,inst,ifield));
    print(savename,'-dpng');%close    
end
%% see what's wrong w/ 27th mode 
ifield=8;
figure
for i=25:29
plot(Cld_arr(:,i)./Cldavg(i),'o-');hold on
end
title(sprintf('%s,N=%d, diff',fields(ifield).name,fields(ifield).nfr));
legend({'ell 25','ell 26','ell 27','ell 28','ell 29'});
ylabel('$C_\ell/<C_\ell>_{lab}$','interpreter','latex','fontsize',18);
xlabel('lab data #');
savename=strcat(sprintf('%sNoisePenalty/mode27diff',savedir));
print(savename,'-dpng');%close 

figure
for i=25:29
plot(Clf_arr(:,i)./Clfavg(i),'o-');hold on
end
title(sprintf('%s,N=%d, full',fields(ifield).name,fields(ifield).nfr));
legend({'ell 25','ell 26','ell 27','ell 28','ell 29'});
ylabel('$C_\ell/<C_\ell>_{lab}$','interpreter','latex','fontsize',18);
xlabel('lab data #');
savename=strcat(sprintf('%sNoisePenalty/mode27full',savedir));
print(savename,'-dpng');%close    
%%
for ifield=1:numel(fields)
Cl2ddiff_arr=Cldata(ifield).Cl2ddiff_arr;
weightd=1./(squeeze(std(Cl2ddiff_arr)));
Cldstd=Cldata(ifield).Cldstd;

rstd_arr=zeros(1,29);
for i=1:size(Cl2ddiff_arr,1)
for iell=9:29
mask=zeros(1024);mask((ell >= binl(iell)) & (ell <= binl(iell+1)))=1;

Cl2d=Cl2ddiff_arr(i,:,:);
Cl2d=Cl2d(find(mask));
w=weightd(find(mask));
w=w.*numel(w)./sum(w);
wCl=w.*Cl2d;
rstd_arr(iell)=std(wCl)./sqrt(numel(w))/Cldstd(iell);
end
loglog(l,rstd_arr,'ko-');hold on
drawnow
end
xlim([1e2,2e5]);
ylabel('$\sigma_{annuli}/\sigma_{lab}$','interpreter','latex','fontsize',18);
xlabel('$\ell$','interpreter','latex','fontsize',18);
plt1=refline(0,1);plt1.Color='k';set(plt1,'LineStyle','--');
title(sprintf('%s,N=%d',fields(ifield).name,fields(ifield).nfr));
savename=strcat(sprintf('%sNoisePenalty/pix_inst%di%d',...
                        savedir,inst,ifield));
print(savename,'-dpng');close    

end