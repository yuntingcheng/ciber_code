%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get read noise 1D PS diff=(1st-2nd)/sqrt(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

band=2;
beam=7;
iter_clip=3;
flight=40030;
cp=get_cal_params('flight',flight);
cal=cp(band).apf2eps.*cp(band).eps2nWpm2ps;

DCdir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/';
load(strcat(DCdir,'band',num2str(band),'_DCtemplate'),'DCtemplate');
alldatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(alldatdir,'band',num2str(band),'_alldat'),'alldat');
load(strcat(alldatdir,'band',num2str(band),'_FF5'),'FF5');

savedir=strcat('/Users/ytcheng/ciber/doc/20160905_NoiseModel/TM',...
                    num2str(band),'/');
weightdir=strcat('/Users/ytcheng/ciber/doc/20160905_NoiseModel/TM',...
                    num2str(band),'/full/PS2D/');

darkdir='/Users/ytcheng/ciber/doc/20150810_DarkCurrent/';
%% diff frames
for fieldnum=1:8
load(strcat(weightdir,'b',num2str(band),'_i',...
      num2str(fieldnum), '_std_noise'),'std_noise');
weight=(fftshift(fftshift(1./std_noise)))';

bigmask=alldat(fieldnum).bigmask;
fieldname=alldat(fieldnum).name;
nfr=alldat(fieldnum).nfr;
Mkk=alldat(fieldnum).wMkk;
scandir1=sprintf('%sTM%d/%s/first/data/',darkdir,band,fieldname);
scan1=dir(strcat(scandir1,'dark*'));
scandir2=sprintf('%sTM%d/%s/second/data/',darkdir,band,fieldname);
scan2=dir(strcat(scandir2,'dark*'));

wnCl_arr=zeros(numel(scan1),29);
nCl_arr=zeros(numel(scan1),29);
wnClF_arr=zeros(numel(scan1),29);
nClF_arr=zeros(numel(scan1),29);
wnClFclip_arr=zeros(numel(scan1),29);

for i=1:numel(scan1)
disp(sprintf('field=%d,count=%d/%d',fieldnum,i,numel(scan1)));

fname1=sprintf('%sTM%d/%s/first/data/%s',...
        darkdir,band,fieldname,scan1(i).name);
fname2=sprintf('%sTM%d/%s/second/data/%s',...
        darkdir,band,fieldname,scan2(i).name);

dark1=load(fname1);dark1=dark1.darkmap;
dark2=load(fname2);dark2=dark2.darkmap;
diff=(dark1-dark2)./sqrt(2);diff=diff.*cal;
[~,mask]=get_skymap(diff,bigmask,iter_clip);

if band==2;diff=dc_offset_remove(diff,mask);end

[wCl,~,~,l]=get_Cl(diff,mask,Mkk,beam,weight);
[Cl,~,Cl2d]=get_Cl(diff,mask,Mkk,beam,ones(1024));

nCl_arr(i,:)=Cl;
wnCl_arr(i,:)=wCl;

%%%%%%%%%% FF corrected %%%%%%%

diff=diff./FF5;
diff(find(diff~=diff))=0;diff(find(diff==inf|diff==-inf))=0;
[~,mask]=get_skymap(diff,bigmask,iter_clip);
[wClF,~,~,l]=get_Cl(diff,mask,Mkk,beam,weight);
[ClF,~,Cl2dF]=get_Cl(diff,mask,Mkk,beam,ones(1024));
[wClFclip,~,~,l]=get_Cl(diff,mask,Mkk,beam,weight,'clipCl2d',1);

nClF_arr(i,:)=ClF;
wnClF_arr(i,:)=wClF;
wnClFclip_arr(i,:)=wClFclip;

%%% make plot %%%%
subplot(2,2,1)
imageclip(Cl2d);
subplot(2,2,2)
imageclip(Cl2d.*weight');
subplot(2,2,3)
imageclip(Cl2dF);
subplot(2,2,4)
imageclip(Cl2dF.*weight');
suptitle(sprintf('%s %d',fieldname,i))
savename=sprintf('%sdiff/PS2Dplot/i%d_%s',...
                savedir,fieldnum,scan1(i).name(6:end-4));
print(savename,'-dpng');close
end
save1Ddir=strcat(savedir,'diff/PS1D/');
save(sprintf('%sb%d_i%d_nCl_arr',save1Ddir,band,fieldnum),'nCl_arr');
save(sprintf('%sb%d_i%d_wnCl_arr',save1Ddir,band,fieldnum),'wnCl_arr');
save(sprintf('%sb%d_i%d_nClF_arr',save1Ddir,band,fieldnum),'nClF_arr');
save(sprintf('%sb%d_i%d_wnClF_arr',save1Ddir,band,fieldnum),'wnClF_arr');
save(sprintf('%sb%d_i%d_wnClFclip_arr',save1Ddir,band,fieldnum),...
                                                         'wnClFclip_arr');
end
%% plot the results

for fieldnum=1:8
fieldname=alldat(fieldnum).name;
save1Ddir=strcat(savedir,'diff/PS1D/');
load(sprintf('%sb%d_i%d_nCl_arr',save1Ddir,band,fieldnum),'nCl_arr');
load(sprintf('%sb%d_i%d_wnCl_arr',save1Ddir,band,fieldnum),'wnCl_arr');
load(sprintf('%sb%d_i%d_nClF_arr',save1Ddir,band,fieldnum),'nClF_arr');
load(sprintf('%sb%d_i%d_wnClF_arr',save1Ddir,band,fieldnum),'wnClF_arr');
load(sprintf('%sb%d_i%d_wnClFclip_arr',save1Ddir,band,fieldnum),...
                                                         'wnClFclip_arr');

fig=figure;
yvalue=(l.*(l+1).*mean(nCl_arr)./2./pi);
ebarlow=(l.*(l+1).*(mean(nCl_arr)-prctile(nCl_arr,16))./2./pi);
ebarup=(l.*(l+1).*(prctile(nCl_arr,84)-mean(nCl_arr))./2./pi);
pltnCl=errorbar(l,yvalue,ebarlow,ebarup,'.','color',[0,0.8,1]...
                            ,'markersize',20);hold on
      
yvalue=(l.*(l+1).*mean(wnCl_arr)./2./pi);
ebarlow=(l.*(l+1).*(mean(wnCl_arr)-prctile(wnCl_arr,16))./2./pi);
ebarup=(l.*(l+1).*(prctile(wnCl_arr,84)-mean(wnCl_arr))./2./pi);
pltwnCl=errorbar(l,yvalue,ebarlow,ebarup,'.','color',[0,0,1]...
                            ,'markersize',20);hold on

yvalue=(l.*(l+1).*mean(wnClF_arr)./2./pi);
ebarlow=(l.*(l+1).*(mean(wnClF_arr)-prctile(wnClF_arr,16))./2./pi);
ebarup=(l.*(l+1).*(prctile(wnClF_arr,84)-mean(wnClF_arr))./2./pi);
pltwnClF=errorbar(l,yvalue,ebarlow,ebarup,'.','color',[0,0,0]...
                            ,'markersize',20);hold on

yvalue=(l.*(l+1).*mean(wnClFclip_arr)./2./pi);
ebarlow=(l.*(l+1).*(mean(wnClFclip_arr)-prctile(wnClFclip_arr,16))./2./pi);
ebarup=(l.*(l+1).*(prctile(wnClFclip_arr,84)-mean(wnClFclip_arr))./2./pi);
pltwnClFclip=errorbar(l,yvalue,ebarlow,ebarup,'.','color',[1,0,0]...
                            ,'markersize',20);hold on

xlim([1e2,2e5]);
if band==1;ylim([1e-1,1e4]);end
if band==2;ylim([1e-3,5e2]);end

ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log')
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
legend([pltnCl,pltwnCl,pltwnClF,pltwnClFclip],...
{'RN','RN-weighted','RN-weigted-FFcorr','RN-weigted-FFcorr-clip'},...
       'Location','southeast','FontSize',15);
legend boxoff
title(fieldname);
imname=sprintf('%sb%d_i%d',save1Ddir,band,fieldnum);
print(imname,'-dpng');close
end

