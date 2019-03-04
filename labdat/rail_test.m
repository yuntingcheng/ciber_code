%{
framedir=('/Users/ytcheng/ciber/data/40030/framedata/');
files = dir(strcat(dd,'TM1*'));

% see railtime40030.m for the frame number

fr1=23;
fr2=63;
nfr=fr2-fr1+1;
frames_arr=zeros(nfr,1024,1024);
for i=fr1:fr2
    load(strcat(dd,files(i).name));
    frames_arr(i-fr1+1,:,:)=arraymap;   
end

map1=linfit_map(frames_arr,'verbose',0);

fr1=83;
fr2=93;
nfr=fr2-fr1+1;
frames_arr=zeros(nfr,1024,1024);
for i=fr1:fr2
    load(strcat(dd,files(i).name));
    frames_arr(i-fr1+1,:,:)=arraymap;   
end

map2=linfit_map(frames_arr,'verbose',0);

fr1=98;
fr2=119;
nfr=fr2-fr1+1;
frames_arr=zeros(nfr,1024,1024);
for i=fr1:fr2
    load(strcat(dd,files(i).name));
    frames_arr(i-fr1+1,:,:)=arraymap;   
end

map3=linfit_map(frames_arr,'verbose',0);
%}
%% load field half data and mask
load('/Users/ytcheng/ciber/data/40030/slopedata/TM1_SWIRE_1st.mat');
load('/Users/ytcheng/ciber/data/40030/slopedata/TM1_SWIRE_2nd.mat');

load('/Users/ytcheng/ciber/doc/20160919_FlightDat/band1_alldat');
mkk=alldat(8).wMkk;
mask=alldat(8).bigmask;
clear alldat
%%
framedir=('/Users/ytcheng/ciber/data/40030/framedata/');
files = dir(strcat(framedir,'TM1*'));

fr1=23;
fr2=35;
nfr=fr2-fr1+1;
frames_arr=zeros(nfr,1024,1024);
for i=fr1:fr2
    load(strcat(framedir,files(i).name));
    frames_arr(i-fr1+1,:,:)=arraymap;   
end

rail1=linfit_map(frames_arr,'verbose',0);

fr1=36;
fr2=48;
nfr=fr2-fr1+1;
frames_arr=zeros(nfr,1024,1024);
for i=fr1:fr2
    load(strcat(framedir,files(i).name));
    frames_arr(i-fr1+1,:,:)=arraymap;   
end

rail2=linfit_map(frames_arr,'verbose',0);
%% load other darks
dd=('/Users/ytcheng/ciber/doc/20150810_DarkCurrent/40030/TM1/SWIRE/');
files1=dir(strcat(dd,'first/data/dark*'));
files2=dir(strcat(dd,'second/data/dark*'));
lab1=load(strcat(dd,'first/data/',files1(1).name));lab1=lab1.darkmap;
lab2=load(strcat(dd,'second/data/',files2(1).name));lab2=lab2.darkmap;
%%
dfield=(rawmap1-rawmap2)./sqrt(2);

drail=(rail1-rail2)./sqrt(2);
[~,maskr]=get_skymap(drail,mask,3);

dlab=(lab1-lab2)./sqrt(2);
[~,maskl]=get_skymap(dlab,mask,3);

[Clf,~,Cl2df,l]=get_Cl(dfield,mask,mkk,7,ones(1024));
[Clr,~,Cl2dr]=get_Cl(drail,maskr,mkk,7,ones(1024));
[Cll,~,Cl2dl]=get_Cl(dlab,maskl,mkk,7,ones(1024));
%%
ell = get_l(1024,1024,7,1);
midpoint=513;

figure
setwinsize(gcf,1000,300)
subplot(2,3,1)
imagesc(log10(abs(Cl2df)));
colorbar
caxis([-12,-9])
title(sprintf('SWIREdiff'));
subplot(2,3,4)
imagesc(log10(abs(Cl2df)));
xlim([midpoint-50,midpoint+50]); ylim([midpoint-50,midpoint+50]);
colorbar
caxis([-12,-9])
hold on
[c,h]=contour(ell,[500,1000,3000,5000,8000]);
clabel(c,h)

setwinsize(gcf,1000,300)
subplot(2,3,2)
imagesc(log10(abs(Cl2dr)));
colorbar
caxis([-12,-9])
title(sprintf('rail'));
subplot(2,3,5)
imagesc(log10(abs(Cl2dr)));
xlim([midpoint-50,midpoint+50]); ylim([midpoint-50,midpoint+50]);
colorbar
caxis([-12,-9])
hold on
[c,h]=contour(ell,[500,1000,3000,5000,8000]);
clabel(c,h)

setwinsize(gcf,1000,300)
subplot(2,3,3)
imagesc(log10(abs(Cl2dl)));
colorbar
caxis([-12,-9])
title(sprintf('dark'));
subplot(2,3,6)
imagesc(log10(abs(Cl2dl)));
xlim([midpoint-50,midpoint+50]); ylim([midpoint-50,midpoint+50]);
colorbar
caxis([-12,-9])
hold on
[c,h]=contour(ell,[500,1000,3000,5000,8000]);
clabel(c,h)
%%
figure
cp=get_cal_params('flight',40030);
cal=cp(1).apf2eps.*cp(1).eps2nWpm2ps;
pltf=loglog(l,l.*(l+1).*Clf./2./pi.*cal.*cal);hold on
pltr=loglog(l,l.*(l+1).*Clr./2./pi.*cal.*cal);
pltl=loglog(l,l.*(l+1).*Cll./2./pi.*cal.*cal);


xlim([2e2,2e5]);ylim([1e-2,5e3]);

xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
legend([pltf,pltr,pltl],{'SWIRE','rail','dark'},...
       'Location','southeast','FontSize',15);
legend boxoff
